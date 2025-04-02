!
!+ Various checks to be performed in monitoring and assimilation
!
MODULE mo_fg_checks
!
! Description:
!
!   Various checks to be performed in monitoring and assimilation:
!
!     identifier    routine         check
!
!     CHK_BLACKLIST check_black     check for blacklisted reports
!     CHK_FG        check_fg        simple first guess check
!     CHK_OBS_ERR   check_fg        nominal observation error
!     CHK_SUFF      check_suff      final check for sufficient data in report
!     CHK_RULE      check_rule      check for specific rules
!     CHK_OPERATOR  check_operator  validity of obs.operator applications
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
!  fixes for zero number of observations in a box
! V1_5         2009/05/25 Andreas Rhodin
!  new check for minimum surface temperature
! V1_7         2009/08/24 Andreas Rhodin
!  new option for artificial data generation: add 1 to obs-fg
! V1_8         2009/12/09 Harald Anlauf
!  check_fg: implement rejection of low-wind/null-wind observations
!            fix FG check for wind observations (now isotropic)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  DISMISS negative bending angles
! V1_13        2011/11/01 Andreas Rhodin
!  whitelist; new rules:fr_land,_sea,zenith_angle,lon,lat,bound on obsv.vars
! V1_15        2011/12/06 Andreas Rhodin
!  option to remove invalid observations in the assimilation step
! V1_16        2011/12/09 Andreas Rhodin
!  account for satid, phase (or fov) for RADiances and wind (SCATTerometer)
! V1_20        2012-06-18 Andreas Rhodin
!  changes for KENDA: use plev instead of olev; add varno codes
! V1_22        2013-02-13 Andreas Rhodin
!  changes for RADAR, LETKF; add comment for dismissed reports
! V1_23        2013-03-26 Andreas Rhodin
!  add entries for /rules/ (for cloud top analysis): VN_CTH, VN_TRH, VN_LM/H
! V1_26        2013/06/27 Andreas Rhodin
!  option to fix misbehaviour of check_no_obs
! V1_27        2013-11-08 Robin Faulwetter
!  Rewrite of radiance flagging
! V1_28        2014/02/26 Robin Faulwetter
!  Introduced tsurf in /RULES/, change usage for radiances
! V1_31        2014-08-21 Andreas Rhodin
!  adapt evaluation of /RULES/ to RADAR observation operator
! V1_37        2014-12-23 Robin Faulwetter
!  Option for rules,that shall only be applied, if an instrument is missing
!  activate check for obs_err >> fg_err in case of fg_err==0
! V1_42        2015-06-08 Andreas Rhodin
!  minor cleanup
! V1_44        2015-09-30 Andreas Rhodin
!  check_rules: only pass instrument information if present (not for COSMO)
! V1_45        2015-12-15 Andreas Rhodin
!  adaptions to MEC and LETKF
! V1_46        2016-02-05 Harald Anlauf
!  check_rule: call get_rule for surface obs. passing station height
! V1_47        2016-06-06 Andreas Rhodin
!  add obstype as selection criterium in namelist /rules/.
!  namelist /RULES/: use rules applied to U,V,FF also for DD.
!  new subroutine check_redo: restore report status checks from 'flags'.
! V1_48        2016-10-06 Andreas Rhodin
!  handle T2M,TD2M,U10M,V10M,DD,FF,PS
!  fix for radiances in COSMO LETKF (Valerio Cardinaly)
! V1_49        2016-10-25 Harald Anlauf
!  bugfix for FF/DD monitoring
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured/unified cloud detection for radiances.
!  changes to run 3dvar with COSMO/COMET data              (A.Rhodin)
!  check_fg: don't check DD for 'too large observed value' (H.Anlauf)
! V1_51        2017-02-24 Andreas Rhodin
!  handle VN_PRH VN_PRH2M VN_Q2M VN_LWC (for COMET)
!              2018-06-14 Christoph Schraff
!  Estimation (in subr. estimate_fgchk) and use of an additional, systematic
!  error for fg check which depends on the observed atmospheric stability.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004  original source
!------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: wp, sp           ! precision kind parameters
  use mo_exception,  only: finish           ! abort on error conditions
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_bcast,        &! generic MPI_BCAST routine
                           p_sum            ! generic MPI sum routine
! use mo_run_params, only: method           ! 'PSAS', 'LETKF' etc
  use mo_dec_matrix, only: t_vector,       &! decomposed vector data type
                           assignment(=)    ! t_vector = real
  use mo_t_table,    only: name_value       ! find name of table entry
  use mo_fdbk_tables,only: VN_U,           &! u wind component     code
                           VN_V,           &! v wind component     code
                           VN_W,           &! w wind component     code
                           VN_U10M,        &! 10m u wind component code
                           VN_V10M,        &! 10m v wind component code
                           VN_FF,          &! wind speed           code
                           VN_DD,          &! wind direction       code
                           VN_HLOS,        &! horizontal line of sight wind
                           VN_PRH,         &! pseudo relative humidity code
                           VN_RH,          &! relative humidity    code
                           VN_PRH,         &! pseudo relative hum  code
                           VN_PRH2M,       &! 2m pseudo relative hum code
                           VN_RH2M,        &! 2m relative humidity code
                           VN_TD,          &! dewpoint temperature code
                           VN_TD2M,        &! 2m dewpoint temper.  code
                           VN_T,           &! temperature          code
                           VN_T2M,         &! 2m temperature       code
                           VN_VT,          &! virtual temperature  code
                           VN_Q,           &! specific humidity    code
                           VN_Q2M,         &! 2m specific humidity code
                           VN_Z,           &! geopotential height  code
                           VN_HEIGHT,      &! height               code
                           VN_P,           &! pressure             code
                           VN_PS,          &! surface pressure     code
                           VN_BENDANG,     &! bending angle        code
                           VN_RADVEL,      &! radial velocity      code
                           VN_RREFL,       &! radar reflectivity   code
                           VN_REFL,        &!       reflectivity   code
                           VN_RAD_GL,      &! global radiation     code
                           VN_CTH,         &! cloud top height     code
                           VN_TRH,         &! transformed humidity code
                           VN_N_L,         &! low  cloud fraction  code
                           VN_N_M,         &! med. cloud fraction  code
                           VN_N_H,         &! high cloud fraction  code
                           VN_NH,          &! cloud base height    code
                           VN_PWSOL,       &! solar power data     code
                           VN_PWIND,       &! wind  power data     code
                           VN_LWC,         &! liquid water content code
                           VN_NSOILM,      &! normalized soil moisture code
                           VN_SOILM,       &! volumetric soil moisture code
                           VN_SPD,         &! slant path delay     code
                           VN_ZPD,         &! zenith path delay    code
                           VN_FLEV,        &! nominal flight level code
                           VN_HOSAG,       &! height above ground  code
                           VN_MIXR,        &! mixing ratio
                           OT_SYNOP,       &! SYNOP    report type
                           OT_TEMP,        &! TEMP     report type
                           OT_DRIBU,       &! BUOY     report type
                           OT_SATEM,       &! SATEM    report type
                           OT_RAD,         &! ATOVS    report type
                           OT_SATOB,       &! SATOB    report type
!                          OT_AIREP,       &! Aircraft report type
                           OT_GPSRO,       &! GPS radio occultation report type
                           OT_GPSGB,       &! GPS ground based (ZTD, STD)
                           OT_RADAR,       &! RADAR    report type
                           OT_SCATT,       &! SCATTerometer report type
                           OT_PILOT,       &! PILOT    report type
                           obstype,        &! list of observation types
                           n_ot,           &! number of observation types
                           SUR_SEA,        &! surftype: sea
                           SUR_LAND,       &! surftype: land
                           SUR_HIGHLAND,   &! surftype: highland
                           SUR_MISSING,    &! surftype: missing
                           SUR_MISMATCH,   &! surftype: mismatch
!                          SUR_MWSURF,     &! surftype: bad stype detected by MW sounder
                           SUR_BLK_PP,     &! surftype: blacklisted by sat_pp
                           TF_SURF_TYPE,   &! tovs_flag: wrong surface type
                           TF_SURF_INFL,   &! tovs_flag: surface influence
                           TF_SURF_MODEL,  &! tovs_flag: surface properties (model)
                           TF_SURF_RETR,   &! tovs_flag: surface type mismatch (retrieved and model)
                           TF_EMIS_FAILED   ! tovs_flag: emiss./refl. calculation failed
  use mo_t_obs,      only: t_obs,          &! generic observation data type
                           t_spot,         &! report meta data type
                           t_head,         &! report meta data type
                           fix_no_obs,     &! .false. to revert no_obs-check
                           oq_name,        &! convert observed to mnemonic
                           varno_oq,       &! convert varno to observed quantity
                           COSMO            ! observation module numbers
  use mo_obs_set,    only: t_obs_set        ! observation data type set
  use mo_temp,       only: fgchk_inv_t,    &! stability-dep fg-check limit (T)
                           fgchk_inv_q      ! stability-dep fg-check limit (RH)
  use mo_amv,        only: check_amv_surf   ! check AMVs for sea/land surface
  use mo_rad,        only: USE_LAND,       &! Use radiances over land
                           USE_SEA,        &! use radiances over sea
                           USE_SEAICE,     &! use radiances over seaice
                           USE_MINTSURF,   &! use radiances over areas with tsurf<min_tsurf
                           USE_HIGHLAND,   &! use radiances over high land
                           USE_MISMATCH,   &!
                           USE_BLK_PP,     &!
                           USE_QCNOTUSE,   &! do noe use for QC
                           m_chan,         &!
                           t_rad_set,      &
                           n_set
  use mo_tovs,       only: decr_tovs_use    ! flagging routine for tovs data
  use mo_tovs_prof,  only: set_surface_type,&! set/check surface type for radiances
                           check_recalc_stype! check recalculated surface type (RAD)
  use mo_t_tovs,     only: t_tovs,         &! TOVS specific information
                           t_tovs_instr,   &! info on instruments in t_tovs
                           load,           &! load t_tovs
                           destruct,       &! destruct t_tovs
                           mx_nlev,        &! max. number of
                           TTOVS_CI,       &! flag: load t_tovs%ci
                           TTOVS_FLAG,     &! flag: load t_tovs%flag
                           TTOVS_EMIS       ! flag: load t_tovs%emis
  use mo_t_datum,    only: t_datum,        &! body entry derived data type
                           rvind            ! invalid value
  use mo_t_use,      only: decr_use,       &! decrease state of report or datum
                           CHK_FG,         &! first guess       check id
                           CHK_OBS_ERR,    &! observation error check id
                           CHK_SURF,       &! correct surface   check id
                           CHK_NO_OBS,     &! insufficient data check id
                           CHK_NONE,       &! valid data        check id
                           CHK_RULE,       &! specific rule     check id
                           CHK_DOMAIN,     &! valid area        check id
                           CHK_BLACKLIST,  &! blacklist         check id
                           CHK_WHITELIST,  &! missing whitelist check id
                           CHK_OPERATOR,   &! operator validity check id
                           CHK_FINAL,      &! final             check id
!                          CHK_CLOUD,      &! cloud             check id
                           CHK_INSDAT,     &! insufficient data check id
                           CHK_NOTUSED,    &! data not used     check id
                           STAT_ABORT,     &! abort
                           STAT_DISMISS,   &! dismiss
!                          STAT_FORGET,    &! forget
                           STAT_OBS_ONLY,  &! observation only
                           STAT_REJECTED,  &! rejected
                           STAT_PASSIVE,   &! passive
                           STAT_NOTACTIVE, &! passive or rejected
                           STAT_ACTIVE_0I, &! (deprecated)
                           STAT_ACTIVE,    &! active
!                          STAT_ACTIVE_1I, &! (deprecated)
                           STAT_ACTIVE_1,  &!
                           STAT_ACCEPTED,  &! accepted
!                          stat_mnem,      &!
                           chk,            &! mnemonics ... for 'checks'
                           n_chk            ! max. number of check
  use mo_obs_tables, only: decr_rpt_use,   &! change use-flags of report
                           t_rept_use,     &! use table entry data type
                           rept_use,       &! use table entry
                           obstyp,         &! table of obs. operator types
                           write_report,   &! write report to status file
                           write_pending    ! write pending reports to file
  use mo_obs_rules,  only: get_rule,       &! get observation processing rules
                           iu1,            &! undefined integer value
                           empty_set,      &! empty data type t_set
                           t_set,          &! data type
                           t_ilev,         &! Instrument "level" information
                           PRC_ADDFG,      &! flag value: add first guess
                           PRC_ADD_1,      &!             add 1 to observation
                           PRC_ADD_R,      &!             add obs.err. to obs.
                           PRC_SUB_R,      &!             subtract obs.err.
                           PRC_ADDOBSER,   &!             add random obs.error
                           PRC_ADDOBSE2,   &!             add random obs.error
                           PRC_OBS_FG       !             set obs to fg
  use mo_algorithms, only: random_gauss     ! normal distribution
  use mo_random,     only: random_state_t   ! Derived type:random generator state
  use mo_blacklist,  only: t_black,        &! blacklist entry data type
                           black_entry,    &! find entry in blacklist
                           verbose          ! verbosity level of blacklisting
  use mo_p_output,   only: oline,          &! output line buffer
                           iol,            &! number of next line to write
                           nextline,       &! routine to increment line number
                           flush_buf,      &! routine to write buffer
                           add_line_pio     ! routine to write string on I/O PE
  use mo_physics,    only: t0c,            &! 273.15 K
                           gacc             ! gravity acceleration
  use mo_physical_constants, &
                     only: rd               ! gas constant (287.04) [J/K/kg]

  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: check_rule     ! check for specific rules
  public :: check_black    ! check for blacklisting
  public :: check_fg       ! perform generic first guess check
  public :: check_gross    ! generic gross check (to be called from MEC)
  public :: check_eo       ! modify observational error
  public :: check_suff     ! final check for sufficient data in report
  public :: check_obs      ! check for valid report
  public :: check_cons     ! consistency check
  public :: check_nofgscan ! to be called if the first guess scan is skipped
  public :: check_operator ! check applicability of observation operator
  public :: check_accept   ! final check for acceptance
  public :: check_redo     ! redo report check based on 'flags'
  public :: check_obstype  ! check that obstype is in list (for MEC)
  public :: check_var      ! reject observations not handled by var scheme
  public :: estimate_fgchk ! estimate additional error for fg check

contains
!------------------------------------------------------------------------------
  subroutine check_rule (obs)
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data
  !-------------------------
  ! check for specific rules
  !-------------------------

    integer                        :: ib, is, i ! indices
    integer                        :: iinstr    ! RAD instrument index
    type(t_spot)      ,pointer     :: si        ! pointer to report meta data
    type(t_head)      ,pointer     :: hd        ! pointer to report meta data
    type(t_datum)     ,pointer     :: d
    integer                        :: u         ! use flag
    integer                        :: u_aux     ! use flag
    type(t_rept_use)  ,pointer     :: use       ! use table entry
    integer                        :: state
    integer                        :: nc
    type(t_set)                    :: o
    real(wp)                       :: pobs
    real(wp)                       :: zobs     ! observation height above surface
    type(t_tovs)                   :: tovs
    type(t_tovs_instr)             :: ti
    type(t_ilev)      ,allocatable :: levs(:)
    type(t_rad_set)   ,pointer     :: rs => null()
    real(wp)                       :: z_bg     ! model orography
    logical                        :: lchk
    logical                        :: l_rs_flag

    !----------------
    ! loop over boxes
    !----------------
    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      !------------------------
      ! loop over reports, FOVs
      !------------------------
      do is = 1, obs(ib)% n_spot
        si  => obs(ib)% spot(is)
        hd  => si% hd
        !-----------------
        ! handle RADiances
        !-----------------
        if (hd% obstype == OT_RAD) then
          u_aux = STAT_ACCEPTED
!         z_bg  = si%gp_bg / gacc
          l_rs_flag = .false.
          if (si%p%n > 0 .and. n_set > 0) then
            call load(obs(ib), si, tovs, rs=rs, ti=ti, tovs_io=TTOVS_CI+TTOVS_FLAG+TTOVS_EMIS)
            l_rs_flag = associated(rs%flag)
            nc = tovs%nchan
            if (.not.allocated(levs)) allocate(levs(m_chan))
            levs(1:nc)%value = rs%chan(tovs%ci(1:nc))
            do i = 1, ti%n_instr
              levs(ti%o_ch_i(i)+1:ti%o_ch_i(i)+ti%n_ch_i(i))%instr = rs%instr(ti%ii(i))
            end do
            if (l_rs_flag) then
              nc = nc  - count(btest(rs%flag(tovs%ci(1:tovs%nchan)), USE_QCNOTUSE))
              levs(1:nc) = pack(levs(1:tovs%nchan), mask=.not.btest(rs%flag(tovs%ci(1:tovs%nchan)), USE_QCNOTUSE))
            end if
!Comment: To restore old behaviour, uncomment lines starting with !res
!         and comment lines between !removemaker1 and !removemarker2
!res          else
!res            call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz, &
!res                           si% statid,                                       &
!res                           obstype=hd%obstype, codetype=hd%codetype,         &
!res                           psurf=si%ps_bg, satid= int (hd% satid),           &
!res                           phase= int (si% phase), zsurf=z_bg,               &
!res                           lat= si% col% c% dlat, lon= si% col% c% dlon,     &
!res                           sol_zenith= si% sozen, tsurf=si%ts_bg,o=o,        &
!res                           use = u_aux)
          endif
        endif
        !------------------------------
        ! loop over single observations
        !------------------------------
        z_bg = si%gp_bg / gacc
        do i = si%o% i + 1, si%o% i + si%o% n
          d => obs(ib)% body(i)
          u = STAT_ACCEPTED
          !--------------
          ! general rules
          !--------------
          select case (hd% obstype)
          case (OT_RAD)
            if (si%p%n > 0 .and. n_set > 0) then
              iinstr = 1
              if (ti%n_instr > iinstr) then
                if (i-si%o%i > ti%o_ch_i(iinstr+1)) iinstr = iinstr + 1
              end if
              !--------------------------
              ! RADiance observation type
              !--------------------------
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz, &
                             si% statid,                                       &
                             obstype=hd%obstype, codetype=hd%codetype,         &
                             psurf=si%ps_bg,                                   &
                             satid= int (hd% satid), phase= int (si% phase),   &
                             instr=rs%instr(ti%ii(iinstr)), zsurf=z_bg,        &
                             instrs=rs%instr(ti%ii(1:ti%n_instr)),             &
                             levs=levs(1:nc),                                  &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,     &
                             zobs= obs(ib)% olev (i), sol_zenith= si% sozen,   &
                             sat_zenith=si% stzen, fr_land= si% sl_bg,         &
                             tsurf=si%ts_bg, o=d% set, use = u)
            else
!removemarker1 do not remove this line - unless above !Comment is depricated
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz, &
                             si% statid,                                       &
                             obstype=hd%obstype, codetype=hd%codetype,         &
                             psurf=si%ps_bg,                                   &
                             satid= int (hd% satid), phase= int (si% phase),   &
                             instr=si%sttyp, zsurf=z_bg,                       &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,     &
                             zobs= obs(ib)% olev (i),                          &
                             sol_zenith= si% sozen, fr_land= si% sl_bg,        &
                             tsurf=si%ts_bg,o=o, use = u_aux)
!removemarker2 do not remove this line - unless above !Comment is depricated
              u = u_aux
              d% set = o
            endif
          case (OT_GPSRO)
            !-----------------------------------
            ! radio occultation observation type
            !-----------------------------------
            call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                           si% statid, obstype=hd%obstype,                  &
                           codetype=hd%codetype,  center= si% center_id,    &
                           satid= int (hd% satid), phase= int (si% phase),  &
                           lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                           fr_land= si% sl_bg, zobs= obs(ib)% olev (i),     &
                           sol_zenith= si% sozen, tsurf=si%ts_bg, o= d% set,&
                           use = u)
          case (OT_GPSGB)
            !-----------------------------------
            ! GNSS ground based observation type
            !-----------------------------------
            call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                           si% statid, obstype=hd%obstype,                  &
                           codetype=hd%codetype,                            &
                           lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                           fr_land= si% sl_bg, zobs= obs(ib)% olev (i),     &
                           sol_zenith= si% sozen, tsurf=si%ts_bg, o= d% set,&
                           use = u)
          case (OT_RADAR)
            !-----------------------
            ! RADAR observation type
            !-----------------------
            select case (obs(ib)% body(i)% lev_typ)
            case (VN_HEIGHT)
              zobs = obs(ib)% olev (i) - si% z
            case (VN_HOSAG)
              zobs = obs(ib)% olev (i)
            case default
              zobs = -999._wp
            end select
            select case (obs(ib)% varno(i))
            case (VN_RADVEL)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, zobs=zobs,                           &
                             psurf=si%ps_bg, uv=d% set, obstype=hd%obstype,   &
                             codetype=hd%codetype,                            &
                             satid= int (hd% satid), phase= int (si% phase),  &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith= si% sozen, sat_zenith=si% stzen,     &
                             tsurf=si%ts_bg, use = u)
            case (VN_RREFL, VN_REFL)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, zobs=zobs,                           &
                             psurf=si%ps_bg, o=d% set, obstype=hd%obstype,    &
                             codetype=hd%codetype,                            &
                             satid= int (hd% satid), phase= int (si% phase),  &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith= si% sozen, sat_zenith=si% stzen,     &
                             tsurf=si%ts_bg, use = u)
            end select
          case (OT_PILOT)
            !-----------------------
            ! PILOT observation type
            !-----------------------
            pobs = obs(ib)% body(i)% plev
            select case (obs(ib)% body(i)% lev_typ)
            case (VN_HEIGHT)
              ! Wind profilers
              zobs = obs(ib)% olev (i) - si% z
            case (VN_HOSAG)
              zobs = obs(ib)% olev (i)
            case default
              ! PILOT reports (using or converted to pressure levels)
              ! do not have a supported height above ground
              zobs = -999._wp
            end select
            select case (obs(ib)% varno(i))
            case (VN_U, VN_V, VN_FF, VN_DD)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, uv=d% set,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             zobs=zobs, zsurf=z_bg, fr_land=si% sl_bg,        &
                             sol_zenith= si% sozen, tsurf=si%ts_bg, use = u)
              if (abs(si% col% c% dlat) > 89._wp) then
                 state = STAT_DISMISS
                 call decr_use (d% use, state, CHK_DOMAIN)
              endif
            case(VN_T, VN_VT)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, t=d% set, &
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             zobs=zobs, zsurf=z_bg, fr_land=si% sl_bg,        &
                             sol_zenith= si% sozen, tsurf=si%ts_bg, use = u)
            case(VN_RH, VN_Q, VN_MIXR)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, q=d% set, &
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             zobs=zobs, zsurf=z_bg, fr_land=si% sl_bg,        &
                             sol_zenith= si% sozen, tsurf=si%ts_bg, use = u)
            case (VN_W)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, o=d% set, &
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             varno=obs(ib)% varno(i),                         &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             zobs=zobs, zsurf=z_bg, fr_land=si% sl_bg,        &
                             sol_zenith= si% sozen, tsurf=si%ts_bg, use = u)
            case (VN_RADVEL)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, uv=d% set,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             zobs=zobs, zsurf=z_bg, fr_land=si% sl_bg,        &
                             sol_zenith= si% sozen, tsurf=si%ts_bg, use = u)
            case default
              !-----------------------------------
              ! other variables (should not occur)
              !-----------------------------------
              write(0,*) "check_rule(PILOT), not implemented: ", &
                   si%statid,hd%codetype,obs(ib)% varno(i)
              cycle
            end select
          !------------------------
          ! other observation types
          !------------------------
          case default
            pobs = obs(ib)% body(i)% plev
            select case (obs(ib)% varno(i))
            case (VN_T, VN_T2M)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, t=d% set, &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
            case (VN_Z, VN_HEIGHT, VN_CTH)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, gp=d% set,&
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
            case (VN_RH, VN_RH2M, VN_TRH, VN_TD, VN_TD2M, VN_PRH, VN_PRH2M)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, q=d% set, &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
            case (VN_U, VN_V, VN_U10M, VN_V10M, VN_HLOS, VN_FF, VN_DD, VN_PWIND)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, uv=d% set,&
                             satid= int (hd% satid), phase= int (si% phase),  &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
               if (abs(si% col% c% dlat) > 89._wp) then
                 state = rept_use (hd% obstype)% use (CHK_DOMAIN)
                 state = STAT_DISMISS
                 call decr_use (d% use, state, CHK_DOMAIN)
               endif
            case (VN_PS)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg,  p=d% set,&
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
            case (VN_Q, VN_Q2M)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, q=d% set, &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
!           case (OBS_DUM)
!             call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
!                            si% statid, pobs=pobs, psurf=si%ps_bg, o=d% set, &
!                            lat= si% col% c% dlat, lon= si% col% c% dlon,    &
!                            sol_zenith= si% sozen,                           &
!                            obstype=hd%obstype, codetype=hd%codetype, use = u)
            case (VN_N_L, VN_N_M, VN_N_H, VN_PWSOL, VN_RAD_GL, VN_NH)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg, o=d% set, &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             varno=obs(ib)% varno(i),                         &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
            case (VN_LWC, VN_SOILM, VN_NSOILM)
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, o=d% set,                 &
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, zobs=si% z,                &
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             varno=obs(ib)% varno(i),                         &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
            case default
              !----------------
              ! other variables
              !----------------
              call get_rule (hd%modtype, hd%buf_type, hd%buf_subtype, hd%dbkz,&
                             si% statid, pobs=pobs, psurf=si%ps_bg,  o=d% set,&
                             lat= si% col% c% dlat, lon= si% col% c% dlon,    &
                             sol_zenith=si% sozen, tsurf=si%ts_bg, zobs=si% z,&
                             obstype=hd%obstype, codetype=hd%codetype,        &
                             fr_land=si% sl_bg, instr=si% sttyp, use = u)
              cycle
            end select
          end select
          u = min (u, d% set% use)
          call decr_use (d% use, u, CHK_RULE)
          !------------------------------------
          ! abort if unreasonable state was set
          !------------------------------------
          if (d% use% state <= STAT_ABORT) then
            write(0,*) 'check_rule: u, set, state, typ, vn=',          &
              u, d% set% use,d% use% state,hd% obstype,obs(ib)% varno(i)
            call finish('check_rule','state <= STAT_ABORT')
          endif
        end do

        select case (hd% obstype)
        !---------------
        ! specific rules
        !---------------
        !------
        ! SATOB
        !------
        case (OT_SATOB)
          call check_amv_surf (si)
        end select
        !--------------------
        ! other general rules
        !--------------------
        use => rept_use (hd% obstype)
        if (        si% sl_bg < use% fr_land  ) &
             call decr_use_local(comment="fr_land", not_bit=USE_SEA     , tovs_flag=TF_SURF_MODEL)
        if (1._wp - si% sl_bg < use% fr_sea   ) &
             call decr_use_local(comment="fr_sea" , not_bit=USE_LAND    , tovs_flag=TF_SURF_MODEL)
        if (1._wp - si% fi_bg < use% fr_noice ) &
             call decr_use_local(comment="fr_ice" , not_bit=USE_SEAICE  , tovs_flag=TF_SURF_MODEL)
        if (-t0c  + si% ts_bg < use% min_tsurf) &
             call decr_use_local(comment="tsurf"  , not_bit=USE_MINTSURF, tovs_flag=TF_SURF_MODEL)

        if (hd% obstype == OT_RAD .and. l_rs_flag) then
          if (btest(si% stlsf, SUR_SEA).and..not.btest(si% stlsf, SUR_LAND)) then
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 not_bit=USE_SEA, flag_rpt=.false., hint_debug='sur_sea')
          end if
          if (btest(si% stlsf, SUR_LAND)) then
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 not_bit=USE_LAND, flag_rpt=.false., hint_debug='sur_land')
          end if
          if (btest(si% stlsf, SUR_MISSING)) then
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 not_bit=USE_SEA, flag_rpt=.false., hint_debug='sur_missing')
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 not_bit=USE_LAND, flag_rpt=.false., hint_debug='sur_missing')
          end if
          if (btest(si% stlsf, SUR_HIGHLAND)) then
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 not_bit=USE_HIGHLAND, flag_rpt=.false., hint_debug='sur_highland')
          end if
          if (btest(si% stlsf, SUR_MISMATCH)) then
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 not_bit=USE_MISMATCH, flag_rpt=.true., hint_debug='sur_mismatch')
          end if
          if (btest(si% stlsf, SUR_BLK_PP)) then
            call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_TYPE, &
                 bit=USE_BLK_PP, flag_rpt=.false., hint_debug='sur_blk_pp')
          end if
          ! QC for MW / IR emissivity or snow in BRDF climatology
          if (associated(tovs%flag)) then
            if (any(btest(tovs%flag(1:tovs%nchan), TF_EMIS_FAILED))) then
              call decr_tovs_use(si, obs(ib),CHK_INSDAT, set=rs, tovs=tovs, state=STAT_REJECTED,&
                   chan_mask=btest(tovs%flag(1:tovs%nchan), TF_EMIS_FAILED),                    &
                   hint_debug='emis_failed')
            end if
          end if
          ! Check whether tovs%mw_stype makes a difference in surface type classification
          ! compared to original surface type classification (where tovs%mw_stype was not known)
          if (tovs%mw_stype >= 0 .and. check_recalc_stype) then
            lchk = .true.
            call set_surface_type(si, tovs, l_check=lchk)
            if (.not.lchk) then
              call decr_tovs_use(si, obs(ib), CHK_SURF, tovs=tovs, tovs_flag=TF_SURF_RETR, &
                   not_bit=USE_MISMATCH, flag_rpt=.false., hint_debug='recalc_stype')
            end if
          end if
          call destruct(tovs)
        end if

      end do
    end do

  contains

    subroutine decr_use_local(comment, not_bit, tovs_flag)
    !-------------------------------------------------------------------------
    ! for RADIANCES:
    ! more sophisticated flagging taking into account channel characteristics.
    ! Not for 'LETKF' (=KENDA) as channel characteristics are currently
    ! not available.
    !-------------------------------------------------------------------------
      character(len=*), intent(in)           :: comment
      integer,          intent(in), optional :: not_bit
      integer,          intent(in), optional :: tovs_flag
      if (si% hd% obstype == OT_RAD .and. l_rs_flag) then
        call decr_tovs_use(si, obs(ib), CHK_SURF, not_bit=not_bit, &
             hint_debug='check_rule '//trim(comment), tovs_flag=tovs_flag)
      else
        call decr_rpt_use (si, CHK_SURF, comment=comment)
      end if

    end subroutine decr_use_local

   end subroutine check_rule
!------------------------------------------------------------------------------
  subroutine check_black (obs)
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data
  !--------------------------------------------------
  ! Check for blacklisted reports (observation level)
  !-------------------------------------------------
    integer ,parameter         :: mb = 3          ! max.# of blacklist entries
    integer                    :: ib, is, i, n, j ! indices
    type(t_spot)      ,pointer :: si              ! pointer to report meta data
    type(t_head)      ,pointer :: hd              ! pointer to report meta data
    type(t_black)              :: black (mb)      ! blacklist entry
    type(t_datum)     ,pointer :: d
    integer                    :: stat_new(n_ot)  ! status after blacklisting
    real(wp)                   :: pobs            ! "pressure" of observation
    real(wp)                   :: zobs            ! "height"   of observation
    logical                    :: lpobs
    real(wp)                   :: scalf           ! scaling factor for biascor.
    !-------------
    ! write report
    !-------------
    call flush_buf
!   call add_line_pio (repeat('-',79))
    call add_line_pio ('')
    call add_line_pio ('  Processing BLACKLIST')
    call add_line_pio ('')
    !-------------------------------------
    ! New report status after blacklisting
    !-------------------------------------
    stat_new(1:n_ot) = rept_use (1:n_ot)% use (CHK_BLACKLIST)
    !------------------
    ! loop over reports
    !------------------
    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      do is = 1, obs(ib)% n_spot
        si  => obs(ib)% spot(is)
        hd  => si% hd
        !---------------------------------------------------------------
        ! skip blacklist check for obstypes processed by COSMO operators
        !---------------------------------------------------------------
        if (hd% modtype == COSMO) cycle
        !------------------------------------
        ! conventional observation types only
        !------------------------------------
        select case (hd% obstype)
        case (OT_SYNOP:OT_SATEM,OT_GPSGB)
          n = black_entry (si% statid, hd% obstype, hd% codetype, black)
          if (n < 0) call decr_rpt_use (si, CHK_BLACKLIST)
          if (n > 0) then
            if (verbose > 0) then
              !-------------
              ! write report
              !-------------
              do j = 1, n
                call nextline
                write(oline(iol),'(a,1x,a,1x,a,1x,a,1x,l1)') 'blacklisted  ', &
                  obstyp(hd% obstype)% name, si% statid, black(j)% statid,    &
                  black(j)% mwl
                if (verbose > 1) then
                  if (black(j)% lg) then
                    call nextline
                    write(oline(iol),'(a,2i8,1x,2f8.2)') &
                    '  height     ', black(j)% pg(2), black(j)% pg(1), black(j)% zg(2), black(j)% zg(1)
                  endif
                  if (black(j)% lt) then
                    call nextline
                    write(oline(iol),'(a,2i8,1x,2f8.2)') &
                    '  temperature', black(j)% pt(2), black(j)% pt(1), black(j)% zt(2), black(j)% zt(1)
                  endif
                  if (black(j)% ld) then
                    call nextline
                    write(oline(iol),'(a,2i8,1x,2f8.2)') &
                    '  humidity   ', black(j)% pd(2), black(j)% pd(1), black(j)% zd(2), black(j)% zd(1)
                  endif
                  if (black(j)% lw) then
                    call nextline
                    write(oline(iol),'(a,2i8,1x,2f8.2)') &
                    '  wind       ', black(j)% pw(2), black(j)% pw(1), black(j)% zw(2), black(j)% zw(1)
                  endif
                endif
              end do
            end if
            !----------------------------------
            ! check for missing whitelist entry
            !----------------------------------
            if (black(1)% mwl) call decr_rpt_use (si, CHK_WHITELIST)
            !--------------------------------------
            ! check individual observations, levels
            !--------------------------------------
            do j = 1, n
              do i = si%o% i + 1, si%o% i + si%o% n
                d => obs(ib)% body(i)
                if (black(j)% mwl) call decr_use (d% use, check= CHK_WHITELIST)
                !---------------------------------
                ! Use observed or derived pressure
                !---------------------------------
                lpobs = .true.
                select case (obs(ib)% body(i)% lev_typ)
                case default
                   pobs = obs(ib)% olev(i)
                case (VN_FLEV, VN_HEIGHT)
                   zobs = obs(ib)% olev(i)
                   lpobs = .false.
                case(VN_HOSAG)
                   zobs = obs(ib)% olev(i) + si% z
                   lpobs = .false.
                end select
                select case (obs(ib)% varno(i))
                case (VN_Z, VN_PS) !, OBS_HS)
                  if (black(j)% bcg /= 0.) then
                    scalf = merge (si% pz_bg, 1._wp, obs(ib)% varno(i) == VN_PS)
                    d% o  = d% o - d% bc
                    d% bc = black(j)% bcg * scalf
                    d% o  = d% o + d% bc
                    if (verbose > 0) then
                      call nextline
                      write(oline(iol),'(a,a,f8.2)') &
                      '  height     ', 'biascor:',black(j)% bcg
                    end if
                  endif
                  if (.not. black(j)% lg)     cycle
                  if (lpobs) then
                    if (black(j)% pg(1) < pobs .or. &
                        black(j)% pg(2) > pobs) cycle
                  else
                    if (black(j)% zg(1) > zobs .or. &
                        black(j)% zg(2) < zobs) cycle
                  end if
                case (VN_T, VN_T2M)
                  if (.not. black(j)% lt)     cycle
                  if (lpobs) then
                    if (black(j)% pt(1) < pobs .or. &
                        black(j)% pt(2) > pobs) cycle
                  else
                    if (black(j)% zt(1) > zobs .or. &
                        black(j)% zt(2) < zobs) cycle
                  end if
                case (VN_RH, VN_RH2M, VN_TD, VN_TD2M, VN_Q, VN_Q2M)
                  if (.not. black(j)% ld)     cycle
                  if (lpobs) then
                    if (black(j)% pd(1) < pobs .or. &
                        black(j)% pd(2) > pobs) cycle
                  else
                    if (black(j)% zd(1) > zobs .or. &
                        black(j)% zd(2) < zobs) cycle
                  end if
                case (VN_U, VN_V, VN_U10M, VN_V10M, VN_FF, VN_DD)
                  if (.not.black(j)% lw)      cycle
                  if (lpobs) then
                    if (black(j)% pw(1) < pobs .or. &
                        black(j)% pw(2) > pobs) cycle
                  else
                    if (black(j)% zw(1) > zobs .or. &
                        black(j)% zw(2) < zobs) cycle
                  end if
                case (VN_SPD, VN_ZPD)
                case default
                  cycle
                end select
                !---------------------------
                ! mark rejected observations
                !---------------------------
                call decr_use (d% use, state = stat_new(hd% obstype), &
                                       check = CHK_BLACKLIST          )
                if (verbose > 1) then
                   call nextline
                   if (lpobs) then
                     write(oline(iol),'(a,a,f10.0)') &
                     '    rejected ', oq_name (varno_oq (obs(ib)% varno(i))), pobs
                   else
                     write(oline(iol),'(a,a,f10.0)') &
                     '    rejected ', oq_name (varno_oq (obs(ib)% varno(i))), zobs
                   end if
                end if
              end do
            end do
          endif
        end select
      end do
    end do
    !-------------
    ! write report
    !-------------
    call flush_buf
    call add_line_pio ('')
  end subroutine check_black
!------------------------------------------------------------------------------
  subroutine check_gross (obs, fg)
  type (t_obs)    ,intent(in)    :: obs(:) ! observation data
  type (t_vector) ,intent(inout) :: fg     ! first guess
  !----------------------------------------------------------------------
  ! Gross error check. To be called from MEC.
  ! These checks are otherwise done by check_fg in the assimilation cycle
  !----------------------------------------------------------------------

    target                    :: obs
    integer                   :: ib, is, i ! indices
    type(t_obs)      ,pointer :: ob        ! pointer to observation box
    type(t_spot)     ,pointer :: si        ! pointer to report meta data
    type(t_datum)    ,pointer :: bd        ! pointer to report body data
    type(t_set)      ,pointer :: set       ! use flags

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      ob   => obs(ib)
      do is = 1, ob% n_spot
        si  => ob% spot(is)
        do i = si%o% i + 1, si%o% i + si%o% n
          bd  => ob% body(i)
          set => bd% set

          !-------------------------------------------------
          ! 3) gross check not applicable (if fg is invalid)
          !-------------------------------------------------
          if (fg% s(ib)% x(i) == rvind) cycle

          !------------------------------------
          ! 4) bnd_fg: check bounds on fg value
          !------------------------------------
          select case (ob% varno(i))
          case default
            !--------------------------------
            ! DISMISS negative bending angles
            !--------------------------------
            if (VN_BENDANG == ob% varno(i) .and. fg% s(ib)% x(i) <= 0._wp) &
              fg% s(ib)% x(i) = rvind
            !---------------------------------------------
            ! Reject non-wind observation if out of bounds
            !---------------------------------------------
            if (fg% s(ib)% x(i) < set% bnd_fg(1)  .or. &
                fg% s(ib)% x(i) > set% bnd_fg(2))      &
              fg% s(ib)% x(i) = rvind
          case (VN_V, VN_U, VN_V10M, VN_U10M, VN_HLOS, VN_FF)
            !-------------------------------------------------------
            ! Reject low-wind observations for high-wind first-guess
            !-------------------------------------------------------
            if (abs (bd% o)           < set% bnd_fg(1)  .and. &
                abs (fg% s(ib)% x(i)) > set% bnd_fg(2))       &
              fg% s(ib)% x(i) = rvind
          end select

        end do
      end do
    end do

  end subroutine check_gross
!------------------------------------------------------------------------------
  subroutine estimate_fgchk (obs, o, e_fgchk)
  type (t_obs)    ,intent(in)    :: obs(:) ! observation meta data
  type (t_vector) ,intent(in)    :: o      ! observed values
  type (t_vector) ,intent(inout) :: e_fgchk! additional error for fg check
  !------------------------------------------------------------------------
  ! CS: estimate additional error for fg check;
  ! this error depends on the stability of the observed temperature profile
  ! (thus only available for multi-level temperature reports) and is
  ! intended to reflect a systematic error which is not accounted for
  ! in the usual estimation of the fg error e.g. by the ensemble spread
  !------------------------------------------------------------------------
    type(t_obs)      ,pointer :: ob        ! pointer to observation box
    type(t_spot)     ,pointer :: si        ! pointer to report meta data
    type(t_datum)    ,pointer :: bd        ! pointer to report body data
    target                    :: obs
    integer        ,parameter :: maxlev    = 1000
    real(wp)       ,parameter :: zlapsmean = -.0065_wp
    real(wp)       ,parameter :: dlplim    = 0.03_wp  ! thickness (as ln(p))
                                           ! of layer within which enhancement
                                           ! factor computation is iterated
                                           ! (about 30 hPa at 1000 hPa)
    real(wp) :: pf_t(maxlev), pf_rh(maxlev), pf_p(maxlev), pf_lp(maxlev)
    real(wp) :: zdt (maxlev), zdq (maxlev) ! additional T,RH error for fg check
    real(wp) :: cb_t(maxlev), cb_q(maxlev) ! inverse of (f.g.) check bounds
                                           !   (typically = 1/3.)
    integer  :: ix_t(maxlev), ix_q(maxlev) ! index used to assign the obs 'pf_t'
                                           !   to the vector segment 'x'
    real(wp) :: z1ddz,zf! weight factor for vertical interpol. to the obs. point
    real(wp) :: zlapse  ! (approx) lapse rate
    real(wp) :: zddt    ! difference betw. T(ilev+1) and T extrapolated from
                        !   T(ilev) using lapse rate limit 'zlapsmean'(-0.0065K/m)
                        !   this depends on vertical distance over which the
                        !   lapse rate is computed (Delta_beta T)
    real(wp) :: zinvrs  ! inversion in [K] between considered levels
    real(wp) :: zstabf  ! stability factor: = 0 if less stable than '-0.0065 K/m'
                        !                   = 1 if inversion ; linear in between
    real(wp) :: zfddt   ! stability factor enhancement depending on 'zddt'
    real(wp) :: zinvla  ! scaled positive lapse rate at inversions plus 1
    real(wp) :: zstabt, zstabq ! term for additional T,RH error due to stable
                               !   layer without inversion
    real(wp) :: zdtinv, zdqinv ! additional T,RH error for fg check
    real(wp) :: dlpmin  ! min. distance (ln(p)) betw. level in current iteration
    integer  :: nlev    ! number of vertical levels in report
    integer  :: ib, is, i, k, ilev, ilev2, idlev ! indices
  !-------------------------------------------------------

    if (fgchk_inv_t == 0._wp .and. &
        fgchk_inv_q == 0._wp       ) return

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      ob => obs(ib)
      ! loop over reports
      do is = 1, ob% n_spot
        si  => ob% spot(is)
!       print *,'### e_fgchk: ib,is,ot',ib, obs(ib)% pe, is, si% hd% obstype
        if (si% hd% obstype /= OT_TEMP) cycle
        if (si% col% nlev   <= 3      ) cycle
        pf_p (:) = 0._wp
        pf_t (:) = 0._wp
        pf_rh(:) = 0._wp
        cb_t (:) = 1._wp
        cb_q (:) = 1._wp
        ilev    = 0
        ! loop over observations of the report
        do k  = 1,  si%o% n
          i   = k + si%o% i
          bd  => ob% body(i)
!         set => bd% set
!         print *,'### e_fgchk: k,name,plev,olev,varno,obs,state', k         &
!                , si% statid, bd% plev, ob% olev(i), ob% varno(i)           &
!                , o% s(ib)% x(i), si% use% state, bd% use% state
!         print *,'### e_fgchk: k,t_index_l',k,si% statid,si% l,si% col% nlev

          ! fill usable T, RH obs in distict arrays 'pf_t', 'pf_rh'
          ! (array index is level index, where level pressure is in 'pf_p')
          if (      (     (ob% varno(i) == VN_T )                            &
                     .or. (ob% varno(i) == VN_RH))                           &
              .and. (     (bd% use% state >  STAT_REJECTED)                  &
                     .or. (bd% use% state == STAT_PASSIVE ))) then
            if (abs( bd% plev - pf_p(max(ilev,1)) ) > 0.1_wp) then
              ilev = ilev + 1
              pf_p (ilev) = bd% plev
              ix_t (ilev) = 0
              ix_q (ilev) = 0
            endif
            if (ob% varno(i) == VN_T ) then
              pf_t (ilev) = o% s(ib)% x(i)
              ix_t (ilev) = i
              cb_t (ilev) = 1._wp/ bd% set% sgm_fg(1)
            elseif (ob% varno(i) == VN_RH) then
              pf_rh(ilev) = o% s(ib)% x(i)
              ix_q (ilev) = i
              cb_q (ilev) = 1._wp/ bd% set% sgm_fg(1)
            endif
          endif
        enddo
        nlev = ilev
        ! discard levels without observed temperature, and compute LOG( p )
        idlev = 0
        do ilev = 1 , nlev
          if (pf_t(ilev) > 1._wp) then
            pf_p (ilev-idlev) = pf_p (ilev)
            pf_t (ilev-idlev) = pf_t (ilev)
            pf_rh(ilev-idlev) = pf_rh(ilev)
            ix_t (ilev-idlev) = ix_t (ilev)
            ix_q (ilev-idlev) = ix_q (ilev)
            cb_t (ilev-idlev) = cb_t (ilev)
            cb_q (ilev-idlev) = cb_q (ilev)
            pf_lp(ilev-idlev) = log( pf_p(ilev-idlev) )
            zdq  (ilev-idlev) = 0._wp
            zdt  (ilev-idlev) = 0._wp
          else
            idlev = idlev + 1
          endif
        enddo
        nlev = nlev - idlev
        iterate_gap: do idlev = 1 , 7
          dlpmin = 2._wp * dlplim
          do ilev = 1 , nlev-idlev
            ilev2 = ilev + idlev
            IF ((idlev == 1) .OR. (pf_lp(ilev)-pf_lp(ilev2) <= dlplim)) THEN
              z1ddz  = - gacc/rd *2._wp /(pf_t (ilev2) + pf_t (ilev))          &
                                        /(pf_lp(ilev2) - pf_lp(ilev))
              zlapse = z1ddz            *(pf_t (ilev2) - pf_t (ilev))
              zddt   = max( 0._wp ,   pf_t(ilev2)                              &
                                   - (pf_t(ilev ) + zlapsmean /z1ddz) )
              zinvrs = max( 0._wp, pf_t(ilev2) - pf_t(ilev) )
              zstabf = max( 0._wp, min( 1._wp , 1._wp - zlapse / zlapsmean ) )

              ! humidity: see COSMO Newsletter 6, slightly revised
              ! --------
              zfddt  = 1._wp + zddt/(1._wp + zddt)
              !   0 <= zstabq < 0.5  for thin  isothermal layer
              !        zstabq = 0.25 for thick isothermal layer
              zstabq = 0.25_wp * zstabf * zfddt

              !   isothermal = 1 <= zinvla <= 2 = sharp inversions
              zinvla = max( 0._wp, min( 1._wp , zlapse / 0.05_wp ) ) + 1._wp
              !     for fgchk_inv_q == 0.6:
              !   zdqinv = 15% RH for thin  isothermal layer
              !   zdqinv = 30% RH for thick isothermal layer
              !   zdqinv = 75% RH for inversion >= 5K or sharp inversion >= 2.5K
              zdqinv = min( 0.75_wp, fgchk_inv_q *(  zstabq                    &
                                                   + 0.2_wp* zinvrs* zinvla) )

              ! temperature
              ! -----------
              zfddt  = 1._wp + 2._wp* zddt/(5._wp + zddt)
              !   values for zero lapse rate (zstabf = 1):
              !     zstabt = 1.5K for Delta_beta T small,i.e. thin layer
              !     zstabt =  2 K for Delta_beta T = 1K, i.e. Delta_z = 150m
              !     zstabt =  3 K for Delta_beta T = 5K, i.e. Delta_z = 750m
              !     zstabt = 3.5K for Delta_beta T =10K, i.e. Delta_z =1500m
              zstabt = 1.5_wp * zstabf * zfddt
              !     for fgchk_inv_t == 0.8:
              !   zdtinv = 1.5K for thin  isothermal layer or inversion <~ 2 K
              !   zdtinv =  3 K for thick isothermal layer or inversion = 3.75 K
              !   zdtinv = 12 K for inversion >= 15 K
!             zdtinv = MIN( 12._wp , MAX( zstabt , 0.8_wp* zinvrs ) )
              zdtinv = min( 12._wp,  fgchk_inv_t * max( zstabt , zinvrs ) )

              !   (later (in subr. 'check_fg'),
              !    the threshold for the first guess check is set to
              !    (set% sgm_fg(1) = 1/cb_x = 3)*sigma (sigma: estimated error);
              !    therefore multply zdx by cb_x here to get meaningful values
              !    in cases of strong inversions)
              zdq (ilev ) = max( zdq(ilev ) , zdqinv * cb_q(ilev ) )
              zdq (ilev2) = max( zdq(ilev2) , zdqinv * cb_q(ilev2) )
              zdt (ilev ) = max( zdt(ilev ) , zdtinv * cb_t(ilev ) )
              zdt (ilev2) = max( zdt(ilev2) , zdtinv * cb_t(ilev2) )
!             if (zdt(ilev) > 1.0_wp)                                          &
!               print *,'  e_fgchk,ZZ4: k,p,T,dTinv, ...',k, si% statid, idlev &
!                      , pf_p(ilev), pf_t(ilev), zdt(ilev)                     &
!                      , pf_p(ilev2),pf_t(ilev2),zdt(ilev2),zdtinv,zstabt,zinvrs
            endif
            dlpmin = min( dlpmin , pf_lp(ilev) - pf_lp(ilev2) )
          enddo
          !   up to 7 iterations as long as minimum vertical distance < 25 hPa
          if (dlpmin > dlplim)   exit iterate_gap
        enddo iterate_gap
        do ilev = 1 , nlev
          if (pf_p(ilev) < 60000._wp) then
          !   above 600 hPa, reduce additional contributions by 50% for safety
          !     (and impose an upper limit of 1K instead of 4K for T-obs,
          !      because it's potentially dangerous to use T-obs themselves
          !      to enhance the threshold used for QC of T-obs:
          !      truly erroneous T-obs may introduce spurious inversions)
            zdq (ilev) =      zdq(ilev) * 0.5_wp
            zdt (ilev) = min( zdt(ilev) * 0.5_wp , cb_t(ilev) * 4._wp )
          elseif (pf_p(ilev) < 80000._wp) then
          !   within 600 - 800 hPa, taper towards full values
          !   below 800 hPa, use full additional contribution
            zf = (80000._wp - pf_p(ilev)) / 20000._wp
            zdq (ilev) =      zdq (ilev) * (zf *0.5_wp + (1._wp-zf))
            zdt (ilev) = min( zdt (ilev) * (zf *0.5_wp + (1._wp-zf))           &
                            , cb_t(ilev) * (zf *4.0_wp + (1._wp-zf)*12.0_wp ) )
          endif
        enddo

        ! set additional contribution to fg error
        do ilev = 1 , nlev
          ! e_fgchk and o are both type t_vector, therefore the same index 'i'
          ! is used for s(ib)% x(i)
          if (ix_t(ilev) > 0)  e_fgchk% s(ib)% x(ix_t(ilev)) = zdt(ilev)
          if (ix_q(ilev) > 0)  e_fgchk% s(ib)% x(ix_q(ilev)) = zdq(ilev)
          ! diagnostic print   (only large values, zdtinv>2.5K or zdqinv>0.45)
          if (      (ix_t(ilev) > 0) .and. (ix_q(ilev) > 0)                    &
              .and. (     (zdt(ilev) > 2.5_wp *cb_t(ilev))                     &
                     .or. (zdq(ilev) > 0.45_wp*cb_q(ilev)))) then
            print '("  e_fgchk: k,name,p,2x(varno,obs,e_fgchk)",i4,2x,a,f9.0   &
                  &,i4,f8.2,f6.2,i4,2f6.3)', k, si% statid                     &
                   , ob% olev(max(ix_t(ilev),ix_q(ilev)))                      &
                   , ob% varno(ix_t(ilev)), o% s(ib)% x(ix_t(ilev))            &
                                    , e_fgchk% s(ib)% x(ix_t(ilev))            &
                   , ob% varno(ix_q(ilev)), o% s(ib)% x(ix_q(ilev))            &
                                    , e_fgchk% s(ib)% x(ix_q(ilev))
!                  , ob% varno(ix_q(ilev)), o% s(ib)% x(ix_q(ilev)), zdq(ilev)
          elseif ((ix_t(ilev) > 0) .and. (zdt(ilev) > 2.5_wp*cb_t(ilev))) then
            print '("  e_fgchk: k,name,p,   varno,obs,e_fgchk ",i4,2x,a,f9.0   &
                  &,i4,f8.2,f6.2         )', k, si% statid                     &
                   , ob% olev(max(ix_t(ilev),ix_q(ilev)))                      &
                   , ob% varno(ix_t(ilev)), o% s(ib)% x(ix_t(ilev)), zdt(ilev)
          elseif ((ix_q(ilev) > 0) .and. (zdq(ilev) > 0.45_wp*cb_q(ilev))) then
            print '("  e_fgchk: k,name,p,   varno,obs,e_fgchk ",i4,2x,a,f9.0   &
                  &,18x         ,i4,2f6.3)', k, si% statid                     &
                   , ob% olev(max(ix_t(ilev),ix_q(ilev)))                      &
                   , ob% varno(ix_q(ilev)), o% s(ib)% x(ix_q(ilev)), zdq(ilev)
          endif
        end do
      end do
    end do

  end subroutine estimate_fgchk
!------------------------------------------------------------------------------
  subroutine check_fg (obs, o_fg, e_fg, e_o, fg, e_fgchk)
  type (t_obs)    ,intent(inout) :: obs(:)  ! observation data
  type (t_vector) ,intent(inout) :: o_fg    ! observation minus forecast
  type (t_vector) ,intent(in)    :: e_fg    ! first guess error
  type (t_vector) ,intent(in)    :: e_o     ! observational error
  type (t_vector) ,intent(in)    :: fg      ! first guess
  type (t_vector) ,intent(in)    :: e_fgchk ! additional term for fg check
  target                         :: obs
  optional                       :: e_fgchk
  !----------------------------------------------------------------------------
  ! Generic first guess check:
  ! Apply checks on first guess and observational values
  ! based on the settings in namelist /rules/.
  ! If a test fail, remaining checks are skipped.
  !
  !    /rules/ :
  ! 1) prc     : disturb observational data (artificial data)
  ! 2) bnd_obs : reject observations with too large observed values
  ! 3)         : first guess check not applicable (if fg is invalid)
  ! 4) bnd_fg  : check bounds on fg value
  ! 5)         : first guess check not applicable (fg-error <0 or obs-error<=0)
  ! 6) sigm_of : actual first guess check
  ! 7)         : GPS-RO specific: upper limit for e_fg
  ! 8) sgm_oc  : checks on observation error
  !----------------------------------------------------------------------------
    real(wp)                  :: sigm_of   ! |o-fg| / sqrt (e_o**2 + e_fg**2)
    real(wp)                  :: sigm_o    ! |o-fg| / e_o
    real(wp)                  :: abso      ! |o|
    real(wp)                  :: absfg     ! |fg|
    real(wp)                  :: abs_ofg   ! |o-fg|
    real(wp)                  :: eo        ! e_o        observational error
    real(wp)                  :: efg       ! e_fg       fg error or spread
    real(wp)                  :: eo_o      ! e_o / |o|
    real(wp)                  :: eo_efg    ! e_o / e_fg
    real(wp)                  :: eps       ! normal distributed random variable
    integer                   :: ib, is, i ! indices
    type(t_obs)      ,pointer :: ob        ! pointer to observation box
    type(t_spot)     ,pointer :: si        ! pointer to report meta data
    type(t_datum)    ,pointer :: bd        ! pointer to report body data
    type(random_state_t),pointer :: seed   ! pointer to random generator state
    integer                   :: status    ! active, passive, etc.
    integer                   :: stat_new  ! status for obs.error check
!   integer                   :: stat_old  ! old status
    type(t_set)      ,pointer :: set       ! use flags

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      ob   => obs(ib)
      seed => fg% info% b(ib)% seed
      do is = 1, ob% n_spot
        si  => ob% spot(is)
        do i = si%o% i + 1, si%o% i + si%o% n
          bd  => ob% body(i)
          set => bd% set
          !-----------------------------------------------
          ! check if any rule matched that sets the status
          !-----------------------------------------------
          if (set% use == empty_set% use) then
             call decr_use (bd% use, check=CHK_NOTUSED)
             cycle
          end if
          !-------------------------------------------
          ! for wind direction take obs-fg modular 360
          !-------------------------------------------
          select case (ob% varno(i))
          case (VN_DD)
            if (o_fg %s(ib)%x(i) < -180._wp)  o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) + 360._wp
            if (o_fg %s(ib)%x(i) >  180._wp)  o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) - 360._wp
          end select
          !-------------------------------------------------------
          ! 1) prc  : disturb observational data (artificial data)
          !           add first guess or observational error
          !-------------------------------------------------------
          if (set% prc /= iu1) then
            if (iand (set% prc, PRC_ADDFG) /= 0) then
              o_fg %s(ib)%x(i) = bd% o
              bd% o            = bd% o + fg %s(ib)%x(i)
            endif
            if (iand (set% prc, PRC_OBS_FG) /= 0) then
              o_fg %s(ib)%x(i) = 0._wp
              bd% o            = fg %s(ib)%x(i)
            endif
            if (iand (set% prc, PRC_ADD_1) /= 0) then
              o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) + 1._wp
              bd% o            = bd% o            + 1._sp
            endif
            if (iand (set% prc, PRC_ADD_R) /= 0) then
              o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) + e_o %s(ib)%x(i)
              bd% o            = bd% o            + e_o %s(ib)%x(i)
            endif
            if (iand (set% prc, PRC_SUB_R) /= 0) then
              o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) - e_o %s(ib)%x(i)
              bd% o            = bd% o            - e_o %s(ib)%x(i)
            endif
            if (iand (set% prc, PRC_ADDOBSER) /= 0) then
              call random_gauss (eps, seed)
              o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) + e_o %s(ib)%x(i) * eps
              bd% o            = bd% o            + e_o %s(ib)%x(i) * eps
            endif
            if (iand (set% prc, PRC_ADDOBSE2) /= 0) then
              call random_gauss (eps, seed) ;eps = eps / sqrt(2._wp)
              o_fg %s(ib)%x(i) = o_fg %s(ib)%x(i) + e_o %s(ib)%x(i) * eps
              bd% o            = bd% o            + e_o %s(ib)%x(i) * eps
            endif
          endif

          !----------------------------
          ! pre-compile some quantities
          !----------------------------
          eo         =         e_o % s(ib)% x(i)
          select case (ob% varno(i))
          case default
             abs_ofg = abs  (  o_fg% s(ib)% x(i))
             absfg   = abs  (    fg% s(ib)% x(i))
             abso    = abs  (  o_fg% s(ib)% x(i)   + fg% s(ib)% x(i)  )
             efg     =         e_fg% s(ib)% x(i)
          case (VN_V, VN_V10M) ! same as VN_U
          case (VN_U, VN_U10M)
             !----------------------------------------------
             ! Wind observation: apply checks to wind vector
             !----------------------------------------------
             abs_ofg = sqrt (  o_fg% s(ib)% x(i+1)**2 &
                            +  o_fg% s(ib)% x(i)  **2 )
             absfg   = sqrt (    fg% s(ib)% x(i+1)**2 &
                            +    fg% s(ib)% x(i)  **2 )
             abso    = sqrt ( (o_fg% s(ib)% x(i+1) + fg% s(ib)% x(i+1))**2 &
                            + (o_fg% s(ib)% x(i)   + fg% s(ib)% x(i)  )**2 )
             efg     = sqrt (  e_fg% s(ib)% x(i+1)**2 & ! for EnKF (efg=spread)
                            +  e_fg% s(ib)% x(i)  **2 ) / sqrt(2._wp)
          end select

          !----------------------------------------------------------------
          ! 2) bnd_obs : reject observations with too large observed values
          !              (does not apply to wind direction)
          !----------------------------------------------------------------
          if (abso > set% bnd_obs .and. ob% varno(i) /= VN_DD) then
             call decr_use (bd% use, check=CHK_FG)
             cycle
          end if

          !-------------------------------------------------------
          ! 3) first guess check not applicable (if fg is invalid)
          !-------------------------------------------------------
          if (fg% s(ib)% x(i) == rvind) then
            call decr_use (bd% use, STAT_PASSIVE, CHK_OPERATOR)
            cycle
          endif

          !------------------------------------
          ! 4) bnd_fg: check bounds on fg value
          !------------------------------------
          select case (ob% varno(i))
          case default
            !--------------------------------
            ! DISMISS negative bending angles
            !--------------------------------
            if (VN_BENDANG == ob% varno(i) .and. fg% s(ib)% x(i) <= 0._wp) then
              call decr_use (bd% use, STAT_DISMISS, CHK_FG)
            endif
            !---------------------------------------------
            ! Reject non-wind observation if out of bounds
            !---------------------------------------------
            if (fg% s(ib)% x(i) < set% bnd_fg(1)  .or. &
                fg% s(ib)% x(i) > set% bnd_fg(2)) then
              call decr_use (bd% use, check=CHK_FG)
              cycle
            endif
          case (VN_DD)
            !---------------------------------
            ! not applicable to wind direction
            !---------------------------------
          case (VN_V, VN_V10M, VN_U, VN_U10M, VN_HLOS, VN_FF)
            !-------------------------------------------------------
            ! Reject low-wind observations for high-wind first-guess
            !-------------------------------------------------------
            if (abso  < set% bnd_fg(1)  .and. &
                absfg > set% bnd_fg(2)) then
              call decr_use (bd% use, check=CHK_FG)
              cycle
            endif
          end select

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! for ADM experiments modify observation error in fg check
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! select case (ob% varno(i))
! case (VN_U, VN_V)
! if (si% hd% obstype /= OT_AIREP .or. ob% olev(i) > 30000._wp) eo = eo / sqrt(2._wp)
! end select

          !-------------------------------------
          ! 5) first guess check not applicable
          !    if fg-error < 0 or obs-error <= 0
          !-------------------------------------
          if (eo <= 0._wp) then
            call decr_use (bd% use, STAT_PASSIVE, CHK_FG)
            cycle
          endif
          if (efg <  0._wp) then
            call decr_use (bd% use, STAT_PASSIVE, CHK_FG)
            cycle
          endif

          !---------------------------------
          ! pre-compile some more quantities
          !---------------------------------
          if (present (e_fgchk)) then
            ! 'efg' estimated by ensemble spread accounts only for random errors;
            ! additional term 'e_fgchk' meant to account for systematic errors
            !   in T, RH in presence of (observed) inversions or stable layers
            sigm_of  = abs_ofg / sqrt (efg ** 2 + eo ** 2 +                   &
                                       e_fgchk% s(ib)% x(i) ** 2)
          else
            sigm_of  = abs_ofg / sqrt (efg ** 2 + eo ** 2)
          endif
          sigm_o   = abs_ofg / eo
          if (efg>0._wp) then
            eo_efg = eo / efg
          else                   ! first guess error may be zero for:
            eo_efg = huge(1._sp) ! wind at the poles, efg derived by randomisation
          endif                  ! flag adjoint=F in first guess scan
          eo_o     = abso / eo
          if (eo_o /= 0) then
            eo_o   = 1._wp / eo_o
          else
            eo_o   = sqrt(huge(1._sp))
          endif

          !-------------------------------------
          ! 6) sigm_of: actual first guess check
          !-------------------------------------
          status = STAT_ACTIVE_1
          stat_new = rept_use (si% hd% obstype)% use (CHK_FG)
          if (sigm_of >= set% sgm_fg(1)) status = min (status, stat_new)
!         if (ob% varno(i) /= VN_LWC) then
!           if (sigm_of >= set% sgm_fg(1)) status = min (status, stat_new)
!         elseif (ob% varno(i) == VN_LWC) then
!           if (abs_ofg >= 2.5*set% stdv_oi(si% stlsf)) status = min (status, stat_new)
!         endif
          if (sigm_o  >= set% sgm_fg(2)) status = min (status, STAT_ACTIVE)
          if (sigm_o  >= set% sgm_fg(3)) status = min (status, stat_new)
          if (status < STAT_ACTIVE_0I) then
            call decr_use (bd% use, status, CHK_FG)
          else
            call decr_use (bd% use, status, CHK_NONE)
          endif

          !-----------------------------------------
          ! 7) GPS-RO specific: upper limit for e_fg
          !-----------------------------------------
          if ((si% hd% obstype == OT_GPSRO) .and. (efg > 0.1_wp)) then
            call decr_use (bd% use, STAT_REJECTED, CHK_FG)
          endif

          !----------------------------------------
          ! 8) sgm_oc : checks on observation error
          !----------------------------------------
          if (eo     >= set% sgm_oc(1) .or. &
              eo_o   >= set% sgm_oc(2) .or. &
              eo_efg >= set% sgm_oc(3)) then
             status   = bd% use% state
             stat_new = rept_use (si% hd% obstype)% use (CHK_OBS_ERR)
             if (stat_new /= status) &
               call decr_use (bd% use, stat_new, CHK_OBS_ERR)
          end if

        end do
      end do
    end do
  end subroutine check_fg
!------------------------------------------------------------------------------
  subroutine check_eo (obs, e_fg, e_o)
  !---------------------------
  ! modify observational error
  !---------------------------
  type (t_obs_set) ,intent(inout) :: obs    ! observation data
  type (t_vector)  ,intent(in)    :: e_fg   ! first guess error
  type (t_vector)  ,intent(inout) :: e_o    ! observational error

    integer                   :: ib, is, i ! indices
    integer                   :: ii, j, jj ! indices
    type(t_obs)      ,pointer :: ob        ! pointer to observation box
    type(t_spot)     ,pointer :: si        ! pointer to report meta data
    type(t_datum)    ,pointer :: bd        ! pointer to report body data
    type(t_set)      ,pointer :: set       ! use flags
    real(wp)                  :: eo        ! e_o
    real(wp)                  :: efg       ! e_fg
    real(wp)                  :: eo_efg    ! e_o / e_fg
    real(wp)                  :: scale     ! scale factor for row/column of R
    integer                   :: i_, j_    ! loops over R-matrix (scaling R)
    integer                   :: i0, i1    ! required for scaling R
    integer                   :: ii_, jj_  ! required for scaling R

    do ib = 1, size(obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      ob => obs% o(ib)
      do is = 1, ob% n_spot
        si  => ob% spot(is)
        do i  = si%o% i + 1, si%o% i + si%o% n
          bd  => ob% body(i)
          set => bd% set
          ii  = i - si%o% i                          ! ++++++++ !
          do j=obs% R% b(ib,ib)% ia (i),      &
               obs% R% b(ib,ib)% ia (i + 1) - 1
             jj = obs% R% b(ib,ib)% ja (j) - si%o% i
             if (ii==jj) then
               eo     = sqrt (real (obs% R% b(ib,ib)% packed(j), wp))
!              eo     = e_o % s(ib)% x(i)
               select case (ob% varno(i))
               case (VN_U, VN_U10M)
                 efg = sqrt (  e_fg% s(ib)% x(i+1)**2 & ! for EnKF (efg=spread)
                            +  e_fg% s(ib)% x(i)  **2 ) / sqrt(2._wp)
               case (VN_V, VN_V10M) ! same as VN_U
               case default
                 efg  = e_fg% s(ib)% x(i)
               end select
               if (efg > 0._wp) then
                 eo_efg = eo / efg
                 if (eo_efg <  set% sgm_oc(4)) then
                   scale = (set% sgm_oc(4) * efg) / eo
                   eo = set% sgm_oc(4) * efg
                   if (eo >= set% sgm_oc(1)) &
                     call decr_use (bd% use, check = CHK_OBS_ERR)
                   e_o% s(ib)% x(i) = eo
                   obs% o(ib)% body(i)% eo     = eo
                   ! Scale corresponding row/column in R-matrix
                   ! (works for repr=csr/csc)
                   do i_  = si%o% i + 1, si%o% i + si%o% n
                     ii_  = i_ - si%o% i
                     i0 = obs% R% b(ib,ib)% ia (i_)
                     i1 = obs% R% b(ib,ib)% ia (i_ + 1) - 1
                     if (ii_ == ii) then
                       obs% R% b(ib,ib)% packed(i0:i1) = scale * obs% R% b(ib,ib)% packed(i0:i1)
                     end if
                     do j_ = i0, i1
                       jj_ = obs% R% b(ib,ib)% ja (j_) - si%o% i
                       if (jj_ == ii) then
                         obs% R% b(ib,ib)% packed(j_) = scale * obs% R% b(ib,ib)% packed(j_)
                       end if
                     end do
                   end do
                 endif !eo_efg < set%sgm_oc(4)
               endif ! efg > 0
             endif
          end do
        end do
      end do
    end do
  end subroutine check_eo
!------------------------------------------------------------------------------
  subroutine check_cons (obs)
  !------------------
  ! consistency check
  !------------------
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data

    target                    :: obs
    integer                   :: ib, is, i, l ! indices
    integer                   :: m(1)         ! result from minloc
    type(t_spot)     ,pointer :: si           ! pointer to report meta data
    type (t_obs)     ,pointer :: oi           ! pointer to observation box
    integer                   :: status       ! active, passive, etc.
    integer                   :: check        !
    integer                   :: i_u

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      oi => obs(ib)
      do is = 1, obs(ib)% n_spot
        si  => oi% spot(is)
        !-------------------------
        ! pairs of wind components
        !-------------------------
        i_u = 0
        do l = 2, si%o% n
          i = l + si%o% i
          if (((oi% varno(i-1) == VN_U    .and. oi% varno(i) == VN_V   )  .or.  &
               (oi% varno(i-1) == VN_U10M .and. oi% varno(i) == VN_V10M)) .and. &
              ( oi% olev (i-1) ==               oi% olev (i)            )       ) then
            if (oi% body (i-1)% use% state /=   oi% body (i)% use% state) then
               m      = minloc (oi% body(i-1:i)% use% state) -2
               status = oi% body(i+m(1))% use% state
               check  = oi% body(i+m(1))% use% check
               m(1)   = -1 - m(1)
               call decr_use (oi% body(i+m(1))% use, status, check)
            else
               status = oi% body(i)% use% state
               check  = oi% body(i)% use% check
            end if
            i_u = i-1
            cycle
          endif
          !--------------------------------
          ! applies also to adjacent FF, DD
          !--------------------------------
          if (i_u /= 0) then
            select case (oi% varno(i))
            case (VN_FF, VN_DD)
              if (oi% olev (i_u) == oi% olev (i)) then
                call decr_use (oi% body(i)% use, status, check)
              endif
            case default
              i_u = 0
            end select
          endif
        end do
      end do
    end do

  end subroutine check_cons
!------------------------------------------------------------------------------
  subroutine check_suff (obs)
  !------------------------------------------
  ! final check for sufficient data in report
  !------------------------------------------
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data

    integer               :: ib, is, i       ! indices
    type(t_spot) ,pointer :: si              ! pointer to report meta data
    integer               :: status, ostatus ! active, passive, etc.
    integer               :: check, chkb
    logical               :: passive         ! passive observation present

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      do is = 1, obs(ib)% n_spot
        si      => obs(ib)% spot(is)
        check   =  CHK_NO_OBS
        status  =  0
        passive = .false.
        do i = si%o%i+1, si%o%i+si%o%n
          ostatus = obs(ib)% body(i)% use% state
          passive = passive .or. ostatus == STAT_PASSIVE
          !----------------------------------
          ! derive highest observation status
          !----------------------------------
          if (status < ostatus) then
              status = ostatus
              chkb   = obs(ib)% body(i)% use% check
              if (chkb /= CHK_NO_OBS) check = chkb
          endif
        end do
        !---------------------------------------------------------------
        ! do not set to 'rejected' if 'passive' observations are present
        !---------------------------------------------------------------
        if (status == STAT_REJECTED .and. &
            passive .and. fix_no_obs      ) status = STAT_NOTACTIVE
        !------------------------------
        ! update report status variable
        !------------------------------
        call decr_rpt_use (si, CHK_NO_OBS, status, chk(check)% mnem)
      end do
    end do

  end subroutine check_suff
!------------------------------------------------------------------------------
  subroutine check_redo (obs)
  !---------------------------------------------------------------
  ! redo the checks on the report based on the bits set in 'flags'
  !---------------------------------------------------------------
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data

    integer               :: ib, is, i ! indices
    type(t_spot) ,pointer :: si        ! pointer to report meta data
!   integer               :: state     ! active, passive, etc.
    integer               :: check

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      do is = 1, obs(ib)% n_spot
        si => obs(ib)% spot(is)
        !-----------------------------------
        ! update observation status variable
        !-----------------------------------
        check   = si% use% check
!       state   = si% use% state

        do i = 1, n_chk
          if (i == CHK_NO_OBS) cycle
          if (btest (si% use% flags, i)) then
            call decr_rpt_use (si, check, comment='redo check')
          endif
        end do

!if (si% use% state /= state) then
!print *,dace% pe,'### check_redo: ',                    &
!stat_mnem (state),' -> ',stat_mnem (si% use% state),' , ',&
!chk(check)% mnem, ' -> ',       chk(si% use% check)% mnem
!endif

      end do
    end do

  end subroutine check_redo
!------------------------------------------------------------------------------
  subroutine check_obstype (obs, list)
  !-------------------------------------------
  ! check that obstype is in list given by MEC
  !-------------------------------------------
  type(t_obs)      ,intent(inout) :: obs(:) ! observation data
  character(len=*) ,intent(in)    :: list   ! list of obstypes
    integer               :: ib, is, io   ! indices
    type(t_spot) ,pointer :: si           ! pointer to report meta data

    if (list == '') return ! return for empty list

    do io = 1, n_ot
      if (index (list, trim(name_value (obstype, io))) <= 0) &
        rept_use(io)% use(:) = min (rept_use(io)% use(:), STAT_DISMISS)
    end do

    do ib = 1, size(obs)
      do is = 1, obs(ib)% n_spot
        si => obs(ib)% spot(is)
        !-----------------------------------
        ! update observation status variable
        !-----------------------------------
        if (index (list, trim(name_value (obstype, si% hd% obstype))) > 0) cycle
          call decr_rpt_use (si, CHK_OPERATOR, use= STAT_DISMISS, &
                                           comment= 'not in list' )
      end do
    end do

  end subroutine check_obstype
!------------------------------------------------------------------------------
  subroutine check_var  (obs)
  !------------------------------------------------------------------
  ! reject observations not handled by the variational scheme
  ! +++ preliminary code to handle observations provided by COMET +++
  !------------------------------------------------------------------
  type (t_obs) ,intent(inout) :: obs  ! observation data

    integer               :: is, i    ! indices
    type(t_spot) ,pointer :: si       ! pointer to report meta data

      do is = 1, obs% n_spot
        si => obs% spot(is)
        do i = si%o%i+1, si%o%i+si%o%n

          select case (obs% varno(i))
          !-------------------------
          ! ensure positive pressure
          !-------------------------
          case (VN_P, VN_PS)
            if (obs% body(i)% o <= 0._sp) &
              call decr_use (obs% body(i)% use, STAT_DISMISS, CHK_OPERATOR)
          !------------------------------------------------------
          ! replace pseudo relative humidity by relative humidity
          !------------------------------------------------------
          case (VN_PRH)
            call decr_use (obs% body(i)% use, STAT_OBS_ONLY, CHK_OPERATOR)
            obs% varno(i) = VN_RH
          !------------------------------
          ! remove vertical wind velocity
          !------------------------------
          case (VN_W)
            call decr_use (obs% body(i)% use, STAT_DISMISS, CHK_OPERATOR)
          case (VN_T, VN_T2M)
            select case (si% hd% obstype)
            !-----------------------------------------------
            ! remove temperature observations in AMV reports
            !-----------------------------------------------
            case (OT_SATOB)
!             call decr_use (obs% body(i)% use, STAT_DISMISS, CHK_OPERATOR)
              if (i /= si%o%i+3) call finish('check_var','SATOB T not at position 3')
              if (3 /= si%o%n  ) call finish('check_var','SATOB size /= 3')
              si%o%n = 2
              cycle
            !------------------------------------------------------
            ! set temperature and dewpoint passive in SYNOP reports
            ! temperature not implemented
            ! use RH instead of Td
            !------------------------------------------------------
            case (OT_SYNOP)
              if (si% hd% buf_type      /= 0 .or.    &! not SYNOP land
                  obs% body(i)% lev_typ /= VN_HOSAG) &! sensor height not set
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_OPERATOR)
              cycle
            case (OT_DRIBU)
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_OPERATOR)
              cycle
            end select
          case (VN_TD, VN_TD2M)
            select case (si% hd% obstype)
            case (OT_SYNOP, OT_DRIBU)
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_OPERATOR)
              cycle
            end select
          !-------------------------------------------------
          ! set FF, DD passive in SYNOP, BUOY, SCATT reports
          ! U, V should be used
          !-------------------------------------------------
          case (VN_FF, VN_DD)
            select case (si% hd% obstype)
            case (OT_SYNOP, OT_DRIBU)
              if (obs% body(i)% use% state == STAT_PASSIVE) cycle
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_OPERATOR)
              cycle
            case (OT_SCATT)
              !---------------------------------------------------------------
              ! SCATT: keep FF if there is only this observation in the report
              ! retain "check" if already passive
              !---------------------------------------------------------------
              if (obs% body(i)% use% state == STAT_PASSIVE) cycle
              if (obs% varno(i) == VN_DD) then
                call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_OPERATOR)
                cycle
              else if (si%o%n > 1) then
                call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_OPERATOR)
                cycle
              endif
            end select
          end select
        end do
      end do

  end subroutine check_var
!------------------------------------------------------------------------------
  subroutine check_obs (obs)
  !-----------------------
  ! check for valid report
  !-----------------------
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data

    integer               :: ib, is, i ! indices
    type(t_spot) ,pointer :: si        ! pointer to report meta data
    integer               :: status    ! active, passive, etc.
    integer               :: check

    do ib = 1, size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      do is = 1, obs(ib)% n_spot
        si => obs(ib)% spot(is)
        !-----------------------------------
        ! update observation status variable
        !-----------------------------------
        status = si% use% state
        check  = si% use% check
        do i = si%o%i+1, si%o%i+si%o%n
          call decr_use (obs(ib)% body(i)% use, status, check)
        end do
      end do
    end do

  end subroutine check_obs
!------------------------------------------------------------------------------
  subroutine check_accept (obs, w_qc)
  !-----------------------------------
  ! final check of acceptance of datum
  !-----------------------------------
  type (t_obs)    ,intent(inout) :: obs(:) ! observation data
  type (t_vector) ,intent(inout) :: w_qc   ! VQC weight

    integer               :: ib, is, i     ! indices
    type(t_spot) ,pointer :: si            ! pointer to report meta data
    integer               :: status        ! active, passive, etc.
    integer               :: check, chkb   ! check id
    integer               :: nact, nacc    ! count vqc-rejected, accepted
    character(len=32)     :: comment

    do ib = 1, size(obs)
      !------------------------------------------------
      ! set acceptance flag in dependence on VQC weight
      !------------------------------------------------
      if (obs(ib)% pe == dace% pe) then
        do i = 1, obs(ib)% n_obs
          if (obs(ib)% body(i)% use% state >= STAT_ACTIVE_0I .and. &
              w_qc% s(ib)% x(i) > 0.5_wp)  then
            obs(ib)% body(i)% use% state = STAT_ACCEPTED
            obs(ib)% body(i)% use% flags = &
              ibset (obs(ib)% body(i)% use% flags, CHK_FINAL)
          endif
        end do
      endif
      !------------------
      ! broadcast results
      !------------------
      if (obs(ib)% n_obs > 0) then
        call p_bcast (obs(ib)% body% use% flags, obs(ib)% pe)
        call p_bcast (obs(ib)% body% use% state, obs(ib)% pe)
      endif
      !--------------------------------------------------
      ! derive highest observation status for each report
      !--------------------------------------------------
      if (obs(ib)% pe == dace% pe) then
        do is = 1, obs(ib)% n_spot
          si => obs(ib)% spot(is)
          check  = CHK_NO_OBS
          status = 0
          nact   = 0; nacc = 0;
          do i = si%o%i+1, si%o%i+si%o%n
            if (status < obs(ib)% body(i)% use% state) then
                status = obs(ib)% body(i)% use% state
                chkb   = obs(ib)% body(i)% use% check
                if (chkb /= CHK_NO_OBS) check = chkb
            endif
            if (obs(ib)% body(i)% use% state >= STAT_ACTIVE_0I) nact=nact+1
            if (obs(ib)% body(i)% use% state >= STAT_ACCEPTED)  nacc=nacc+1
          end do
          !------------------------------
          ! update report status variable
          !------------------------------
!         call decr_rpt_use (si, check, status)
          call decr_rpt_use (si, CHK_NO_OBS, status, "chk="//chk(check)% mnem)
          write (comment, '(a,3i4)') 'pas. vqcrej. acc.:', &
            si%o%n-nact, nact-nacc, nacc
          call write_report (si, comment)
        end do
      endif
      !------------------------------
      ! write messages to report file
      !------------------------------
      call write_pending
    end do
  end subroutine check_accept
!==============================================================================
  subroutine check_nofgscan (obs)
  !---------------------------------------------------------------
  ! This routine is called in case the first guess scan is skipped
  !---------------------------------------------------------------
  type (t_obs)    ,intent(inout) :: obs ! observation data

    call decr_use (obs% body(:obs%n_obs)% use, STAT_ACTIVE, CHK_NONE)

  end subroutine check_nofgscan
!==============================================================================
  subroutine check_operator (obs, pass, n_invalid, H_atm, n_inval)
  type(t_obs_set) ,intent(inout) :: obs       ! observation data
  integer         ,intent(in)    :: pass      ! 1:fg 2:analysis 3:verification
  integer         ,intent(out)   :: n_invalid ! # observations to be removed
  type(t_vector)  ,intent(inout) :: H_atm     ! H appl. to atm
  integer         ,intent(out)   :: n_inval(:)! # invalid observations / type
  optional :: n_invalid, H_atm, n_inval
  !------------------------------------------------------------------------
  ! Check the applicability of observation operators.
  ! If not applicable the flag 'op_na' has been set to a nonzero value.
  ! In this case the following action is taken:
  !
  ! pass == 0 : do nothing
  !
  ! pass == 1 : First Guess Check
  !             Decrease the status flag of the observation.
  !
  ! pass == 2 : Analysis (to be checked within each outer PSAS iteration)
  !             Decrease the status flag of the observation.
  !             Count number 'n_invalid' of observations removed from
  !             the set of active (assimilated) variables.
  !             Unset mask for observations not used any more.
  !
  ! pass == 3 : Verification
  !             Set model equivalent 'H_atm' to NetCDF missing value
  !
  ! pass == 4 : COSMO LETKF
  !             If (H_atm == rvind) decrease the observation status flag.
  !------------------------------------------------------------------------

    integer               :: ib, is, i, n, o ! indices
    type(t_spot) ,pointer :: si              ! pointer to report meta data
    type(t_datum),pointer :: d(:)            ! pointer to body   meta data
    integer               :: status          ! active, passive, etc.
    integer               :: inv_op(6)       ! observation operator validity check
    integer               :: n_inv_op(n_ot)  ! invalid obs.operator applications
    integer               :: ot              ! temporary

    select case (pass)
    case (0)
       return
    case (2)
      if (.not. present (n_invalid))                                      &
        call finish ("check_operator","pass=2 requires argument n_invalid")
      n_invalid = 0
      n_inv_op  = 0
    case (3)
      if (.not. present (H_atm))                                      &
        call finish ("check_operator","pass=3 requires argument H_atm")
    case (4)
     if (.not. present (H_atm))                                      &
        call finish ("check_operator","pass=4 requires argument H_atm")
      do ib = 1, size(obs% o)
        if (obs% o (ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          si       => obs% o(ib)% spot(is)
          o        =  si% o% i
          n        =  si% o% n
          d        => obs% o(ib)% body (o+1 : o+n)
          do i = 1, n
            if (H_atm% s(ib)% x(o+i) == rvind) &
              call decr_use (d(i)% use, check=CHK_OPERATOR)
          end do
        end do
      end do
    end select

    obs% vc% active = 1._wp
    do ib = 1, size(obs% o)
      if (obs% o (ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        si       => obs% o(ib)% spot(is)
        ot       = si% hd% obstype
        inv_op   = rept_use (ot)% inv_op
        if (inv_op(1) < 0) cycle
        o        =  si% o% i
        n        =  si% o% n
        d        => obs% o(ib)% body (o+1 : o+n)
        select case (pass)
        case (1,2)
          status   = rept_use (ot)% use (CHK_OPERATOR)
          do i = 1, n
            if (    status <  d(i)% use% state .and. &
                any(inv_op == d(i)% op_na)           ) then
              if (pass == 2) n_inv_op(ot) = n_inv_op(ot) + 1
              if (pass == 2) write (0,*)                      &
                   "check_operator: statid,obs_id,box,i,op_na= ", &
                   si% statid, si% hd% id, ib, i, d(i)% op_na
              call decr_use (d(i)% use, check=CHK_OPERATOR)
            endif
            if (pass == 2 .and. d(i)% use% state < STAT_ACTIVE_0I) &
              obs% vc% active% s(ib)% x(o+i) = 0._wp
          end do
        case (3)
          do i = 1, n
            if (any(inv_op == d(i)% op_na)) then
              H_atm% s(ib)% x(o+i) = rvind
            endif
          end do
        end select
      end do
    end do
    if (pass == 2) then
       n_inv_op  = p_sum (n_inv_op)
       n_invalid = sum (n_inv_op)
       if (present (n_inval)) then
          n_inval = n_inv_op
       end if
    end if

  end subroutine check_operator
!==============================================================================
end module  mo_fg_checks
