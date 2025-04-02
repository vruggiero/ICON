!
!+ radiances observation operator specific routines
!
MODULE mo_tovs
!
! Description:
!   ATOVS (RAD) observation operator specific routines.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  bug fix in write_rttov_prof: write channel number only if flag 8 (H) is set
! V1_4         2009/03/26 Andreas Rhodin
!  fixes for fillvalues in 1dvar-files, for less than 5 AMSU-B channels present
! V1_7         2009/08/24 Andreas Rhodin
!  Changes for GME_L60
! V1_8         2009/12/09 Detlef Pingel
!  interface to RTTOV9
! V1_9         2010/04/20 Andreas Rhodin
!  changes for online bias-correction
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Detlef Pingel
!  changes for IASI, HIRS, SATPP, biascorrection
! V1_15        2011/12/06 Andreas Rhodin
!  option to interpolate stratospheric humidity in tv,q instead of tv,rh
! V1_19        2012-04-16 Andreas Rhodin
!  use flags flg_cld/flg_cld_... for steering of AMSU-A cloud flag
!  dismiss satellite zenith angles > 85 degree
!  bugfix for selection of minimal set of channels (H.Anlauf)
!  remove support for RTTOV-6                      (H.Anlauf)
!  reading coefficients on 1PE                     (H.Anlauf)
!  changes for RTTOV10                             (A.Messer)
! V1_20        2012-06-18 Harald Anlauf
!  Framework for arbitrary number of input feedback files
! V1_22        2013-02-13 Robin Faulwetter
!  work-around for discontinuities in RTTOV10,
!  fix bug in adjoint t2m,q2m interpolation for RTTOV7/10,
!  disable calculation of generalized inverse in cofRTOVP.         (A.Rhodin)
!  adjustments for ICON,
!  FTRACE-instrumentation of slow, non-parallelized part.          (H.Anlauf)
!  Merged some changes from Alexandre Lanciani into rtifc,
!  Adaptation to 50/51 levels for RTTOV,
!  Fix level count in for parallelization, (rttov_mult_prof=T)
!  Implemented loading of new style SATPP files, (Superobing)
!  Added ATMS instrument and NPP satellite to tables.              (A.Messer)
!  Fixed a bug caused by nonzero values in the H-operator,
!  Write biascor files also for pre-assimilation,
!  Improved selection of FOVs to be bias corrected, monitored,
!  implementation of vectorized K-mode                         (R.Faulwetter)
! V1_23        2013-03-26 Robin Faulwetter
!  Implemented processing of CrIS data
!  remove unused components p_surf,u_10m,v_10m (1dvar background)
!  bound constrained optimisation (variable transform, cloud top and fraction)
! V1_26        2013/06/27 Robin Faulwetter
!  Generalized the 'chan'-entry in TOVS_OBS_CHAN namelists.
!  Introduced a check on the influence of the surface onto radiances.
!  Introduced USE_MWSURF bit. Corrected the usage of other USE_* bits.
! V1_27        2013-11-08 Robin Faulwetter
!  Implemented FOV-dependent obserrors; Rewrite of radiance flagging
! V1_28        2014/02/26 Andreas Rhodin
!  changes for control variables: cloud-top-height/fraction, IR emmisivity PCs
! V1_29        2014/04/02 Harald Anlauf
!  process_tovs_mult: array bound violation bugfix
!  Hack to allow complete disabling of RTTOV
!  use PC fg coefficients for IR emissivity calculations
! V1_31        2014-08-21 Andreas Rhodin
!  /tovs_obs/ set bound_constr default to 1 (quadratic transform).
!  changes for IR emissivity PC approach.
!  Unify mo_rad with COSMO. Improve mo_rtifc. New write_rttov_prof routine.
! V1_35        2014-11-07 Harald Anlauf
!  fill_rad: bugfix for no radiances (detected by Intel v15)
! V1_37        2014-12-23 Andreas Rhodin
!  changes for Variational Ensemble Kalman Filter (VarEnKF)
! V1_42        2015-06-08 Andreas Rhodin
!  temporal interpolation for MEC; pass op_na flag from RTTOV to 3dvar
! V1_43        2015-08-19 Robin Faulwetter
!  Added features required for high peaking radiances over land/clouds
! V1_44        2015-09-30 Andreas Rhodin
!  load_tovs: abort if data is not present (with correct size)
! V1_45        2015-12-15 Andreas Rhodin
!  revised determination of localisation pressure level
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_47        2016-06-06 Robin Faulwetter
!  Many improvements for radiances.
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented RTTOV12
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured cloud detection for radiances. Added features to rttov12.
! V1_51        2017-02-24 Andreas Rhodin
!  change meaning of sink variable mode
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin    DWD  2001-2008  original source
! Detlef Pingel     DWD  2004       update to RTTOV7, consistent forward,Jacobi
!                        2004       write artificial data, Ts dummy variable
!                        2004       dummy sink variables for upper levels
! Oliver Schmid     DWD  2005       new obs data type
! Andreas Rhodin    DWD  2005       Write profiles for monitoring
! Christina Koepken DWD  2005       corrections
! Christina Koepken DWD  2007       introduce NOAA18,METOP
! Christina Koepken DWD  2008       introduce AMSU-B,MHS
!==============================================================================
!
! TODO
! Future support of individual RTTOV options per instrument:
! The array iopts will provide the instrument settings but currently just one
! instrument is supported and the settings are given in iopts(1). This needs to
! be replaced/checked as soon as options for several instruments are supported.
! => check iopts, iopts(1), ...
!==============================================================================
#if defined (_FTRACE)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use iso_fortran_env,  only: stderr => error_unit               ! Standard error
  use mo_cpu_time,      only: stop_time                          ! determine cpu and wall time
  use mo_exception,     only: finish,                           &! abort routine
                              message                            ! write message to stderr
  use mo_kind,          only: wp, sp, i8                         ! kind parameters
  use mo_run_params,    only: nex,                              &! experiment number
                              model,                            &! forecast model
                              ana_time,                         &! analysis time
                              run_time,                         &! time of run of this program
                              obsinput,                         &! obs. input path
                              input,                            &! input path
                              data,                             &! constant input data directory
                              aux,                              &! auxiliary output path
                              path_file,                        &! concatenate path/filename
                              interp_strato,                    &! flag:interpolate IFS stratosphere
                              method
  use mo_mpi_dace,      only: dace,                             &! MPI group info
                              p_bcast,                          &! broadcast routine
                              p_max,                            &! MPI generic max routine
                              p_min,                            &! MPI generic min routine
                              p_ior,                            &! MPI generic ior routine
                              p_or,                             &! MPI generic or  routine
                              p_sum,                            &! MPI generic sum routine
                              p_gatherv,                        &! MPI gather routine
                              p_alltoall,                       &! alltoall communcation
                              p_scatterv
  use mo_namelist,      only: position_nml,                     &! position namelist
                              nnml,                             &! namelist Fortran unit number
                              POSITIONED                         ! ok    code from position_nml
  use mo_obs_set,       only: t_obs_set,                        &!
                              t_obs_block,                      &!
                              obs_block                          !
  use mo_time,          only: t_time,                           &! Data type to hold date and time
                              cyyyymmddhh,                      &! get time as character string
                              cyyyymmddhhmm,                    &! get time as character string
                              init_time,                        &!
                              imm                                ! get month from t-time
  use mo_grid_intpol,   only: idx_init                           ! determine model grid indices
  use mo_obs_rules,     only: get_rule,                         &! routine to get a rule
                              t_set,                            &! result data type
                              iud                                ! undefined integer value
  use mo_vqc,           only: svqc,                             &! default var.qual.cntrl stdev.
                              vqc_form                           ! formulation of obs-costfunction
  use mo_instrid,       only: rttov_instr,                      &! RTTOV number from WMO instrument id
                              instr_rttov,                      &! WMO instrument id from RTTOV number
                              hss_instr,                        &! Hyperspectral sounder?
                              mw_instr                           ! MW instrument or not
  use mo_dace_string,   only: char3,                            &! transform integer to string
                              tolower      
  use mo_p_output,      only: add_line,                         &! add line to output buffer
                              add_line_pio,                     &! add line from IO processor
                              flush_buf                          ! write output buffer
#if (_RTTOV_VERSION >= 12)
  use mo_fortran_units, only: get_unit_number,                  &! reserve a unit number
                              return_unit_number                 ! release the unit number
#endif
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_grid,      only: VCT_P_ISO,                        &! isobaric pressure coordinate
                              VCT_P_HYB,                        &! hybrid   pressure coordinate
                              VCT_Z_HYB,                        &! hybrid      z     coordinate
                              VCT_Z_GEN,                        &! generalized z     coordinate
                              MO_ICON,                          &! ICON  model
                              MO_COSMO                           ! COSMO model
  use mo_atm_state,     only: t_atm                              ! atm. state data type
  use mo_t_col,         only: t_cols,                           &! model columns data type
                              t_col,                            &! component of t_cols
                              COL_TV,                           &! internal Temperature code
                              COL_T,                            &!
                              COL_RH,                           &!          rel.hum.    code
                              COL_Q,                            &!
                              COL_P,                            &!
                              COL_OZONE,                        &!
                              COL_CO2,                          &!
                              COL_CLC,                          &! cloud cover
                              COL_QCDIA,                        &! spec. cloud water, diag.
                              COL_QIDIA,                        &! spec. cloud ice, diag.
                              COL_REFF_QI,                      &! eff. radius cloud water
                              COL_REFF_QC                        ! eff. radius cloud ice
  use mo_physics,       only: rh_q,                             &! relative from specific humidity
                              gacc,                             &! g
                              rearth,                           &! radius of the earth
                              RDRD,                             &! gas const. (dry air)/gas const. (water)
#if (_RTTOV_VERSION >= 12)
                              ppmv_dry_q,                       &!         ppmv (dry air)    <- specific humidity
#endif
                              dppmv_dry_dq,                     &!deriv.:  ppmv (dry air)    <- specific humidity
                              q_ppmv_dry,                       &!         ppmv (dry air)    -> specific humidity
                              dppmv_moist_dq,                   &!deriv.:  ppmv (moist air)  <- specific humidity
                              q_ppmv_moist,                     &!         ppmv (moist air)  -> specific humidity
                              ngases,                           &!
                              d2r                                ! factor degree -> radians
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_t_obs,         only: t_obs,                            &!
                              t_spot,                           &!
                              t_head,                           &!
                              new_spot,                         &! reserve memory
                              new_obs,                          &! reserve memory
                              new_int,                          &! reserve memory
                              new_sink,                         &! reserve memory
                              set_xuv,                          &! set unit vector,solar zenith angle
                              set_sink,                         &! set sink variable parameters
                              destruct,                         &! destruct data type
                              TOVS,                             &! tovs data type id
                              ITY_ICOL,                         &! interpolation type: column
                              TSK_INIT,                         &! FLAGS: initialize module
                              TSK_READ,                         &!  read observation data
                              TSK_SETUP_COLS,                   &!  setup model columns
                              TSK_SET_CHR,                      &!  set observation characteristics
                              TSK_SETUP_FULL,                   &!  setup description of PSAS-space
                              TSK_SETUP_FUL0,                   &!  setup description of PSAS-space
                              TSK_SETUP_OP,                     &!  setup RTTOV
                              TSK_SHRINK,                       &!  release unused obs. in report
                              TSK_K,                            &!  setup H (Jacobi matrix)
                              TSK_R,                            &!   setup observational errors
                              TSK_Y,                            &!  run forward model
                              CHR_NONL,                         &! characteristics: nonlinear
                              CHR_EXP,                          &!                  expensive
                              OBS_TV,                           &! internal Temperature code
                              OBS_RH,                           &!          rel.hum.    code
                              OBS_Q,                            &!
                              OBS_HS,                           &!          geopotential height
                              OBS_DUM,                          &!          dummy variable
                              shrink_report,                    &! release unused channels
                              source,                           &! list   of source files
                              n_source,                         &! number of source files
                              m_source,                         &! max. number of source files
                              add_source,                       &! add source file to list
                              bufr_inv,                         &! observation file inventory
                              p_bcast,                          &!
                              FT_SATPP,                         &! SATPP file type
                              n_filetype,                       &! # of source files of specific type
                              t_ilev,                           &! level-dep. thinning information
                              match_ilev,                       &! check, whether data matches with t_thlev info
                              po_context,                       &! process_obs context
                              po_ilns,                          &! process_obs context: number of call within linesearch
                              po_lmon,                          &! process_obs context: monitoring call?
                              POC_FG,                           &! process_obs context: first guess
                              POC_ANA,                          &! process_obs context: analysis
                              POC_TOP,                          &! process_obs context: test obs.operator
                              POC_NN,                           &! process_obs context: not known
                              nwv_rad,                          &! type of vertical interpolation for radiances
                              vint_rttov,                       &! RTTOV interp. mode
                              int_nn,                           &! nearst neighbor horiz. interp.
                              npmax,                            &!
                              spt_hd_debug,                     &! debug selected spot(s)
                              spt_debug,                        &! ...
                              ldeb,usd,dpref,                   &! ...
                              ldeb_spot,                        &! ...
                              debug_spot                         ! ...
  use mo_fdbk_tables,   only: VN_RAWBT,                         &! brightness temperature
                              VN_REFL,                          &! reflectance
                              OT_RAD,                           &! Radiances report type ID
                              OF_BT_CLEAR_SKY,                  &! operator flag for clear-sky Tbs
                              SUR_SEA,                          &! surftype: sea
                              SUR_LAND,                         &! surftype: land
                              SUR_HIGHLAND,                     &! surftype: highland
                              SUR_MISMATCH,                     &! surftype: mismatch
                              SUR_MISSING,                      &! surftype: classification missing
                              SUR_MWSURF,                       &! surftype: problem detected by MW
                              SUR_BLK_PP,                       &! surftype: blacklisted by sat_pp
                              TF_EMIS_FASTEM,                   &! FASTEM identifier
                              TF_SURF_INFL,                     &! tovs_flag: surface influence
                              TF_DESERT_DUST,                   &! tovs_flag: desert dust
                              TF_VOLCANIC_ASH,                  &! tovs_flag: volcanic ash
                              TF_BC_FAILED,                     &! tovs_flag: BC failed
                              n_optv
  use mo_obs_sndrcv,    only: p_send,                           &! mpi-send variable(s) of type obs
                              p_recv                             ! mpi-recv variable(s) of type obs
  use mo_obs_tables,    only: check_report_0,                   &! simple checks on report validity
                              check_report_1,                   &!
                              decr_rpt_use,                     &! dismiss report
                              rept_use,                         &! status flag table
                              write_pending                      ! write pending dismissed reports
  use mo_t_use,         only: t_use,                            &! status flag data type definition
                              decr_use,                         &! dismiss report or datum
                              incr_use,                         &! increase state of a datum
                              STAT_DISMISS,                     &! status flag: dismiss report
                              STAT_PASSIVE,                     &!              passive report
                              STAT_REJECTED,                    &!              rejected report
                              STAT_PAS_REJ,                     &!              rejected passive report
                              STAT_ACTIVE,                      &!              active report
                              CHK_SURF,                         &! ! correct surface ?
                              CHK_DATASET,                      &! dataset flags (e.g.1D-Var flag)
                              CHK_CLOUD,                        &! rejected due to cloud flag
                              CHK_DOMAIN,                       &! out of model domain
                              CHK_CONSIST,                      &! inconsistent data flag
                              CHK_NOTUSED,                      &! data not used
                              CHK_BIASCOR,                      &! no bias correction possible
                              CHK_OBS_ERR,                      &! obs/fg error check
                              CHK_OPERATOR,                     &! Operator not applicable
                              CHK_MERGE,                        &! Merged reports
                              n_stat
  use mo_hum_ana,       only: hum_ana_top
  !--------------------------------------------------------
  ! access radiance observation tata types, RTTOV interface
  !--------------------------------------------------------
  use mo_rttov,         only: set_rttov_vers,                   &! set RTTOV version to be used
                              destruct,                         &! destruct radiances data type
                              call_rttvi,                       &! call RTTOV initialisation routine
                              t_rttov_prof,                     &! rttov variable input data type
                              construct,                        &! construct t_rttov_prof
                              jplev,                            &! no. of pressure levels
                              preslev,                          &! reference pressure levels (Pa)
                              preshpa,                          &! reference pressure levels (hPa)
                              lnp,                              &! log (pressure levels)
                              call_rttov,                       &! interface to RTTOV routine
                              nsv,                              &! # of single level input variables
                              nbs,                              &! # of bound constrained sink var.
                              rt_humi,                          &! Switch for humidity input to RTTOV
                              rt_levs_exact,                    &! Use exact rttov levels from RTTOV directly
                              rt_hard_limits                     ! Check RTTOV hard limits
  use mo_rtifc,         only: rtifc_set_opts,                   &! Set RTTOV options
                              rtifc_get_opts,                   &! Get RTTOV options
                              rtifc_fill_input,                 &! Supply profiles to RTTOV
                              rtifc_direct,                     &! RTTOV forward operator
                              rtifc_k,                          &! RTTOV jacobian+forward operator
                              rtifc_errmsg,                     &! rttov interface version
                              rtifc_cleanup,                    &! rttov deallocation routine
#if (_RTTOV_VERSION >= 12)
                              rtifc_init_atlas,                 &! initialize atlases
                              rtifc_init_brdf_atlas,            &!
                              rtifc_emis_atlas,                 &! get atlas emissivity
                              rtifc_tskin_retrieve,             &! retrieve tskin
#endif
                              rtifc_l2c_god,                    &! god-corrected l2c
                              rtifc_coef_prop,                  &! coef properties
                              rtifc_check_nlevs,                &!
                              rtifc_alloc_mode,                 &!
                              rtifc_vers,                       &!
                              read1pe,                          &! Read/distrib. coeffs
                              nlevs_top,                        &! Additional RTTOV levels above user levels
                              OUT_ASB,                          &! write all-sky brighntess temp. in fdbk files
                              OUT_CSB,                          &! clear-sky
                              OUT_CSR,                          &! write clear-sky radiances in fdbk files
                              OUT_VIS,                          &! output of reflectances
                              default_salinity,                 &! default value for salinity in RTTOV
                              default_gas_units,                &! default value for gas unit used in RTTOV
                              god_par_file,                     &! generalized optical depth (god) transformation parameters
                              wr_god,                           &! write plot files for god transformation
                              out_path,                         &! path to write files from mo_rtifc
                              chk_god_ifc => chk_god,           &! flag to control god check
                              god_thresh_ifc => god_thresh,     &! Threshold For god check
                              chk_reg_lims_ifc => chk_reg_lims, &! Apply (RTTOV) regularization limits
                              chk_plim_t,                       &! p-limit to avoid application of app_reg_lims on T
                              chk_plim_q,                       &! p-limit to avoid application of app_reg_lims on q
                              min_od,                           &! minimum optical depth in RTTOV
                              atlas_single_inst,                &! load atlas for single instr. (IR)
                              rts_land,rts_sea,rts_ice           ! surface types (RTTOV definition)
  use mo_cloud_params,  only: param_amsub                        ! AMSUB cloud detection parameters
#if (_RTTOV_VERSION >= 12)
  use rttov_const,      only: gas_unit_specconc,                &! rttov gas unit: specific concentration
                              gas_unit_ppmvdry,                 &! rttov gas unit: ppmv over dry air
                              gas_unit_ppmv                      ! rttov gas unit: ppmv over moist air
#endif
  use mo_rad,           only: t_rad_set,                        &! specific. of instr.
                              t_radv,                           &! Complete satellite dataset
                              t_rad_iopt,                       &!
                              t_emis_opt,                       &! Emissivity options
                              read_satpp_feedbk,                &! read sat_pp files
                              set_indx,                         &!
                              chan_indx,                        &!
                              rinvalid,                         &!
                              m_instr,                          &! maximum number of instruments in sat dataset
                              m_chan,                           &! maximum number of channels
                              m_bd,                             &! maximum number of bands
                              destruct,                         &! destruct t_radv and t_rad_set
                              assignment(=),                    &! assignment of t_rad_set
                              print_rad_set,                    &! debug
                              read_tovs_obs_chan_nml,           &
                              link_rad,                         &! link radiance datasets (t_radv) to rad_set
                              reduced_rad_set,                  &! routine to eliminate channels from t_rad_set
                              rad_set,                          &! Radiance dataset descriptions
                              n_set,                            &! Number of radiance datasets
                              m_rad_set,                        &! Maximum number of radiance datasets
                              USE_PASSIVE,                      &! Use only passively
                              USE_BCOR,                         &! Keep channel always, since it is required for bias correction
                              USE_SURFINFL,                     &! Apply surface influence check
                              USE_L2C_CLD,                      &! Use only if l2c above tropopause
                              USE_L2C_SURF,                     &! Use onyl if l2c above surface
                              USE_MINTSURF,                     &! Use also of min(t_surf) criterion is not fulfilled
                              default_flags,                    &! default flags (used if NO tovs_obs_chan nml is available)
                              sat_sun_azimuth_angle,            &! calc. the azimuth angle between sun and satellite
                              ! defaults for t_rad_opts to be modified in /tovs_obs/:
                              n_max_calc_k_def,                 &! max.# of forward  calculations
                              n_max_prof_k_def,                 &! max.# of K-matrix calculations
                              n_max_calc_y_def,                 &! max.# of profiles
                              n_max_prof_y_def,                 &! max.# of profiles
                              quality_mode_def,                 &! default for quality mode
                              l2c_type_def,                     &! default for l2c_type
                              l2c_rel_lim_def,                  &! default for l2c_rel_lim
                              l2c_abs_lim_def,                  &! default for l2c_abs_lim
                              l2c_use_rad_def,                  &! default for l2c_use_rad
                              l2c_max_def,                      &! default for l2c_max
                              l2c_max_trop_def,                 &! default for l2c_max_trop
                              l2c_max_midlat_def,               &! default for l2c_max_midlat
                              l2c_max_polar_def,                &! default for l2c_max_polar
                              surf_infl_mode_def,               &! default for surf_infl_mode
                              max_surf_infl_def,                &! default for max_surf_infl
                              d_stemp_def,                      &! default for d_stemp
                              d_emiss_def,                      &! default for d_emiss
                              l_max_surf_hgt_def,               &! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
                              write_rtovp,                      &! routine for writing of *RTOVP.nc files
                              WRM_OPEN,                         &! flag for write_rtovp
                              WRM_WRITE,                        &! flag for write_rtovp
                              WRM_CLOSE,                        &! flag for write_rtovp
                              t_nc_attr,                        &! attribute for write_rtovp
                              err_msg,                          &! Error message
                              p_bcast,                          &
                              lev2chan,                         &! level to channel assignment
                              thin_superob_radv,                &!
                              instr_type,                       &!
                              ITYP_VIS, ITYP_IR,                &!
                              usd_rad => usd,                   &!
                              ! ltest,                          &!
                              OPTV_L2C,                         &!
                              OPTV_NWC_FLG,                     &!
                              OPTV_ORB_PH ,                     &!
                              OPTV_INS_TMP,                     &!
                              OPTV_CLD_FRC,                     &!
                              OPTV_CLD_FLG,                     &!
                              t_tskin_opt,                      &!
                              n_styp
  use mo_t_tovs,        only: t_tovs,                           &! observation operator specific type
                              t_tovs_instr,                     &! information on instruments in t_tovs
                              tpp,                              &!
                              set_size,                         &! set '..._size' module variables
                              construct,                        &! t_tovs constructor routine
                              destruct,                         &! t_tovs  destructor routine
                              load,                             &! load  t_tovs
                              store,                            &! store t_tovs
                              get_tovs_rs,                      &! get rad_set and instrument info for t_tovs
                              link_tovs_rs,                     &! link t_tovs to the corresponding rad_set
                              add_av_cont,                      &! add entries to t_tovs%av_cont
                              get_im_ch_ind,                    &! get index in t_tovs%im_ch_v
                              TTOVS_BASE, TTOVS_CI, TTOVS_AV,   &! Constants that determine, which parts of
                              TTOVS_L2C, TTOVS_EMIS,            &! t_tovs are to be stored/loaded
                              TTOVS_FLAG,TTOVS_SPEC,TTOVS_BTCS, &!
                              TTOVS_TR, TTOVS_SINFL,            &!
                              TTOVS_IMCL, TTOVS_IMCH,           &!
                              TTOVS_CLDLEV,                     &!
                              IMCL_FRAC, IMCL_MEAN, IMCL_STDV,  &!
                              IMCH_MEAN_B, IMCH_EMIS,           &!
                              mx_nlev, mx_nav, mx_imch           ! maximum dimension of t_tovs%av
  use mo_ir_emis,       only: n_pc                               ! # of PC to use; -1: off
  use mo_emis,          only: read_nml_emis,                    &
                              update_emis_opt,                  &
                              m_emis,                           &
                              MODE_ATLAS,                       &
                              ATLS_BRDF,                        &
                              dynret_w_pref, dynret_avg
  use mo_tskin,         only: read_nml_tskin
  use mo_satid,         only: satname,                          &! Satellite name
                              satid_bufr2rttov                   ! convert WMO satid do RTTOV sat. identifiers
  use mo_obserr_rad,    only: read_obserr_rad_nml,              &! read namelist /OBSERR_RAD/
                              obs_corr_rad,                     &! full R matrix for different instruments
                              n_obserr_rad                       ! # of entries in obs_corr_rad
  !------------------------
  ! access matrix data type
  !------------------------
  use mo_dec_matrix,    only: t_vector,                         &! vector data type
                              t_vector_segm,                    &! vector segment data type
                              t_matrix,                         &! matrix data type
                              insert,                           &! insert sub-matrix
                              insert_full_block_2,              &
                              construct,                        &! allocate data type components
                              destruct,                         &! deallocate data type components
                              gather,                           &!
                              operator(*),                      &! matrix-vector multiply
                              operator(+),                      &! add vectors
                              assignment(=)                      ! assign matrices or vectors
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,           only: NF90_FILL_INT
  use mo_t_netcdf,      only: stanc,           &! NetCDF start  parameter
                              counc,           &! NetCDF count  parameter
                              strnc,           &! NetCDF stride parameter
                              ncid,           &!
                              get_var           ! read variable
  !-----------------------
  ! background error model
  !-----------------------
  use mo_set_matrix,    only: set_Pb                             !
  use mo_t_bg_err_op,   only: compress_cov,                     &! store covm only once
                              uncompress_cov                     ! store covm on each PE
  use mo_tovs_prof,     only: fill_rad,                         &! fill t_radv
                              prep_rttov_prof,                  &! prepare t_rttov_prof type with profile info
                              t_jac,                            &!
                              t_jac_arr,                        &!
                              destruct,                         &! t_jac, t_jac_arr
                              construct,                        &! t_jac
                              valid,                            &!
                              no_ps_dep,                        &!
!                              no_t_dep,                         &!
                              force_fastem,                     &!
                              surf_class_vers,                  &!
                              check_recalc_stype,               &!
                              scale_qi,                         &!
                              min_dia_qc,                       &!
                              min_dia_qi,                       &!
                              e_bg_ts_sea,                      &!
                              e_bg_ts_land,                     &!
                              e_bg_drts_land,                   &
                              e_bg_ts_ice ,                     &!
                              e_bg_t_top  ,                     &!
                              e_bg_lnq_top,                     &!
                              use_hum_top,                      &!
                              fg_prof_top,                      &!
                              p_top,                            &!
                              warn_gen_hum,                     &!
                              snw_frc_mode,                     &!
                              snw_frc_par,                      &!
                              c_var,                            &!
                              get_tovs_var,                     &!
                              errmsg
  use mo_range_fparse,  only: t_range_fparse,                   &!
                              init,                             &!
                              evaluate,                         &!
                              destruct
  use utilities,        only: sortrx

  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  !----------------------------------------
  ! general operations on TOVS observations
  !----------------------------------------
  public :: process_tovs      ! calculate costfunction, derivatives, ...
  public :: process_tovs_mult ! call rttov_direct / rttov_k in parallel
  public :: read_fdbk_rad     ! read RADIANCES specific stuff from fdbk-file
  public :: write_rttov_prof  ! write profiles on rttov levels
  public :: f_ins             ! 'level' = chan + f_ins * instrument
  public :: rttov_mult_prof   ! RTTOV call with multiple profiles
  public :: valid             ! returns true, if observation for given state
                              ! passed all 3dvar-checks
  public :: lBii              ! flag value to write B in interpolation space

  public :: decr_tovs_use     ! decrease status/set flags for TOVS
  public :: calc_tovs_quality ! calculate pcc values for TOVS
  public :: p_bcast
  public :: destruct_tovs_oe  ! destruct the t_range_fparse structures in rad_set
  public :: thin_superob_tovs ! Superobbing on TOVS observations
  public :: use_reff          ! Use effective cloud droplet radii from ICON
  public :: print_tovs_levs   ! print stratistics on relevant levels
  public :: feedbk_files      ! sat_pp input files
  public :: prep_H_btcs       ! gather clear-sky Tb in t_vector
  !---------------------------------------------------------------------
  ! AMSU-B cloud detection parameters (shall be moved to mo_cloud later)
  !---------------------------------------------------------------------
  public :: param_file_amsub
  public :: amsub_ec_bnd
  public :: amsub_delta_ch20_ch18
  public :: clchk_amsub
  !---------------------------------------------------------------------
  ! estimation of representative pressure level (for LETKF localisation)
  !---------------------------------------------------------------------
  public :: pl_e_bg_t         ! typical temp. background error
  public :: pl_e_bg_rh        ! typical rel.hum. background error
  public :: pl_method         ! method
  public :: pl_log            ! use log(p) for mean, stdev
  !---------------------------------------------------------------------
  ! required for mec
  !---------------------------------------------------------------------
  public :: read_tovs_nml
  public :: filter_ch_stat
  public :: n_filter_ch
  !---------------------
  ! tracegases for RTTOV
  !---------------------
  public :: trg_file         ! File with atmospheric variables from external sources
  public :: trg_clim_file    ! File with atmospheric variables from external sources (climatology)
  public :: trg_clim_trend   ! climate trend for tracegases
  public :: t_trg_clim_trend ! type of trg_clim_trend
  public :: trg_hist_file    ! History file for trace gases
  public :: trg_hist_inidate ! Initialization date for trace gas history file
  public :: glob_use_o3      ! ior of all use_o3 options
  public :: glob_use_co2     ! ior of all use_co2 options
  !---------------------
  ! extrap. above model top
  !---------------------
  public :: nlev_ext         ! number of extra levels
  public :: p_top_ext        ! File with atmospheric variables from external sources (climatology)
  public :: p_blend_ext      ! Top of blending between model and external profile
  public :: nld_t_max        ! max. number of T-dummy levels
  public :: nld_h_max        ! max. number of hum.-dummy levels
  !--------------
  ! interpolation
  !--------------
  public :: iatm_tovs_all    ! all COL_* required for all radiances
  public :: COL_CLD_TOVS     ! all possible COL_* values that might be required for cloudy RTTOV

!------------------------------------------------------------------------------
  !================================
  ! Module variables and data types
  !================================

  !-------------------------
  ! private module variables
  !-------------------------
  real(wp),parameter :: highland_thresh = 1000._wp          ! height parameter [m] defining highland surface type
  integer ,parameter :: f_ins           = 100     ! 'level'=chan+f_ins*instrument
  integer ,parameter :: mf              = m_source! max number of i/outp. files
  integer            :: format_fdbk(mf) = 0      ! format of input/outp. files


  !==================
  ! Namelist TOVS_OBS
  !==================
  !------
  ! Input
  !------
  character(len=128) :: netcdf_path     = ''       ! path to read input files
  character(len=128) :: data_path       = ''       ! path to coefficient files
  character(len=128) :: feedbk_files(mf)= ''       ! names of input/output files
  !---------------
  ! Data selection
  !---------------
  integer            :: flg_sur         =  0       ! 1dvar surface flag
  logical            :: use_ch_blcklst  = .true.   ! use channel blacklist from preprocessing
  logical            :: req_nmlst       = .true.   ! require namelist for every input file
  integer            :: n_filter_ch     = 0
  integer            :: filter_ch_stat(n_stat) = -1
  !-------------------
  ! AMSU-B cloud check
  !-------------------
  character(len=128) :: param_file_amsub = ''      ! AMSU-B clchk parameter file
  character(len=128) :: param_fil_amsub  = ''      ! old, obsolete name (kept for backwards compatibility)
  integer            :: clchk_amsub     = 0        ! AMSU-B cloud check
  real(wp)           :: amsub_ec_bnd (2)           ! ECMWF cloud check bounds
  real(wp)           :: amsub_delta_ch20_ch18      ! AMSU-B ch(20)-ch(18) threshold
  !------------
  ! Computation
  !------------
  integer            :: rttov_version   = -1       ! rttov version to use (7, 10, or 12)
  integer            :: rttov_levels    = 43       ! number of rttov levels used
  logical            :: lhum_dum_ana    = .true.   ! dummy humidity above an.lev
#if (_RTTOV_VERSION >= 12)
  integer            :: rt_ir_emis_mod  = 2        ! IREMIS
#else
  integer            :: rt_ir_emis_mod  = 1        ! ISEM
#endif
  integer            :: rt_mw_emis_mod  = 5        ! FASTEM version
  logical            :: rt_do_lambertian= .false.  ! Lambertian reflection in RTTOV
  logical            :: rt_use_q2m      = .false.  ! Use q2m in in RTTOV
  real(wp)           :: rt_salinity     = 0._wp    ! Salinity for emiss. calc. in RTTOV (PSU)

  integer            :: fix_hgpl        = 0        ! Flag to control bugfix for psurf<->zurf-Bug in RTTOV
  logical            :: app_reg_lims    = .false.  ! Apply (RTTOV) regularization limits
  integer            :: qflag_bitmask   = 0        ! CHK_OPERATOR, if iand(radiance%quality,qflag_bitmask) /= 0
  logical            :: use_reff        = .false.  ! Use effective radii from ICON
  logical            :: rt13_clip_gas_opdep = .true. ! Crucial for convergence in our variational system,
                                                   ! which requires .false.
  integer            :: rt_alloc_mode   = 1        ! (De)allocation of array required by RTTOV
                                                   ! A higher number means: less (de)allocations,
                                                   ! but temporarily more memory required.
                                                   ! 0: (de)allocate for each RTTOV call
                                                   ! 1: (de)allocate for each new instrument
                                                   ! 2: (de)allocate when necessary and at
                                                   !    start/end of TSK_K/Y call
                                                   ! 3: (de)allocate when necessary, i.e. keep arrays
                                                   !    as long as possible

  !---------------------
  ! Tracegases for RTTOV
  !---------------------
  character(len=128) :: trg_file       = ''        ! file with tracegas conc. from external sources
  character(len=128) :: trg_clim_file  = ''        ! file with tracegas conc. from external sources (clim.)
  type t_trg_clim_trend
    character(len=10) :: gas_name      = ''
    character(len=14) :: date          = ''
    real(kind=wp)     :: val           = -1._wp  ! value at date [kg/kg]
    real(kind=wp)     :: trend         =  0._wp  ! annual trend [kg/(kg*year)]
  end type t_trg_clim_trend
  type(t_trg_clim_trend), parameter :: empty_trg_clim_trend = t_trg_clim_trend('','',-1._wp,0._wp)
  type(t_trg_clim_trend), target    :: trg_clim_trend(ngases)
  character(len=128)  :: trg_hist_file = 'bias_RAD_TRGHIST._YYYYMMDDHHMM_'
  character(len=14)   :: trg_hist_inidate = '2099010100'
  integer             :: glob_use_o3      ! ior of all use_o3 options
  integer             :: glob_use_co2     ! ior of all use_co2 options


  !------------------------------------------------------------
  ! pressure level estimation (for LETKF vertical localisation)
  !------------------------------------------------------------
  real(wp)          :: pl_e_bg_t       = 0.5_wp ! nominal temp. background error
  real(wp)          :: pl_e_bg_rh      = 0.1_wp ! nominal rel.hum. background error
  logical           :: pl_log          = .true. ! use log(p)
  integer           :: pl_method       = 2      ! method :
  !                    In nominal_height in mo_obs:
  !                    1: take maximum value of temperature Jacobian
  !                    2: take mean value of Jacobian with respect to
  !                       temperature and relative humidity,
  !                       weighted by the nominal background errors
  !                    3: as 2 but use sensitivities squared
  !                       instead of absolute values
  !                    4: as 2 but use actual background errors
  !                    In get_hgt@process_tovs_mult (because the transmission is required
  !                    5: take the mean value and std.dev. of transm * opdep
  !                    6: take the mean value and std.dev. of the weighting fnct., i.e.
  !                       the derivative of the transmission
  !                    7: take the height, where the transmission is 0.5,
  !                       as width use half the vertical distance between the levels,
  !                       where the transmission is pl_wd_thresh and 1.-pl_wd_thresh
  real(kind=wp)     :: pl_vis          = -1._wp ! plevel for vis radiances
  real(kind=wp)     :: pl_wd_thresh    = 0.2_wp !
  real(kind=wp)     :: plev_min        = 0.0_wp
  !-----------------------------------------------
  ! multiple RTTOV profiles processed in one chunk
  !-----------------------------------------------
  logical           :: rttov_mult_prof = .true.   ! multiple prof. RTTOV call
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! temporary work arounds for discontinuities in RTTOV10
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(wp)          :: rt_min_od       = 1.e-5_wp ! parameter min_od in RTTOV10
  real(wp)          :: rtt10_min_od    = 1.e-5_wp ! obsolete, use variable above
  !------------
  ! Diagnostics
  !------------
  integer           :: monitor_prof          = 0                ! monitor profiles (fg-scan)
  integer           :: mon_ana_prof          = 0                ! monitor profiles (an-scan)
  integer           :: monitor_prof_split    = 0                ! split monRTOVP files
  integer           :: mon_ana_prof_split    = 0                ! split cofRTOVP files
                                                                !   1: separate file for each satellite
                                                                !   2: separate file for each instrument grid
  real(wp)          :: monitor_prof_wrbxf    = 1.0_wp           ! Fraction of boxes, that are written at once
  real(wp)          :: mon_ana_prof_wrbxf    = 1.0_wp           ! Fraction of boxes, that are written at once
                                                                ! Values < 1.0 reduce the required memory but
                                                                ! slow down the program
  integer           :: monitor_prof_sat(m_rad_set)  = -1        ! satellites in monRTOVP  (-1: all satellites)
  integer           :: mon_ana_prof_sat(m_rad_set)  = -1        ! satellites in cofRTOVP  (-1: all satellites)
  integer           :: monitor_prof_grid(m_rad_set) = -1        ! instr.grids in monRTOVP (-1: all grids)
  integer           :: mon_ana_prof_grid(m_rad_set) = -1        ! instr.grids in cofRTOVP (-1: all grids)

  integer           :: hd_id_vers            = 2                ! Version of spot%hd%id
  logical           :: lBii                  = .false.          ! write Bii
  logical           :: trace_mem_time        = .false.          ! print memory & time usage
  integer           :: flag_instr            = -1               ! what to do with flag mismatch within instr.
  !-----------------------------------
  ! t_rad_gopt and t_rad_iopt defaults
  !-----------------------------------
  integer           :: quality_mode        = 0            ! mode for tovs quality calculation
  integer           :: l2c_type            = 2            ! Type of level to channel assignment
  logical           :: l2c_god_corr        = .false.      ! god-corrected l2c
  real(wp)          :: l2c_rel_lim         = 0.01_wp      ! relative limit used in level to channel assignment
  real(wp)          :: l2c_abs_lim         = -1._wp       ! absolute limit used in level to channel assignment
  logical           :: l2c_use_rad         = .true.       ! Use radiances or BTs in level to channel assignment
  real(wp)          :: l2c_max             = -1._wp       ! maximum allowed level in l2c check (default)
  real(wp)          :: l2c_max_trop        = -1._wp       ! maximum allowed level in l2c check
  real(wp)          :: l2c_max_midlat      = -1._wp       ! maximum allowed level in l2c check
  real(wp)          :: l2c_max_polar       = -1._wp       ! maximum allowed level in l2c check
  real(wp)          :: l2c_max_tr2ml       = 30._wp       ! border between tropics and midlatitudes for l2c check
  real(wp)          :: l2c_max_tr2ml_width = 10._wp       ! border width between tropics and midlatitudes for l2c check
  real(wp)          :: l2c_max_ml2pl       = 60._wp       ! border between midlatitudes and polar for l2c check
  real(wp)          :: l2c_max_ml2pl_width = 10._wp       ! border width between midlatitudes and polar for l2c check
  integer           :: surf_infl_mode      = 0            ! surface influence test mode
  real(wp)          :: max_surf_infl       = 0.1_wp       ! Maximum influence by surface
  real(wp)          :: d_stemp             = 3._wp        ! delta(t_skin)
  real(wp)          :: d_emiss             = 0.1_wp       ! delta(emissivity)
  logical           :: l_max_surf_hgt      = .false.      ! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT

  !------------
  ! Debugging
  !------------
  integer           :: spt_debug_prof   =  0        ! write profiles in hdf5
  integer           :: debug_opdep_chan = -1        ! 0: all chans, >0: selected channel
  !------------
  ! Obsolete
  !------------
  integer            :: max_scan        = 1000     ! max. number of scans
  integer            :: flg_prc         = 15       ! 1=data,2=min,4=sur,8=cld
  integer            :: flg_prc_amsua   =  0       ! 1=data,2=min,4=sur,8=cld
  integer            :: flg_prc_amsub   =  0       ! 1=data,2=min,4=sur,8=cld
  integer            :: flg_prc_hirs    =  0       ! 1=data,2=min,3=sur,8=cld

  !----------------------------------------
  ! range of bound constrained optimisation
  !----------------------------------------
  real(wp) ,parameter :: ctpmax         = 1050._wp ! maximum cloud top pressure
  real(wp) ,parameter :: ctpmax_r       =   50._wp ! nonlinear range at maximum
  real(wp) ,parameter :: ctpmin_r       =   50._wp ! nonlinear range at minimum
  real(wp) ,parameter :: ctpmin         =   50._wp ! minimum cloud top pressure
  real(wp) ,parameter :: clfmax         =  1.00_wp ! maximum cloud fraction
  real(wp) ,parameter :: clfmax_r       =  0.25_wp ! nonlinear range at maximum
  real(wp) ,parameter :: clfmin_r       =  0.25_wp ! nonlinear range at minimum
  real(wp) ,parameter :: clfmin         =  0.00_wp ! minimum cloud fraction
  !---------------------------------------------
  ! cloud top and cloud fraction
  ! and derived transformed dummy sink variables
  !---------------------------------------------
  real(wp)          :: cld_top_def   = 500._wp  ! first guess cloud top
  real(wp)          :: cld_top_def_v = 500._wp  ! transformed sink variable
  real(wp)          :: cld_frc_def   =   0._wp  ! first guess cloud fraction
  real(wp)          :: cld_frc_def_v = -.25_wp  ! first guess cloud fraction
  real(wp)          :: cld_top_e     =   0._wp  ! cloud top  background error
  real(wp)          :: cld_frc_e     =   0._wp  ! cloud frc. background error
  real(wp)          :: cld_top_e_v   =   0._wp  ! cloud top  background error
  real(wp)          :: cld_frc_e_v   =   0._wp  ! cloud frc. background error
  integer           :: bound_constr  =   1      ! cld_.. boundary constraints
                                                !        0: no
                                                !        1: logarithm.transform
                                                !        2: quadratic transform
  !---------------------------------------------
  ! surface influence test
  !---------------------------------------------
  integer, parameter :: SINFL_ABS      = 0        ! Use abs influence on T_b
  integer, parameter :: SINFL_REL      = 1        ! Use relative influence on T_b
  integer, parameter :: SINFL_LAND     = 2        ! Apply test only over land
  integer, parameter :: SINFL_HIGHLAND = 3        ! Apply test only over highland
  integer, parameter :: SINFL_TRANSM   = 4        ! Calc. infl. on basis of transmittances


  !-----------------------------------------------------
  ! possible values for flags monitor_prof, mon_ana_prof
  !-----------------------------------------------------
  integer, parameter :: MP_OLD      =  1 ! use write_rttov_prof instead of write_rttov_prof2
  integer, parameter :: MP_NETCDF   =  2 ! write NetCDF file
  integer, parameter :: MP_PSAS_ANA =  4 ! write info on PSAS analysis
  integer, parameter :: MP_H        =  8 ! write H matrices
  integer, parameter :: MP_B        = 16 ! write B in interpolation space
  integer, parameter :: MP_NOSORT   = 32 ! skip ordering in *RTOVP.nc
  integer, parameter :: MP_HBH      = 64 ! write covariances in obs-space

!   type t_jac_arr
!     type(t_jac), pointer  :: a(:) => NULL()
!   end type t_jac_arr

  type t_keep_rtstat_spt
    integer          :: id         =  -1
    integer          :: n          =  0 ! Number of integers required to keep information
    integer          :: nbits(0:4) = 0
    integer, pointer :: i(:)       => null()
  end type t_keep_rtstat_spt
  type t_keep_rtstat
    integer                          :: nspot  =  0
    type(t_keep_rtstat_spt), pointer :: spt(:) => null()
  end type t_keep_rtstat
  type(t_keep_rtstat), save, pointer :: keep_rtstat(:) => null()
  integer,  parameter :: nbits_int = bit_size(0)
  type t_rtstat_opts
    integer          :: action_fg  = 0     ! action in first guess checks if test is not passed:
                                           !   Bit 0: print, Bit 1: set flag/op_na for obs, Bit 2: set flag/op_na for report
                                           !                             + check specific bits
    integer          :: action_ana = 0     ! action in analyses if test is not passed
    integer          :: action(50) = 0     ! action in minimization if test is not passed
    integer          :: action_mon = 0
    integer          :: op_na_bit  = -1    ! bit set in op_na
    integer          :: tsk        = TSK_Y ! Task(s) for which test is active
    integer          :: nitout     = 0     ! Minimum number of outer loops with RTTOV problem before op_na is set.
  end type t_rtstat_opts
  type(t_rtstat_opts), parameter    :: empty_rtstat_opts = t_rtstat_opts(0,0,0,0,-1,TSK_Y,0)
  type(t_rtstat_opts), save, target :: chk_reg_lims
                              ! action bit4: distinguish vars., bit5: distinguish levels, bit6: distinguish upper/lower limits
  type(t_rtstat_opts), save, target :: chk_opdep
!  real(wp), parameter :: opdep_thresh = 1.E-03
  real(wp)            :: opdep_thresh(50) = 1._wp  ! opdep is NEGATIVE
  real(wp)            :: opdep_min_transm(50) = 0._wp

  type(t_rtstat_opts), save, target :: chk_ps

  type(t_rtstat_opts), save, target :: chk_god
  real(wp)            :: god_thresh(50) = 99._wp

  integer, save :: wr_opdep = 0            ! Write optical depths from RTTOV for humidity sensitive channels
                                           ! bit 0 : first guess
                                           ! bit 1 : analyses
                                           ! bit 2 : tsk_k in minimisation
                                           ! bit 3 : tsk_y in minimisation (i.e. calls from line search)
                                           ! bit 4 : monitoring
                                           ! bit 5 : output also for passive channels
                                           ! bit 6 : output only for chans with discontinuity
                                           !         (only with bit 4 and chk_opdep active in minimization)
  logical, save :: wr_opd_disc       = .false.
#if (_RTTOV_VERSION >= 12)
  type t_wr_opd
    integer  :: chan
    integer  :: lev
    integer  :: spt
    real(wp) :: opdep   ! raw optical depth
    real(wp) :: opdep_  ! effective optical depth
    real(wp) :: transm(2) ! transmittances from bottom and top of layer to space
    real(wp) :: q
    real(wp) :: rh
    real(wp) :: t
  end type t_wr_opd
#endif
  type t_keep_disc
    integer :: iset
    integer :: id
    integer :: chan
    integer :: lev
  end type t_keep_disc
  type(t_keep_disc), pointer   :: keep_disc(:) => null()
  integer                      :: n_keep_disc  =  0
  integer,           parameter :: n_alloc_keep_disc = 10

  ! Obserror calculation
  integer,  save              :: m_var = 10
  real(wp), save, allocatable :: v(:)        ! variable values required for obserror calc.

  ! Dummy levels above model top
  integer                     :: nlev_ext    = 10    ! extra levels above model top for lev_mode > 0
  real(kind=wp)               :: p_top_ext   = 0._wp ! highest extra pressure level
  real(kind=wp)               :: p_blend_ext = huge(0._wp) ! top of blending model and file_strato
  integer                     :: nld_t_max   = 0     ! maximum value of t_tovs%nld_t
  integer                     :: nld_h_max   = 0     ! maximum value of t_tovs%nld_h

  ! Interpolation to "interpolation space" (interpolate@mo_psasutil)
  integer(i8)                 :: iatm_tovs_all    = 0
  integer(i8), parameter      :: COL_CLD_TOVS     = COL_CLC + COL_QCDIA + COL_QIDIA + &
                                                    COL_REFF_QC + COL_REFF_QI

  namelist /TOVS_OBS/ max_scan, netcdf_path, data_path, read1pe,         &
                      feedbk_files, rttov_version,                       &
                      e_bg_ts_sea, e_bg_ts_land, e_bg_ts_ice,            &
                      e_bg_t_top,e_bg_lnq_top, p_top,                    &
                      flg_prc, flg_prc_amsua, flg_prc_amsub,             &
                      flg_prc_hirs, flg_sur,                             &
                      clchk_amsub, param_fil_amsub, amsub_ec_bnd,        & ! obsolete: should be set in TOVS_CLOUD
                      amsub_delta_ch20_ch18,lhum_dum_ana,                &
                      monitor_prof,       mon_ana_prof,                  &
                      monitor_prof_split, mon_ana_prof_split,            &
                      monitor_prof_sat,   mon_ana_prof_sat,              &
                      monitor_prof_grid,  mon_ana_prof_grid,             &
                      monitor_prof_wrbxf, mon_ana_prof_wrbxf,            &
                      hd_id_vers, rttov_mult_prof, req_nmlst,            &
                      use_ch_blcklst, rttov_levels,                      &
                      l2c_type, l2c_god_corr, no_ps_dep,                 & !no_t_dep,
                      rt_min_od, rtt10_min_od, trace_mem_time,           &
                      cld_top_def, cld_frc_def, cld_top_e, cld_frc_e,    &
                      max_surf_infl, d_stemp, d_emiss, surf_infl_mode,   &
                      spt_debug, spt_hd_debug, usd, spt_debug_prof,      &
                      bound_constr,                                      &
                      n_max_calc_k_def, n_max_prof_k_def,                &
                      n_max_calc_y_def, n_max_prof_y_def,                &
                      quality_mode, l2c_max, l2c_max_trop,               &
                      l2c_max_midlat, l2c_max_polar, l2c_rel_lim,        &
                      l2c_abs_lim, l2c_use_rad, l2c_max_tr2ml,           &
                      l2c_max_tr2ml_width, l2c_max_ml2pl,                &
                      l2c_max_ml2pl_width, l_max_surf_hgt,               &
                      qflag_bitmask, pl_e_bg_t, pl_e_bg_rh, pl_method,   &
                      pl_log, pl_wd_thresh, pl_vis,plev_min,             &
                      rt_ir_emis_mod, rt_mw_emis_mod, force_fastem,      &
                      rt_humi, rt_do_lambertian, rt_use_q2m,             &
                      rt_salinity, rt_levs_exact, rt_hard_limits,        &
                      chk_reg_lims, chk_opdep, chk_ps,                   &
                      opdep_thresh, opdep_min_transm, god_thresh,        &
                      app_reg_lims, fix_hgpl, snw_frc_mode, snw_frc_par, &
                      god_par_file, wr_god, wr_opdep, chk_god,           &
                      use_hum_top, fg_prof_top, warn_gen_hum, filter_ch_stat,&
                      surf_class_vers, trg_file, trg_clim_file,          &
                      trg_clim_trend, trg_hist_file, trg_hist_inidate,   &
                      use_reff, min_dia_qc, min_dia_qi, scale_qi,        &
                      rt13_clip_gas_opdep, atlas_single_inst,            &
                      nlev_ext, p_top_ext, p_blend_ext, rt_alloc_mode,   &
                      dynret_w_pref, dynret_avg, check_recalc_stype,     &
                      debug_opdep_chan, e_bg_drts_land

  !-----------
  ! interfaces
  !-----------
  interface load
    module procedure load_rtstat
  end interface load

  interface store
    module procedure store_rtstat
  end interface store

  interface destruct
    module procedure destruct_rtstat
    module procedure destruct_rtstat_spt
  end interface destruct


  interface p_bcast
    module procedure bcast_rtstat_opts
    module procedure bcast_trg_clim_trend
  end interface

#if (_RTTOV_VERSION >= 12)
  interface p_alltoall
    module procedure p_alltoall_t_wr_opd
  end interface
#endif

!==============================================================================
contains
!==============================================================================
  subroutine process_tovs (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
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
    character(len=*), parameter   :: proc = 'process_tovs'
    integer                       :: tsk            ! task (local copy)
    integer                       :: i,j,k,l,m      ! indices
    integer                       :: ii,ll          ! input data index offset
    integer                       :: jj, kk
    integer                       :: k0, i0, i1, ll0, ll1
    integer                       :: n, no
    integer                       :: iset           ! rad_set entry that is linked to rad. dataset
    character(len=128)            :: file
    character(len=1000)           :: msg
    integer                       :: ncv            ! size of interpolation space
    integer                       :: npv            ! number of profile variables
    integer                       :: stat
    type(t_radv)                  :: radv
    type (t_rad_set)              :: rad_set_tmp(m_rad_set)
    type (t_rad_set), pointer     :: rs => null()
    type (t_tovs)                 :: ttovs          ! tovs data type
    integer,          allocatable :: ci(:)          ! t_tovs channel indices
    type (t_tovs_instr)           :: ti             ! info on instruments in ttovs
    type (t_set)                  :: set            ! set of characteristics for channels
    integer                       :: iii            ! work around NEC SX bug
    logical                       :: change         ! argument to shrink_report
    integer(i8)                   :: iatm           ! parameters required as input
    logical,          pointer     :: msk(:)=>NULL() ! argument to shrink_report
    logical                       :: lvalid         ! Whether a read file contains valid data
    logical                       :: lproc
    logical                       :: ld
    integer                       :: ierr
    integer                       :: recs, recl, rec1, records
    integer                       :: subset1, subset0, subsets, subsetl
                                                    ! instruments in t_tovs
    integer                       :: nlev           ! total number of levels (with extra levels)
    integer                       :: nlev_mod       ! number of model levels
    integer,          allocatable :: rad_set_pe(:)  ! pe that modified the rad_set entry.
    real(wp),         allocatable :: Hnew(:,:)
    integer,          allocatable :: ix(:)          ! channel indices
    real(wp),         allocatable :: ee(:)          ! errors ( sqrt(variances))
    real(sp),         allocatable :: var(:)         ! errors (      variances )

    ! required for (single profile) RTTOV calls
    type(t_jac)                   :: tjac
    type(t_rttov_prof)            :: rtp            ! variables passed         to RTTOV
    type(t_rttov_prof)            :: rtp_k          ! Jacobi matrix provided   by RTTOV
    integer                       :: ifail
    real(wp),         allocatable :: tb(:)

    ! Obserror calc.
    type(t_range_fparse), pointer :: trf => null()
    integer                       :: instr

    ! Setup of av t_tovs%av:
    logical                       :: setup_tt, setup_tt_tr
    logical                       :: l_o3, l_co2, l_cloud, l_vint_rt, l_t, l_q
    integer                       :: nld_t, nld_h

    ! Pre-interpolation of P
    type(t_col),      pointer     :: c
    real(wp),         allocatable :: p(:)           ! pressure profile
    real(wp)                      :: p_top_mod      ! model top pressure
    real(wp)                      :: p_top_         ! higest pressure for T-dummy variables
    real(wp)                      :: p1, p0, dlnp
    real(wp)                      :: wgt(npmax), wt     ! interpolation weights
    integer                       :: it, ic         ! loop indices
    integer                       :: nl_int         ! number of levels to be interpolated
    integer                       :: nt, ncol, ml(1)

    !======================
    ! executable statements
    !======================
!   if (rept_use(OT_RAD)% use(CHK_NONE) <= STAT_DISMISS) return

    if (rept_use(OT_RAD)% init == 0) return
    tsk = task
    if(tsk==0) return


    if (rttov_mult_prof) then
      if (task==TSK_Y) return
      if (task==TSK_K) then
        !-------------------------------------------------------
        ! Reserve space in H matrix for results calculated later
        ! in process_tovs_mult
        ! TODO this can be done more efficiently
        !-------------------------------------------------------
        allocate (Hnew(spot% o% n, spot% i% n))
        Hnew = 1._wp
        call insert (obs% H, Hnew, spot% o% i, spot% i% i)
        return
      end if
    end if

    !======================
    ! module initialization
    !======================
    if (iand (TSK_INIT,tsk) /= 0) then
      !------------------------------------
      ! read namelist, set module variables
      !------------------------------------
      call read_tovs_nml        ! read namelists /TOVS_OBS/ and /TOVS_OBS_CHAN/
      if (rttov_version > 0) then
        call read_obserr_rad_nml   ! read namelist /OBSERR_RAD/
        call scan_satpp_feedback ()
        call set_size
      end if

      !-------
      ! return
      !-------
      tsk = tsk - TSK_INIT
      if (tsk == 0) return
    endif
    !==================
    ! read NetCDF files
    !==================
    if (iand (TSK_READ,tsk) /= 0) then
      !---------------------
      ! new input file style
      !---------------------
FTRACE_BEGIN("process_tovs:read")

      allocate(rad_set_pe(size(rad_set)))
      rad_set_pe = -1

      call flush_buf
      call add_line_pio(repeat('-',79))
      do i=1, n_source
        if (source(i)% obstype  /= OT_RAD)   cycle
        if (source(i)% filetype /= FT_SATPP) cycle
        file = path_file (source(i)% path, source(i)% file)

        lproc = source(i)% used
        subset0 = 0
        if (i > 1) subset0 = bufr_inv(OT_RAD)% subseto(i-1)
        subsetl = bufr_inv(OT_RAD)% subsetl
        subsets = bufr_inv(OT_RAD)% subsets
        subset1 = bufr_inv(OT_RAD)% subseto(i)
        lproc   = lproc .and. (subsetl >  subset0)
        lproc   = lproc .and. (subsets <= subset1)
        recs    = max (subsets-subset0, 0)
        recl    = max (subsetl-subset0, 0)
        records = min(rept_use(OT_RAD)% max_proc, recl - recs)
        rec1    = recs + 1
        lproc   = lproc .and. (records-1 > 0)
        if (lproc) then
          !-------------------
          ! read feedback file
          !-------------------
          FTRACE_BEGIN("process_tovs:read_satpp")
          call read_satpp_feedbk (file, radv, lvalid, status=ierr,        &
            lprint=.false., lread=.true., istart=rec1, iend=rec1+records-1, pe=dace% pe)
          if (ierr /= 0) call finish(proc//'(TSK_READ)', 'Failed to read radiance dataset: '//trim(err_msg))
          FTRACE_END("process_tovs:read_satpp")
          if (lvalid) then
            radv% file_id = i
            !----------------
            ! Print some info
            !----------------
            write (msg,'(4x,"sat. ",I3.3,", grid ",I3.3,", pe ",I3.3,":")', iostat=stat)&
                 radv%i%satid, radv%i%grid, dace% pe
            call add_line(msg)
            write (msg,'(6x,"read  : ",I7," FOVs from ",A,1x,I2.2,2(1x,I8))', iostat=stat) &
                 &radv%n_rec, trim(file), radv% file_id,rec1,rec1+records-1
            call add_line(msg)
            !--------------------------------------------
            ! unify t_rad_set in radv with global rad_set
            !--------------------------------------------
            FTRACE_BEGIN("process_tovs:link_rad")
            call link_rad (radv, status=ierr, i_set=iset)
            if ((ierr /= 0) .or. (iset <=0)) call finish(proc//'(TSK_READ)', &
                 &'Failed to link radiance dataset with TOVS_OBS_CHAN namelists.')
            if (rad_set(iset)% id > 0) rad_set_pe(iset) = dace% pe
            FTRACE_END("process_tovs:link_rad")
            !-------------------------------
            ! Set default emissivity (mw)
            !-------------------------------
            allocate(radv% emiss(radv% i% n_chan, radv%n_rec))
            radv% emiss = 0. ! FASTEM
            !-------------------------------
            ! now fill observation data type
            !-------------------------------
            FTRACE_BEGIN("process_tovs:radv_obs")
            call radv_obs (radv, obs% o, sum(source(1:i-1)% entries)+recs)
            FTRACE_END("process_tovs:radv_obs")
          else
            !-------------------------------------
            ! print warning for invalid satpp-file
            !-------------------------------------
            write (6,'(a)') repeat('*',79)
            write (0,'(a)') repeat('*',79)
            write (6,'(/a,a/)') 'WARNING, invalid satpp-file :',trim(file)
            write (0,'(/a,a/)') 'WARNING, invalid satpp-file :',trim(file)
            write (6,'(a)') repeat('*',79)
            write (0,'(a)') repeat('*',79)
          end if
          call destruct(radv)
        end if
        !-------------------------------
        ! write pending output to stdout
        !-------------------------------
        call write_pending
      end do
      call flush_buf
FTRACE_END("process_tovs:read")

FTRACE_BEGIN("process_tovs:rad_set")
      !------------------------------------------------
      ! Send the modified rad_set entries to other pe's
      !------------------------------------------------
      do i = 1, size(rad_set_pe)
        rad_set_pe(i) = p_max(rad_set_pe(i))
      end do
      do i = 1, size(rad_set_pe)
        if (rad_set_pe(i) >= 0) then
          call p_bcast(rad_set(i), rad_set_pe(i))
        end if
      end do
      !--------------------------
      ! Reorder the rad_set array
      !--------------------------
      n_set = count(rad_set(:)%id > 0)
      if (n_set > 0 .and. n_set < size(rad_set)) then
        rad_set_tmp(1:n_set)  = pack(rad_set(:), rad_set(:)%id >  0) ! Used datasets
        rad_set_tmp(n_set+1:) = pack(rad_set(:), rad_set(:)%id <= 0) ! Not used datasets
        rad_set(:)            = rad_set_tmp(:)
        call destruct(rad_set_tmp)
      end if
FTRACE_END("process_tovs:rad_set")

      call update_emis_opt

      tsk = tsk - TSK_READ
      if (tsk == 0) return
    endif

    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then

      spot% int_type  = ITY_ICOL
      spot% cost      = spot% col% nlev * 5
      spot% char      = CHR_NONL+CHR_EXP
      call set_nr     ! set upper bound for R-matrix size for the current spot
      tsk = tsk - TSK_SET_CHR
      if (tsk == 0) return
    endif

    !=================================================
    ! tsk == TSK_SHRINK:
    ! release unused channels in the assimilation step
    !=================================================
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change, mask=msk)
      if (change) then
        call load  (obs% o, spot, ttovs, rs=rs, ti=ti)
        n            = count(msk)
        no           = ttovs% nchan
        ttovs% nchan = n
        ttovs% ci  (1:n) = pack(ttovs% ci  (1:no), mask=msk)
        ttovs% emis(1:n) = pack(ttovs% emis(1:no), mask=msk)
        if (associated(ttovs%cemi)) &
             ttovs% cemi(1:n) = pack(ttovs% cemi(1:no), mask=msk)
        if (associated(ttovs%sinfl)) &
             ttovs% sinfl(1:n) = pack(ttovs% sinfl(1:no), mask=msk)
        if (associated(ttovs%flag)) &
             ttovs% flag(1:n) = pack(ttovs% flag(1:no), mask=msk)
        if (associated(ttovs%bt_cs)) &
             ttovs% bt_cs(1:n) = pack(ttovs% bt_cs(1:no), mask=msk)
        if (associated(ttovs%tr)) then
          do i = 1, ttovs%ntr
            ttovs% tr(1:n,i) = pack(ttovs% tr(1:no,i), mask=msk)
          end do
        end if
        ! l2c is stored only for instruments with l2c_type > 0
        i = 0
        j = 0
        do ii = 1, ti%n_instr          ! Add to/blend with interpolated model profile

          if (rs% iopts(ti%ii(ii))% l2c_type > 0) then
            i0 = ti%o_ch_i(ii)+1
            i1 = ti%o_ch_i(ii)+ti%n_ch_i(ii)
            m  = count(msk(i0:i1))
            ttovs% l2c(i+1:i+m) = pack(ttovs%l2c (j+1:j+ti%n_ch_i(ii)), &
                                  mask=msk(i0:i1))
            j = j + ti%n_ch_i(ii)
            i = i + m
          end if
        end do
        ttovs% nl2c = i
        call store (obs% o, spot, ttovs)
        call destruct(ttovs)
        deallocate (msk)
        call set_nr(ti) ! set upper bound for R-matrix size for the current spot
      else
        call set_nr     ! set upper bound for R-matrix size for the current spot
      endif
      tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !=================================
    ! determine model columns required
    !=================================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      ! Default, required for all instruments and options
      iatm = COL_P + COL_TV + COL_RH + COL_Q

      call load(obs% o, spot, tovs=ttovs, tovs_io=TTOVS_BASE, rs=rs)
      if (any(rs%iopts(1:rs%n_instr)%cloud_mode > 0) .and. .not.associated(ttovs%av)) then
        ! Additional fields for  cloudy RTTOV calls
        iatm =  iatm + COL_QCDIA
        if (.not.mw_instr(rs%instr(1))) then  !Temp. hack!
          iatm =  iatm + COL_CLC + COL_QIDIA  !Temp. hack!
        end if                                !Temp. hack!
        if ( use_reff ) then
           iatm = iatm + COL_REFF_QC + COL_REFF_QI
        end if
      end if
      if (ldeb(spot)) write(usd,*) dpref//' setup_cols iatm=',iatm
      call destruct(ttovs)
      iatm_tovs_all = ior(iatm_tovs_all, iatm)


      call idx_init (      &
            spot% col% c,  &! <-  column descriptor
            spot% col% h,  &!  -> interpolation coefficients
            obs% o% mc,    &! <-> model column descriptors
            iatm,          &! <-  fields required
            0,             &! <-  tracers required
            atm% grid,     &! <-  model grid
            spot% i_time,  &! <-  time slot
            spot% w_time   )! <-  time interpolation weight

      if (ldeb(spot)) then
        write(usd,*) dpref,'idx_init',size(spot%col%h%imc)
        do i = 1, size(spot%col%h%imc,1)
          ii = spot%col%h%imc(i,1)
          write(usd,*) dpref,'ii',i,ii
          if (ii == 0) exit
          write(usd,*) dpref,'w',spot%col%h%w(i),spot%col%h%ijdp
        end do
      end if


      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN, STAT_DISMISS)

      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif

    !================================
    ! setup description of PSAS-space
    !================================
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then

      !----------------------------
      ! request interpolation space
      !----------------------------
      call load  (obs% o, spot, tovs=ttovs, tovs_io=TTOVS_BASE, rs=rs)
!      if (any(rs%iopts(1:rs%n_instr)%cloud_mode >= 2)) then
      if (rs%gopts%lev_mode > 0) then
        call destruct(ttovs)
        call load  (obs% o, spot, tovs=ttovs)
        ttovs%nlev = atm%grid%nz
        if (interp_strato > 0) ttovs%nlev = ttovs%nlev + nlev_ext
        call store (obs% o, spot, tovs=ttovs)
      end if
      nlev = ttovs%nlev
      call destruct(ttovs)

      rs% n_lev = nlev

      npv  = 2*nlev
      ncv  = 2*nlev + nsv
      !-----------------------------------------------
      ! check for emissivity principal component model
      !-----------------------------------------------
      if ( (any(obs% o% body (spot%o%i+1:spot%o%i+spot%o%n)% lev_sig == 221) .or. &
            any(obs% o% body (spot%o%i+1:spot%o%i+spot%o%n)% lev_sig == 620) ) .and. &
           n_pc > 0) then !
        ncv = ncv + n_pc
      endif
      call new_int (obs% o, spot, ncv)
      tsk = tsk - TSK_SETUP_FUL0
      if (tsk == 0) return
    endif

    !================================
    ! setup description of PSAS-space
    !================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then

      !-------------------------------------------------
      ! set observation id.s and interpolation levels
      !
      ! interpolation space layout:
      !
      !   1             profiles of temperature
      !   ...                   and humidity
      !   npv =2*jplev
      !-------------------------------------------------
      !   npv + 1       t_skin      dummy sink variable
      !   npv + 2       surface geopotential height
      !   npv + 3       cloud top   dummy sink variable
      !   npv + 4       cloud frac. dummy sink variable
      !   npv + 5       snow fraction     sink variable
      !-------------------------------------------------
      !   npv + 6 ...   surf.emmis. dummy sink variables
      !   npv + npc     (principal component model)
      !-------------------------------------------------
      call load  (obs% o, spot, tovs=ttovs, rs=rs, tovs_io=TTOVS_BASE)
      l_vint_rt = (rs%gopts%lev_mode <= 0)
      nlev = ttovs%nlev
      if (nlev <= 0) call finish(proc, 'ttovs%nlev <= 0 !! Maybe, TOVS_OBS namelist is missing?')
      if (.not.l_vint_rt .and. interp_strato > 0) then
        nlev_mod = nlev - nlev_ext
      else
        nlev_mod = nlev
      end if
      allocate(p(nlev))

      setup_tt = (ttovs%nld_t < 0) .or. (ttovs%nld_h < 0) .or. .not.associated(ttovs%av)

      setup_tt_tr = (rs%bc%n_tr > 0) .and. .not.associated(ttovs%tr)

      if (setup_tt) then
        call destruct(ttovs)
        call load  (obs% o, spot, tovs=ttovs)
        ! Determine p-profile
        ! Pre-interpolation of P
        if (.not.l_vint_rt .or. interp_strato > 0) then
          if (.not.present(cols)) call finish(proc//'(TSK_SETUP_FULL)', &
               'cols not present in TSK_SETUP_FULL')
          ncol = size(spot% col% h% imc,1)
          wgt(1:ncol) = spot% col% h% w(1:ncol)
          if ( int_nn .and. spot% col% h% imc(1,1) /= 0) then
            ml = maxloc(wgt(1:ncol))
            wgt        = 0._wp
            wgt(ml(1)) = 1._wp
          endif
          nt = 2; if (spot% w_time==0) nt = 1
          if (.not.l_vint_rt) then
            nl_int = nlev_mod
          else
            nl_int = 1    ! Only model top required
          end if
          p = 0._wp
          do it = 1, nt
            do ic = 1, ncol
              ii = spot% col% h% imc(ic,it)
              if (ii==0) exit
              wt = spot% w_time
              if (it==1) wt = 1._wp - wt
              c => cols% col(ii)
              do i = 1, nl_int
                p(i) = p(i) + wt * wgt(ic) * c% p (i)
              enddo
            end do
          end do
          p(1:nl_int) = exp(p(1:nl_int))
          p_top_mod = p(1)
          if (ldeb(spot)) write(usd,*) dpref,'setup_full p',p_top_mod,p
        end if

        if (.not.l_vint_rt) then
          if (interp_strato > 0) then
            ! Add levels to be filled from interpolate_strato
            p(nlev_ext+1:nlev_ext+nlev_mod) = p(1:nlev_mod)
            p1 = log(p_top_mod)
            p0 = log(p_top_ext)
            if (p0 > p1) p0 = log(0.9_wp*p_top_mod)
            dlnp = (p1 - p0) / nlev_ext
            do i = 1, nlev_ext
              p(i) = exp(p0 + (i-1)*dlnp)
            end do
          end if
        else
          p = preslev
        end if

        if (interp_strato > 0) then
          ttovs%nl_st = count(p < p_top_mod)
          p_top_ = max(p_top, p_top_mod)
        else
          ttovs%nl_st = 0
          p_top_ = p_top
        end if
        nld_t = count(p < p_top_)
        nld_h = count(p < hum_ana_top)
        nld_h = max(nld_h, nld_t)
        if (.not.lhum_dum_ana) nld_h = nld_t
        nld_t_max = max(nld_t_max, nld_t)
        nld_h_max = max(nld_h_max, nld_h)
        ttovs%nld_t = nld_t
        ttovs%nld_h = nld_h
      else
        nld_t = ttovs%nld_t
        nld_h = ttovs%nld_h
        if (.not.l_vint_rt) then
          call destruct(ttovs)
          call load  (obs% o, spot, tovs=ttovs, tovs_io=TTOVS_BASE+TTOVS_AV)
          p = ttovs%av(1:nlev,ttovs%i_p)
          call destruct(ttovs)
        else
          p = preslev
        end if
      end if

      npv  = 2*nlev         ! number of profile                variables
      ncv  = spot%i%n       ! total number of                  variables
      ii   = spot%i%i       ! input data index offset
      obs% o% t_int (ii+1        : ii+2*nld_t  :2) = OBS_DUM
      obs% o% t_int (ii+1+2*nld_t: ii+npv      :2) = OBS_TV
      obs% o% t_int (ii+2        : ii+2*nld_h+1:2) = OBS_DUM
      obs% o% t_int (ii+2+2*nld_h: ii+npv      :2) = OBS_RH
      obs% o% t_int (ii+npv+1                    ) = OBS_DUM
      obs% o% t_int (ii+npv+2                    ) = OBS_HS
      obs% o% lev   (ii+1        : ii+npv      :2) = log(p(1:nlev))
      obs% o% lev   (ii+2        : ii+npv      :2) = obs% o% lev(ii+1 : ii+npv :2)
      obs% o% lev   (ii+npv+1                    ) = 0._wp
      obs% o% lev   (ii+npv+2                    ) = log(100000._wp)
      !-------------------------------------
      ! request emissivity PC sink variables
      !-------------------------------------
      obs% o% t_int (ii+npv+nsv+1 : ii+ncv  ) = OBS_DUM
      obs% o% lev   (ii+npv+nsv+1 : ii+ncv  ) = 0._wp
      !-----------------------------------------------------
      ! request sink variable space for cloud top / fraction
      !                                    and snow fraction
      !-----------------------------------------------------
      obs% o% t_int (ii            +npv +3  ) = OBS_DUM
      obs% o% t_int (ii            +npv +4  ) = OBS_DUM
      obs% o% t_int (ii            +npv +5  ) = OBS_DUM
      if (spot% d% n == 0) then
        call new_sink (obs% o, spot, nbs)
        call set_sink (spot, obs% o, 1, cld_top_def, cld_top_e, npv+3, &
                       ctpmin, ctpmax, ctpmin_r, ctpmax_r, bound_constr)
        call set_sink (spot, obs% o, 2, cld_frc_def, cld_frc_e, npv+4, &
                       clfmin, clfmax, clfmin_r, clfmax_r, bound_constr)
        call set_sink (spot, obs% o, 3, 0._wp      , 0._wp    , npv+5, &
                       clfmin, clfmax, clfmin_r, clfmax_r, bound_constr)
      endif
      obs%   o% lev (ii       + npv+3) = log (        &
        obs% o% sink(spot%d%i +     1) % bg * 100._wp )
      obs% o% lev   (ii       + npv+4) = log (        &
        obs% o% sink(spot%d%i +     1) % bg * 100._wp )
      obs% o% lev   (ii       + npv+5) = log (100000._wp)

      !---------------------------------------------
      ! Set up t_tovs%av array (background profiles)
      !---------------------------------------------
      if (setup_tt) then
        l_cloud   = any(rs%iopts(1:rs%n_instr)%cloud_mode >= 1)
        l_o3      = any(rs%iopts(1:rs%n_instr)%use_o3     >  0)
        l_co2     = any(rs%iopts(1:rs%n_instr)%use_co2    >  0)
        l_t       = (fg_prof_top > 0) .and. (nld_t > 0) .or. (interp_strato > 0)
        l_q       = (fg_prof_top > 0) .and. (nld_h > 0) .or. (interp_strato > 0)
        if (ldeb(spot)) write(usd,*) dpref,'setup_full l',l_cloud,l_vint_rt,l_o3,l_co2,l_t,l_q
        if (ldeb(spot)) write(usd,*) dpref,'setup_full av',ttovs%nav,ttovs%av_cont(1:ttovs%nav)
        if (ldeb(spot)) write(usd,*) dpref,'setup_full levs',ttovs%nl_st,ttovs%nld_t,ttovs%nld_h
        ttovs%nav = 0
        if (l_t) call add_av_cont(ttovs, COL_T, i=ttovs%i_t)
        if (l_q) call add_av_cont(ttovs, COL_Q, i=ttovs%i_q)
        if (.not.l_vint_rt) then
          ! Since av is in single prec. it is better to use it for the p-profile
          ! only if absolutely necessary
          call add_av_cont(ttovs, COL_P, i=ttovs%i_p)
        end if
        if (l_cloud) then
          call add_av_cont(ttovs, COL_QCDIA)
          if (.not.mw_instr(rs%instr(1))) then
            call add_av_cont(ttovs, COL_CLC)
            call add_av_cont(ttovs, COL_QIDIA)
          end if
          if ( use_reff ) then
            call add_av_cont(ttovs, COL_REFF_QC)
            call add_av_cont(ttovs, COL_REFF_QI)
          end if
        end if
        if (l_o3 ) call add_av_cont(ttovs, COL_OZONE, i=ttovs%i_o3)
        if (l_co2) call add_av_cont(ttovs, COL_CO2  , i=ttovs%i_co2)
        allocate(ttovs%av(nlev, ttovs%nav)) ; ttovs%av = -1._tpp
        if (ttovs%i_p > 0) ttovs%av(1:nlev,ttovs%i_p) = p(1:nlev)
        ttovs%init = ior(ttovs%init, TTOVS_AV)
      end if

      if (setup_tt_tr) then
        ttovs%ntr = rs%bc%n_tr
        allocate(ttovs%tr(ttovs%nchan, ttovs%ntr))
        ttovs%tr = -1._sp
        ttovs%init = ior(ttovs%init, TTOVS_TR)
      end if

      if (setup_tt .or. setup_tt_tr) call store  (obs% o, spot, tovs=ttovs)

      if (ldeb(spot)) write(usd,*) dpref,'setup_full av',ttovs%nav,ttovs%av_cont(1:ttovs%nav)
      call destruct(ttovs)

      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) return
    endif

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
      !------------------------------------------
      ! setup observation error covariance matrix
      !   and variational quality control bounds
      !------------------------------------------
      if (obs% o% pe == dace% pe) then

        ld = ldeb(spot)
        allocate(ci(spot%o%n))
        call load  (obs% o, spot, ci=ci, rs=rs, ti=ti)

        ! Initialize obserror evaluation (diagonal)
        if (.not.associated(rs%oe_trf)) then
          allocate(rs%oe_trf(rs%n_chan))
          do i = 1, rs%n_chan
            call init(rs%oe_str(i), c_var, rs%oe_trf(i), stat, used_vnames=.true.)
            if (stat /= 0) then
              write(msg,'("rad_set ",I2," channel ",I5," stat=",I8)') rs%id,i,stat
              call finish(proc//'(TSK_R)', 'Failed to initialize obserror calculation for '//&
                   trim(msg)//' with "'//trim(rs%oe_str(i))//'"')
            end if
            m_var = max(m_var, rs%oe_trf(i)%nvar)
          end do
        end if
        i = -1 ; if (allocated(v)) i = size(v)
        if (m_var > i) then
          if (allocated(v)) deallocate(v)
          allocate(v(m_var+3))
        end if

        ! Calculate obserrors (diagonal)
        allocate(var(spot%o%n))
        instr = 1
        do i = 1, spot%o%n
          trf => rs%oe_trf(ci(i))
          if (i > ti%o_ch_i(instr) + ti%n_ch_i(instr)) then
            instr = instr + 1
            if (instr > size(ti%ii)) call finish(proc//'(TSK_R)', 'Invalid t_tovs_instr')
          end if
          call get_tovs_var(trf%vnames(1:trf%nvar), v(1:trf%nvar), stat, spot, &
               obs=obs% o, fg=y, chan=rs%chan(ci(i)), instr=rs%instr(ti%ii(instr)))
          if (stat /= 0) call finish(proc//'(TSK_R)', 'get_tovs_var failed: '//&
               trim(errmsg))
          var(i) = evaluate(v(1:trf%nvar),trf)
          if (ld) then
            write(usd,*) dpref,'oe_func ',rs%chan(ci(i)),trim(trf%func)
            do j = 1, trf%nvar
              write(usd,*) dpref,'oe_var ',rs%chan(ci(i)),j,trim(trf%vnames(j)),v(j)
            end do
            write(usd,*) dpref,'oe_result ',rs%chan(ci(i)),var(i),sqrt(var(i))
            write(usd,*) dpref,'oe_time ',rs%chan(ci(i)),spot%hd%time,ana_time
          end if              
        end do


        if (any(var(:) <= 0._sp)) then
          do j = 1, spot%o%n
            if (var(j) <= 0._sp) then
              if (.not. btest(rs%flag(ci(j)), USE_PASSIVE)) then
                write(0,*) '*** negative obserr',dace% pe,spot%hd%id,spot% stzen, &
                     spot% phase, spot% col% c% dlat
                write(0,*) j,ci(j),var(j),trim(rs%oe_trf(ci(i))%func)
              else
                ! Negative values might occur for instruments, that are mapped onto instruments with
                ! a wider swath. Since we do not assimilate mapped data, we should not worry about that.
              end if
              var(j) = 999._sp
              call decr_rpt_use(spot, CHK_OBS_ERR, use=STAT_PAS_REJ, comment='negative obserror')
            end if
          end do
        end if

        !--------------------------------
        ! load 'rules' for VQC parameters
        !--------------------------------
        call get_rule (type     = spot% hd% modtype,  &! <- module      type
                       obstype  = spot% hd% obstype,  &! <- observation type
                       codetype = spot% hd% codetype, &! <- code type
                       bf_type  = iud,                &! <- no BUFR     type
                       bf_subt  = iud,                &! <- no BUFR  subtype
                       db_kz    = iud,                &! <- no Datenbankkennzahl
                       stname   = '',                 &! <- no Station Name
                       lat      = spot% col% c% dlat, &! <- latitude
                       lon      = spot% col% c% dlon, &! <- longitude
                       o        = set                 )! -> channel information
        !-------------------------------------------
        ! set devault VQC parameters if not yet done
        !-------------------------------------------
        if (.not. associated (obs% o% s_vqc)) then
          allocate (obs% o% s_vqc (obs% o% n_obs))
          obs% o% s_vqc = svqc
        endif
        if (.not. associated (obs% o% f_vqc)) then
          allocate (obs% o% f_vqc (obs% o% n_obs))
          obs% o% f_vqc = vqc_form
        endif

        !---------------------------------
        ! store R in sparse representation
        !---------------------------------
        l = obs% R% ia (spot% o% i+1)
        !--------------------------------
        ! R is diagonal for the whole FOV
        !--------------------------------
        if (n_obserr_rad == 0) then
          do i=1,spot% o% n
             iii = spot% o% i + i
             obs% o% body  (iii)% eo = sqrt(var(i))
             obs% R% ia    (iii)     = l
             obs% R% packed  (l)     = var(i)
             obs% o% s_vqc (iii)     = set  % sgm_vq
             obs% o% f_vqc (iii)     = set  % frm_vq
             obs% R% ja      (l)     = iii
             l = l + 1
          end do
        !---------------------------------
        ! R may be nondiagonal
        ! loop over instruments in the FOV
        !---------------------------------
        else
l1:       do j = 1, ti%n_instr
            !--------------------------------------------------
            ! loop over pre-compiled R matrices for instruments
            !--------------------------------------------------
            do k = 1, n_obserr_rad
              !-----------------------------------------------
              ! cross-channel correlations for this instrument
              !-----------------------------------------------
              if (obs_corr_rad(k)% instrid == rs% instr_wmo(ti%ii(j))) then
                !------------------------------------------------------------
                ! instrument with nondiag.R found, set up:
                !   - address indeces k0, i0, i1
                !   - channel indices ix(:) to precomputed correlation matrix
                !   - sqrt(variances)
                !-----------------------------------------------------------
                k0 = k
                i0 = ti% o_ch_i(j) + 1
                i1 = ti% o_ch_i(j) + ti% n_ch_i(j)
                allocate ( ix(i0:i1) )
                allocate ( ee(i0:i1) )
                ix(:) = -1

                kk = 1
l2:             do i = i0, i1
                  ee(i) = sqrt (var(i))
                  iii = spot% o% i + i
l3:               do
                    if (kk > size (obs_corr_rad(k)% chan)) exit l3
                    if (obs% o% olev (iii) == obs_corr_rad(k)% chan(kk)) then
                      !---------------------------------------
                      ! corresponding channel found: set index
                      !---------------------------------------
                      if (obs_corr_rad(k)% use_chan(kk)) then
                         ix (i) = kk
                      else
                         ix (i) = -1
                      end if
                      kk = kk + 1
                      cycle l2
                    else if (obs% o% olev (iii) > obs_corr_rad(k)% chan(kk)) then
                      !----------------------------------------------
                      ! corresponding channel not yet found: try next
                      !----------------------------------------------
                      kk = kk + 1
                      cycle l3
                    else  ! obs% o% olev (iii) < obs_corr_rad(k)% chan(kk)
                      !---------------------------------
                      ! no corresponding channel present
                      !---------------------------------
                      exit l3
                    endif
                  end do l3
                  !-------------------------------------------------
                  ! no corresponding channel found:
                  !   abort for active observations,
                  !   treat as uncorrelated for passive observations
                  !-------------------------------------------------
                  if (obs% o% body(iii)% use% state < STAT_ACTIVE) then
                    ix(i) = -1
                  else
                    write (6,*) 'missing channel in correlated R: ',    &
                                 spot% statid,obs_corr_rad(k)% instrid, &
                                 obs% o% olev (iii),                    &
                                 obs% o% body (iii)% use% state
                    write (0,*) 'missing channel in correlated R: ',    &
                                 spot% statid,obs_corr_rad(k)% instrid, &
                                 obs% o% olev (iii),                    &
                                 obs% o% body (iii)% use% state
                    write (0,*) 'obs_corr_rad(k)% chan(:)',k, obs_corr_rad(k)% chan(:)
                    call finish (proc//'(TSK_R)','missing channel in R')
                  endif
                end do l2
              end if ! obs_corr_rad(k)% instrid == rs% instr_wmo(ti%ii(j))
            end do !  k = 1, n_obserr_rad

            if (allocated(ix)) then  ! instrument with full R
              !-------------------------------------------------
              ! correlation matrix found: obs_corr_rad(k0)% corr
              ! set up block-diagonal R in sparse representation
              !-------------------------------------------------
              do i = i0, i1
                iii = spot% o% i + i
                obs% o% body  (iii)% eo = sqrt (var(i))
                obs% R% ia    (iii)     = l
                obs% o% s_vqc (iii)     = set  % sgm_vq
                obs% o% f_vqc (iii)     = set  % frm_vq
                if (obs_corr_rad(k0)% frm_vq >= 0    ) obs% o% f_vqc (iii) = &
                    obs_corr_rad(k0)% frm_vq
                if (obs_corr_rad(k0)% sgm_vq >  0._wp) obs% o% s_vqc (iii) = &
                    obs_corr_rad(k0)% sgm_vq

                if (ldeb(spot)) write(usd,*) dpref,'corr. R',obs% o% olev (iii), ix(i)

                if (ix(i) == -1) then
                  obs% R% packed  (l)     = var(i)
                  obs% R% ja      (l)     = iii
                  l = l + 1
                else
                  do ii = i0, i1
                    jj = spot% o% i + ii
                    if (ix(ii) /= -1 ) then
                      if (ldeb(spot)) write(usd,*) dpref,'corr. R ii',&
                           ii, obs% o% olev (jj), obs_corr_rad(k0)% corr (ix(i),ix(ii))
                      if (obs_corr_rad(k0)% corr (ix(i),ix(ii)) /= 0._wp) then
                        obs% R% packed (l)  = obs_corr_rad(k0)% corr (ix(i),ix(ii)) &
                                           * ee(i) * ee(ii)
                        obs% R% ja     (l)  = jj
                        l = l + 1
                      endif
                    endif
                  end do
                endif
              end do

!--------------------------------------------------------------------------------------------------------
!    print diagonal elements of R matrix (debug)
!--------------------------------------------------------------------------------------------------------
              if (ldeb(spot)) then
                do ii = i0, i1
                  iii = spot% o% i + ii
                  if(ix(ii) /= -1 ) then
                    ll0=obs% R% ia(iii)
                    ll1=obs% R% ia(iii+1)-1
                    do ll = ll0,ll1
                      if(obs% R% ja(ll) == iii) then
                        write(usd,*) dpref,'diag_R',ii,obs% R% packed (ll) , sqrt(var(ii)),var(ii)
                      endif
                    enddo
                  endif
                enddo
              endif
!--------------------------------------------------------------------------------------------------------

              deallocate (ix, ee)

            else  ! ix not allocated
              !-------------------------------
              ! no correlation matrix found
              ! use diagonal R for this instrument
              !-------------------------------

              do i = ti% o_ch_i(j) + 1, ti% o_ch_i(j) + ti% n_ch_i(j)
                iii = spot% o% i + i
                obs% o% body  (iii)% eo = sqrt (var(i))
                obs% R% ia    (iii)     = l
                obs% R% packed  (l)     = var(i)
                obs% o% s_vqc (iii)     = set  % sgm_vq
                obs% o% f_vqc (iii)     = set  % frm_vq
                obs% R% ja      (l)     = iii
                l = l + 1
              end do

            endif

          end do l1
        endif ! n_obserr_rad == 0
        obs% R% ia (spot% o% i + spot% o% n + 1) = l
        ! call destruct(ttovs)
        deallocate(ci, var)
      endif ! obs% o% pe == dace% pe
      tsk = tsk - TSK_R
      if (tsk == 0) return
    endif

    !=========================
    ! run observation operator
    !=========================
    if (iand (TSK_Y,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
FTRACE_BEGIN('process_tovs:calc_y')
        allocate(tb(spot% o% n))
        call prep_rttov_prof(obs, spot, xi, rtp, ifail, ttovs=ttovs)
        if (ifail /= 0) call finish (proc//'(TSK_Y)','prep_rttov_prof failed')
        call call_rttov ('f',ttovs, rtp, ifail, tb=tb, ldebug=ldeb(spot))
        if(ifail/=0) call finish (proc//'(TSK_Y)','call_rttov failed')
        if (present(y)) y% x (spot% o% i+1:spot% o% i+spot% o% n) = tb(:)
        call destruct(ttovs)
        call destruct(rtp)
FTRACE_END('process_tovs:calc_y')
      endif
      tsk = tsk - TSK_Y
      if (tsk == 0) return
    endif

    !=========================
    ! set up H (Jacobi-matrix)
    !=========================
    if (iand (TSK_K,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
FTRACE_BEGIN('process_tovs:calc_k')
!        call run_rttov ('k', obs, spot, xi, y)
        allocate(tb(spot% o% n))
        call prep_rttov_prof(obs, spot, xi, rtp, ifail, ttovs=ttovs, tjac=tjac)
        if (ifail /= 0) call finish (proc//'(TSK_K)','prep_rttov_prof failed')
        call construct (rtp_k, ttovs% nchan, jakobi=1)
        call call_rttov ('k',ttovs, rtp, ifail, tb=tb, x_a=rtp_k, ldebug=ldeb(spot))
        if(ifail/=0) call finish (proc//'(TSK_K)','call_rttov failed')

        allocate (Hnew(spot% o% n, spot% i% n))
        rtp_k% sav (3,:) = rtp_k% sav (3,:) / 100._wp         ! hPa -> Pa
        if (no_ps_dep) then
          rtp_k% sav (3,:) = 0._wp   ! temporary workarount for ..
          tjac% ps_llev    = .true.  ! .. discontinuity in RTTOV
        endif
        call K2H(H   = Hnew,             &
             temp_k  = rtp_k% av(:,1,:), &
             humi_k  = rtp_k% av(:,2,:), &
             t2m_k   = rtp_k% sav( 1,:), &
             q2m_k   = rtp_k% sav( 2,:), &
             psurf_k = rtp_k% sav( 3,:), &
             stemp_k = rtp_k% ssv( 1,:), &
             ctp_k   = rtp_k% cv ( 1,:), &
             cfr_k   = rtp_k% cv ( 2,:), &
             emis_k  = rtp_k% emis(1,:), &
             tj      = tjac,             &
             pz_bg   = spot% pz_bg,      &
             n_chan  = ttovs% nchan,     &
             spot    = spot)
        call insert (obs% H, Hnew, spot% o% i, spot% i% i, use_zero=.true.)

        obs% yi% x (spot% o% i+1:spot% o% i+spot% o% n) = tb(:)
        if (present(y)) y% x (spot% o% i+1:spot% o% i+spot% o% n) = tb(:)

        call store (obs% o, spot, ttovs)
        call destruct(ttovs)
        call destruct(tjac)
        call destruct(rtp_k)
        call destruct(rtp)
FTRACE_END('process_tovs:calc_k')
      endif
      tsk = tsk - TSK_K
      if (tsk == 0) return
    endif

    !==========================
    ! abort if any task is left
    !==========================
    write(0,*) 'process_tovs: unknown task =',tsk
    call finish(proc,'unknown task')

  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_nr(ti)
      type(t_tovs_instr), optional, target :: ti

      type(t_tovs_instr), target  :: ti_dum
      type(t_tovs_instr), pointer :: ti_p
      integer                     :: i
      integer                     :: ci(m_chan)
    !-------------------------------------------------------
    ! set upper bound for R-matrix size for the current spot
    !-------------------------------------------------------
      if (n_obserr_rad == 0) then
        spot% nr        = spot% o%n
      else
        if (present(ti)) then
          ti_p => ti
        else
          call load  (obs% o, spot, rs=rs, ti=ti_dum, ci=ci)
          ti_p => ti_dum
        end if
        spot% nr = 0
        do i = 1, ti_p% n_instr
          if (any (obs_corr_rad(1:n_obserr_rad)% instrid == rs% instr_wmo(ti_p%ii(i)))) then
            spot% nr = spot% nr + ti_p% n_ch_i(i) ** 2
          else
            spot% nr = spot% nr + ti_p% n_ch_i(i)
          endif
        end do
      endif
    end subroutine set_nr

  end subroutine process_tovs


  subroutine print_tovs_levs(obs)
    type(t_obs_set), intent(inout) :: obs
    ! print information on vertical levels
    character(len=15), parameter   :: proc = 'print_tovs_levs'
    type(t_spot),      pointer     :: spt
    type(t_rad_set),   pointer     :: rs
    type(t_tovs)                   :: ttovs
    character(len=6)               :: typ
    integer,           allocatable :: n       (  :)
    real(wp),          allocatable :: nlev    (:,:)
    real(wp),          allocatable :: nl_st(:,:)
    real(wp),          allocatable :: nld_t   (:,:)
    real(wp),          allocatable :: nld_q   (:,:)
    real(wp),          allocatable :: p_ext(:,:), p(:)
    real(wp),          allocatable :: p_top(:,:), p_mod(:,:), p_d_t(:,:), p_d_q(:,:)
    real(tpp),         allocatable :: av(:,:)
    integer                        :: i, ib, is, i_rs, nn, np, nl

    if (n_set <= 0) return
    do i = 1, n_set
      rad_set(i)% n_lev = p_max (rad_set(i)% n_lev)
    end do
    allocate(n(n_set),nlev(n_set,3),nl_st(n_set,3),nld_t(n_set,3),nld_q(n_set,3),&
             p_top(n_set,3), p_mod(n_set,3), p_d_t(n_set,3), p_d_q(n_set,3))
    n=0
    call init_(nlev)  ; call init_(nl_st) ; call init_(nld_t) ; call init_(nld_q)
    call init_(p_top) ; call init_(p_mod) ; call init_(p_d_t) ; call init_(p_d_q)
    mx_nlev   = p_max(mx_nlev)
    mx_nav    = p_max(mx_nav )
    allocate(av(mx_nlev, mx_nav), p(mx_nlev))
    np = 0
    do ib = 1, size(obs%o)
      if (dace%pe /= obs%o(ib)%pe) cycle
      do is = 1, obs%o(ib)%n_spot
        spt => obs%o(ib)%spot(is)
        if (spt% hd% obstype /= OT_RAD) cycle
        call load(obs%o(ib), spt, tovs=ttovs, tovs_io=TTOVS_BASE, i_rs=i_rs, av=av)
        if (i_rs < 1 .or. i_rs > n_set) call finish(proc, 'invalid i_rs')
        rs => rad_set(i_rs)
        nl = ttovs%nlev
        if (ttovs%i_p > 0) then
          p(1:nl) = real(av(1:nl,ttovs%i_p),wp)
        else
          p(1:nl) = preslev(1:nl)
        end if
        call add(real(ttovs%nlev,wp) , nlev (i_rs,:), n(i_rs))
        call add(real(ttovs%nl_st,wp), nl_st(i_rs,:), n(i_rs))
        call add(real(ttovs%nld_t,wp), nld_t(i_rs,:), n(i_rs))
        call add(real(ttovs%nld_h,wp), nld_q(i_rs,:), n(i_rs))
        call add(p(1)                , p_top(i_rs,:), n(i_rs))
        call add(p(ttovs%nl_st+1)    , p_mod(i_rs,:), n(i_rs))
        call add(p(ttovs%nld_t+1)    , p_d_t(i_rs,:), n(i_rs))
        call add(p(ttovs%nld_h+1)    , p_d_q(i_rs,:), n(i_rs))
        n(i_rs) = n(i_rs) + 1
        if (rs%gopts%lev_mode > 0 .and. interp_strato > 0) then
          if (.not.allocated(p_ext)) then
            allocate(p_ext(nlev_ext,3))
            call init_(p_ext)
          end if
          do i = 1, nlev_ext
            call add(p(i), p_ext(i,:), np)
          end do
          np = np + 1
        end if
        call destruct(ttovs)
      end do
    end do
    nn = sum(n(:))
    nn = p_sum(nn)
    if (nn == 0) return

    if (dace%lpio) then
      write(*,*)
      write(*,*) 'Radiance levels info:'
      write(*,*) ' sat|grid|   typ|       #levels| #extra levels| #T.dum.levels| #Q.dum.levels|'
      write(*,*) '    |    |      |  mean|min|max|  mean|min|max|  mean|min|max|  mean|min|max|'
      write(*,*) ' ---|----|------|------|---|---|------|---|---|------|---|---|------|---|---|'
    end if
    do i_rs = 1, n_set
      rs => rad_set(i_rs)
      n(i_rs) = p_sum(n(i_rs))
      if (rs%id > 0 .and. n(i_rs) > 0) then
        call p_avg(nlev (i_rs,:), n(i_rs))
        call p_avg(nl_st(i_rs,:), n(i_rs))
        call p_avg(nld_t(i_rs,:), n(i_rs))
        call p_avg(nld_q(i_rs,:), n(i_rs))
        if (dace%lpio) then
          if (rs%gopts%lev_mode > 0) then
            typ = 'model'
            if (interp_strato > 0) typ = 'model+'
          else
            typ = 'rttov'
          end if
          write(*,'(1x,2(1x,I3.3,"|"),A6,"|",4(1x,F5.1,"|",I3,"|",I3,"|"))') rs%satid,rs%grid,typ,&
               nlev (i_rs,1),int(nlev(i_rs,2:3)) , nl_st(i_rs,1),int(nl_st(i_rs,2:3)),&
               nld_t(i_rs,1),int(nld_t(i_rs,2:3)), nld_q(i_rs,1),int(nld_q(i_rs,2:3))
        end if
      end if
    end do

    if (dace%lpio) then
      write(*,*)
      write(*,*) ' sat|grid|                     top level|               top model level|&
           &     highest non-dummy-T level|     highest non-dummy-Q level|'
      write(*,*) '    |    |      mean|      min|      max|      mean|      min|      max|&
           &      mean|      min|      max|      mean|      min|      max|'
      write(*,*) ' ---|----|----------|---------|---------|----------|---------|---------|&
           &----------|---------|---------|----------|---------|---------|'
    end if
    do i_rs = 1, n_set
      rs => rad_set(i_rs)
      if (rs%id > 0 .and. n(i_rs) > 0) then
        call p_avg(p_top(i_rs,:), n(i_rs))
        call p_avg(p_mod(i_rs,:), n(i_rs))
        call p_avg(p_d_t(i_rs,:), n(i_rs))
        call p_avg(p_d_q(i_rs,:), n(i_rs))
        if (dace%lpio) &
             write(*,'(1x,2(1x,I3.3,"|"),4(1x,3(F9.4,"|")))') &
             rs%satid,rs%grid,p_top(i_rs,:)*0.01_wp,p_mod(i_rs,:)*0.01_wp,&
             p_d_t(i_rs,:)*0.01_wp,p_d_q(i_rs,:)*0.01_wp
      end if
    end do


    np = p_sum(np)
    if (np > 0) then
      if (dace%lpio) then
        write(*,*)
        write(*,*) 'Radiance additional top levels (filled by interp_strato):'
        write(*,*) '  |     mean p|      min p|      max p|'
      end if
      if (.not.allocated(p_ext)) then
        allocate(p_ext(nlev_ext,3))
        call init_(p_ext)
      end if
      do i = 1, nlev_ext
        call p_avg(p_ext(i,:), np)
        if (dace%lpio) write(*,'(1x,I2,"|",3(F11.6,"|"))') i,p_ext(i,:)*0.01_wp
      end do
    end if

  contains

    subroutine add(x, xs, n)
      real(wp), intent(in)    :: x
      real(wp), intent(inout) :: xs(:)
      integer,  intent(in)    :: n
      if (n == 0) then
        xs(:) = x
      else
        xs(1) = xs(1) + x
        xs(2) = min(xs(2), x)
        xs(3) = max(xs(3), x)
      end if
    end subroutine add

    subroutine p_avg(xs, n)
      real(wp), intent(inout) :: xs(3)
      integer,  intent(in)    :: n
      xs(1) = p_sum(xs(1))/real(n,wp)
      xs(2) = p_min(xs(2))
      xs(3) = p_max(xs(3))
    end subroutine p_avg

    subroutine init_(xs)
      real(wp), intent(inout) :: xs(:,:)
      integer :: i
      do i = 1, size(xs,1)
        xs(i,1) = 0._wp
        xs(i,2) = huge(0._wp)
        xs(i,3) = 0._wp
      end do
    end subroutine init_

  end subroutine print_tovs_levs



  subroutine destruct_tovs_oe
    ! Destruct the t_range_fparse structures in rad_set
    integer :: i
    do i = 1, n_set
      if (associated(rad_set(i)%oe_trf)) then
        call destruct(rad_set(i)%oe_trf(:))
        deallocate(rad_set(i)%oe_trf)
      end if
    end do
  end subroutine destruct_tovs_oe

!---------------------------------------------------------------------

  subroutine process_tovs_mult(tsk, obs, atm, xi, y, mask_rs)
  !----------------------------------------
  ! alternative version for TSK_Y, TSK_K:
  ! call rttov_direct / rttov_k in parallel
  !----------------------------------------
  integer            ,intent(in)             :: tsk    ! what to do
  type(t_obs_set)    ,intent(inout)          :: obs    ! observation data type
  type(t_atm)        ,intent(in)   ,optional :: atm    ! atmospheric state
  type(t_vector)     ,intent(in)   ,optional :: xi     ! interpolated values
  type(t_vector)     ,intent(inout),optional :: y      ! observed quantity
  logical            ,intent(in)   ,optional :: mask_rs(:) ! mask rad_set, i.e. exclude some datasets

    character(len=17),   parameter           :: proc = 'process_tovs_mult'
    character(len=256)                       :: msg  = ''
    type(t_rad_set),     pointer             :: set(:)
    type(t_radv),        allocatable, target :: rad(:)
    type(t_radv),        pointer             :: r   => NULL()
    type(t_rad_set),     pointer             :: s   => NULL()
    type(t_spot),        pointer             :: spt => NULL()
    type(t_obs_block)                        :: ob
    type(t_tovs)                             :: ttovs
    ! type(t_tovs_instr)                       :: ti
    type(t_jac_arr),     allocatable         :: tjs(:)
    integer,             allocatable         :: lprofs(:)
    integer,             allocatable         :: lchans(:)
    integer,             allocatable         :: istore(:,:)
    integer,             allocatable         :: ivld(:)
    integer                                  :: ib, is, iset, instr, i, j, k, l, m, ib_prev
    integer                                  :: n_prof
    integer                                  :: n_chan
    integer                                  :: n_calc
    integer                                  :: m_calc  ! max(n_calc) for allocates
    integer                                  :: n_lev, n_lay, n_lev_, nlt
    integer                                  :: n_lev_rt
    integer                                  :: ncv ! sizeof interpolationspace
    integer                                  :: ierr
    integer                                  :: ioff
    integer                                  :: nch
    integer                                  :: ipr(5), ipr_(5), npr, npr_
    logical,             allocatable         :: mask(:)
    logical                                  :: l_l2c, l_l2c_i
    logical                                  :: l_im, l_im_all  ! imager first guess calcs.
    logical                                  :: l_btcs
    logical                                  :: l_debug
    logical                                  :: l_addinterp
    logical                                  :: l_cloud
    logical                                  :: l_aux
    logical                                  :: l_dealloc
    logical                                  :: l_new
    ! Loop over RTTOV calls
    integer                                  :: n_calc_aux
    integer                                  :: n_max_calc
    integer                                  :: i_start
    integer                                  :: i_end
    integer                                  :: n_prof_aux
    integer                                  :: n_prof_block
    integer                                  :: n_max_prof
    integer                                  :: i_prof_start
    integer                                  :: i_prof_end
    integer                                  :: i_prof
    integer                                  :: ip0, ip1
    integer                                  :: npc     ! prof*chan/block
    integer                                  :: rad_out
    ! surface influence
    logical                                  :: l_surf_infl
    logical,             save                :: l_first = .true.
    integer                                  :: mode
    real(sp),            allocatable         :: sinfl(:,:)
    real(sp),            allocatable, target :: sinf(:)
    ! Information on context of RTTOV call and test information
    integer,             save                :: po_context_old = -1
    integer,             save                :: po_ilns_old    = -1
    logical,             save                :: po_lmon_old    = .false.
    ! Checks on operator applicabilty
    integer,             pointer             :: act_opdep   => null()
    integer,             pointer             :: act_ps      => null()
    integer,             pointer             :: act_rl      => null()
    integer,             pointer             :: act_god     => null()
#if (_RTTOV_VERSION > 0)
    logical                                  :: l_chk_opdep =  .false.
    logical                                  :: l_chk_ps    =  .false.
#endif
    logical                                  :: l_chk_rl    =  .false.
    logical                                  :: l_chk_god   =  .false.
    type(t_keep_rtstat), pointer             :: krt         => null()
    logical                                  :: l_set
    integer,             save                :: i_opdep_tsk = -1
!   integer,             save                :: i_god_tsk   = -1
    integer,             save                :: i_ps_tsk    = -1
    integer,             save                :: nout_chk_reg_lims = 0
#if (_RTTOV_VERSION >= 12)
    integer,             save                :: nout_chk_god      = 0
#endif
    integer                                  :: nrl_lev     = 0
    integer                                  :: nrl_var     = 0
    integer                                  :: nrl_dir     = 0
    integer                                  :: ni, n1
    logical                                  :: l_rl_var    = .false.
    logical                                  :: l_rl_lev    = .false.
    logical                                  :: l_rl_dir    = .false.
    logical                                  :: l_new_out   = .false.
    logical                                  :: l_mw        = .false.
    logical,             allocatable         :: l_vis(:)
    logical,             allocatable         :: l_rl(:,:)
    integer,             allocatable         :: reg_lim(:,:,:)
    integer,             allocatable         :: rflag(:,:)
    real(wp),            allocatable         :: radrefl(:,:)
    real(wp),            allocatable         :: refl(:,:)
    integer,             allocatable         :: quality(:,:)
    ! K matrix
    real(wp),            allocatable         :: emiss_k(:,:)
    real(wp),            allocatable         :: temp_k (:,:,:)
    real(wp),            allocatable         :: humi_k (:,:,:)
    real(wp),            allocatable         :: t2m_k  (:,:)
    real(wp),            allocatable         :: q2m_k  (:,:)
    real(wp),            allocatable         :: psurf_k(:,:)
    real(wp),            allocatable         :: ctp_k  (:,:)
    real(wp),            allocatable         :: cfr_k  (:,:)
    real(wp),            allocatable         :: u10m_k (:,:)
    real(wp),            allocatable         :: v10m_k (:,:)
    real(wp),            allocatable         :: stemp_k(:,:)
    real(wp),            allocatable         :: r_clear(:,:)
    real(wp),            allocatable         :: overc  (:,:,:)
    real(wp),            allocatable         :: transm (:,:,:)
    real(wp),            allocatable         :: transmcld (:,:,:)
    real(wp),            allocatable         :: opdep  (:,:,:)
    real(wp),            allocatable         :: dppmvdq(:)
    real(wp)                                 :: dppmvdq_2m
    real(wp),            allocatable         :: Hnew(:,:)
    real(sp),            allocatable         :: l2c_aux(:,:)
    real(sp),            allocatable         :: l2c_corr(:,:)
    real(sp),            allocatable         :: l2c    (:)
    real(wp),            allocatable         :: emis   (:)
    integer                                  :: sign_opdep
    real(sp),            allocatable         :: bt_cs  (:)
    real(wp),            allocatable         :: bt_fg_cs(:,:)
    ! write optical depths for derivation of god parameters
    logical                                  :: wr_opd
#if (_RTTOV_VERSION >= 12)
    type(t_wr_opd),      allocatable         :: wo_snd(:)
    integer                                  :: n_wo, i_wo
#endif
    integer                                  :: lb1_rf, ub1_rf
    ! write profiles for spt_hd_debug
    integer                                  :: n_wr_profs
    character(len=300)                       :: wr_profs_fmt, sdum, descr_call
    ! Emissivity atlases
    integer,             parameter           :: m_opt = m_instr * m_rad_set * 3
    type(t_emis_opt),    pointer             :: eo
    integer                                  :: n_opt = 0
    integer                                  :: i_opt   (m_emis)
    integer                                  :: atlas_id(m_emis)
    logical                                  :: angcorr (m_emis)
    ! BRDF atlas(es)
    integer                                  :: n_brdf
    integer                                  :: i_opt_brdf(m_rad_set)
    ! emissivities
    integer                                  :: ityp, ir_atlas
    integer                                  :: mw_emis_mod
    ! height/width calculation (for LETKF)
    logical                                  :: l_hgt
    logical                                  :: l_vis_
    ! specularity
    real(wp),            pointer             :: specularity(:,:)
    real(wp),            target              :: spec_dum(0,0)
    ! transmission weighted BC predictors
    logical                                  :: l_pr_tr
    ! imager first guess
    integer                                  :: i_imch_v
    real(wp),            allocatable         :: im_fg(:,:), im_emis(:,:)
    real(sp),            allocatable         :: im_ch(:,:)
    integer                                  :: n_im_fg, n_im_emis, n_im_emis_f
    integer                                  :: n_im_tskin, n_im_tskin_f, n_im_tskin_styp(0:n_styp)

    !======================================================
    ! TSK_SETUP_OP: setup RTTOV (and rad_set, if necessary)
    !======================================================
    if (iand (TSK_SETUP_OP,tsk) /= 0) then

      ! Old DEPRECATED way to set RTTOV options. Kept for backwards compatibilty
      call rtifc_set_opts(tmpl              = 'default0',       &!
                          apply_reg_lims    = app_reg_lims,     &!
                          fastem_version    = rt_mw_emis_mod,   &!
                          ir_sea_emis_model = rt_ir_emis_mod,   &!
                          do_lambertian     = rt_do_lambertian, &!
                          addinterp         = .false.,          &!
                          use_q2m           = rt_use_q2m,       &!
                          fix_hgpl          = fix_hgpl,         &!
                          clip_gas_opdep    = rt13_clip_gas_opdep)

      ! Set instrument specific RTTOV options
      ! TODO: idea: use rs%iflag to mark all used channels and load only the required
      !             rtcoef files.
      do iset = 1, n_set
        s => rad_set(iset)
        do instr = 1, s%n_instr
          l_new = .true.
          if (s%rttov_indx(instr) > 0 .and. associated(s%chidx)) then
            l_new = any(s%chidx(s%o_ch_i(instr)+1:s%o_ch_i(instr)+s%n_ch_i(instr)) <= 0)
          end if
          if (l_new) then
            write(msg, '("rttov options for satid ",I3.3,", instr ",I3.2,", grid ",I3.2)') &
                 s%satid, s%instr(instr), s%grid
            l_mw = mw_instr(s%instr(instr))
            if (s%iopts(instr)%rt_mw_emis_mod > 0) then
              mw_emis_mod = s%iopts(instr)%rt_mw_emis_mod
            else
              mw_emis_mod = rt_mw_emis_mod
            end if
            call rtifc_set_opts(iopt       = s%rttov_indx(instr),            &
                                new        = .true.,                         &
                                name       = trim(msg),                      &
                                satid      = s%satid,                        &
                                grid       = s%grid,                         &
                                instr      = s%instr(instr),                 &
                                addclouds  = s%iopts(instr)%cloud_mode > 0 .and. .not.l_mw,  &
                                addinterp  = s%gopts%lev_mode >= 1,          &
                                clw_data   = s%iopts(instr)%cloud_mode > 0 .and. l_mw,  &
                                interp_mode= vint_rttov,                     &
                                ozone_data = s%iopts(instr)%use_o3 > 0,      &
                                co2_data   = s%iopts(instr)%use_co2 > 0,     &
                                conv_overc = .not.s%iopts(instr)%l2c_use_rad,&
                                do_lambertian = s%iopts(instr)%do_lambertian,&
                                fastem_version= mw_emis_mod                  &
                                )
            if (associated(s%iflag)) then
              l_vis_ = any(s%iflag(s%o_ch_i(instr)+1:s%o_ch_i(instr)+s%n_ch_i(instr)) == 1)
            else
              l_vis_ = .false.
            end if
            if (l_vis_) then
              if (rtifc_vers < 13) call finish(proc, 'RTTOV_VERSION >= 13 for visible observations required')
              call rtifc_set_opts(iopt=s%rttov_indx(instr), &
                                  addsolar        = .true., &
                                  vis_scatt_model = 3       )
            end if
          end if
        end do
        if (s%im_rttov_id >= 0) then
          write(msg, '("rttov options for satid ",I3.3,", instr ",I3.2,", grid ",I3.2)') &
               s%satid, s%im_rttov_id, s%grid
          call rtifc_set_opts(iopt       = s%im_rttov_indx,                &
                              new        = .true.,                         &
                              name       = trim(msg),                      &
                              satid      = s%satid,                        &
                              grid       = s%grid,                         &
                              instr      = s%im_rttov_id,                  &
                              addclouds  = .false.,                        &
                              addinterp  = s%gopts%lev_mode >= 1,          &
                              interp_mode= vint_rttov                      &
                             )
        end if
      end do

#if !defined(__ICON__)
      min_od                  = rt_min_od  ! minimum optical density in RTTOV
#endif
      ! Set RTIFC options
      default_gas_units       = rt_humi
      out_path                = aux
      default_salinity        = rt_salinity

      if (n_set > 0) then
        if (rttov_version > 0) then
          call call_rttvi (data_path)
        end if

        !---------------------------------
        ! Initialise MW atlases if needed
        !---------------------------------
#if (_RTTOV_VERSION >= 12)
        n_opt  = 0
        n_brdf = 0
        do iset = 1, n_set
          s => rad_set(iset)
          ir_atlas = -1
          if (any(s%emis_opt(1:s%n_emis_opt)%mode == MODE_ATLAS)) then
            do j = 1, s% n_instr
              ityp = instr_type(s, j)
              ioff = s%o_ch_i(j)
              nch  = s%n_ch_i(j)
              do k = 1,s%n_emis_opt
                eo => s%emis_opt(k)
                if (eo%inst == s%instr(j) .and. eo%mode == MODE_ATLAS .and. &
                    any(s%emis_ind(:,:,ioff+1:ioff+nch) == k)) then
                  if (btest(ityp, ITYP_VIS) .and. eo%atls == ATLS_BRDF) then
                    n_brdf = n_brdf + 1
                    if (n_brdf > size(i_opt_brdf)) call finish(proc, 'i_opt_brdf too small')
                    i_opt_brdf(n_brdf) = s%rttov_indx(j)
                  else
                    n_opt = n_opt + 1
                    if (n_opt > m_opt) call finish(proc, 'm_opt too small')
                    i_opt   (n_opt) = s%rttov_indx(j)
                    atlas_id(n_opt) = eo%atls
                    angcorr (n_opt) = eo%angcorr
                  end if
                  if (btest(ityp, ITYP_IR)) ir_atlas = eo%atls
                end if
              end do
            end do
          end if
          if (s% im_rttov_id >= 0 .and. s%n_im_chan > 0) then
            if (s% gopts% im_atlas_id <= 0) s% gopts% im_atlas_id = ir_atlas
            if (s% gopts% im_atlas_id >  0) then
              ! Atlas for additional imager data
              n_opt           = n_opt + 1
              i_opt   (n_opt) = s%im_rttov_indx
              atlas_id(n_opt) = ir_atlas
              angcorr (n_opt) = .true.
            end if
          end if
        end do
        if (n_opt > 0) then
          call rtifc_init_atlas(i_opt(1:n_opt), atlas_id(1:n_opt), angcorr(1:n_opt),&
               imm(atm%time), data_path, dace%pe, dace%npe, dace%pio, dace%comm)
        end if
        if (n_brdf > 0) then
          call rtifc_init_brdf_atlas(i_opt_brdf(1:n_brdf),imm(atm% time), data_path, &
               dace% pe, dace% npe, dace% pio, dace% comm, ierr)
          if (ierr /= 0 )   call finish ('process_tovs','init_rttov_brdf_atls failed')
        end if
#endif
      end if
    end if

    !----------------------------------------------------
    ! Only TSK_SETUP_OP, TSK_Y and TSK_K implemented here
    !----------------------------------------------------
    if (iand(tsk, TSK_Y+TSK_K) == 0) then
      return
    end if

    !--------------------------------------------
    ! return if no valid satellite data available
    !--------------------------------------------
    if (n_set <= 0) return

    call prep_checks
    wr_opd_disc = btest(wr_opdep, 4) .and. btest(wr_opdep, 6)

    ! Prepare writing of profiles within RTTOV
    n_wr_profs = 0
    if (spt_debug_prof /= 0) then
      n_wr_profs = count(spt_hd_debug > 0)
      if (n_wr_profs > 0) then
        if (tsk == tsk_k) then
          wr_profs_fmt='_K'
        else
          wr_profs_fmt='_Y'
        end if
        select case(po_context)
        case(POC_FG)
          wr_profs_fmt = trim(wr_profs_fmt)//'_FG'
        case(POC_ANA)
          wr_profs_fmt = trim(wr_profs_fmt)//'_ANA'
        case(POC_TOP)
          wr_profs_fmt = trim(wr_profs_fmt)//'_TOP'
        case default
          write(sdum,'("_it",I2.2)') po_context
          wr_profs_fmt = trim(wr_profs_fmt)//trim(sdum)
        end select
        if (po_ilns /= -99) then
          write(sdum,'("_ilns",I3.3)') po_ilns
          wr_profs_fmt = trim(wr_profs_fmt)//trim(sdum)
        end if
        if (po_lmon) &
             wr_profs_fmt = trim(wr_profs_fmt)//'_mon'
        wr_profs_fmt = 'profiles",I12.12,"'//trim(wr_profs_fmt)//'.H5")'
        wr_profs_fmt = path_file(aux, trim(wr_profs_fmt))
        wr_profs_fmt = '("'//trim(wr_profs_fmt)
        if (dace% lpio) write(6,'("    write debug profiles into ",A)') trim(wr_profs_fmt)
      end if
    end if

    if (trace_mem_time) then
      descr_call = ''
      select case(po_context)
      case(POC_FG)
        descr_call = 'fg'
      case(POC_ANA)
        descr_call = 'ana'
      case(POC_TOP)
        descr_call = 'top'
      case default
        write(descr_call,'("iter",I2.2)') po_context
      end select
    end if

    l_hgt = (po_context == POC_FG) .and. pl_method >= 5
    l_im_all = .false.

FTRACE_BEGIN("tovs_mult_prof")
    !--------------------------------------------
    ! fill radiance data structure of type t_radv
    !--------------------------------------------
    if (trace_mem_time) call stop_time ('process_tovs_mult '//trim(descr_call)//': call fill_rad')

    set => rad_set(1:n_set)
FTRACE_BEGIN("tovs_mult_prof:setup")
    ! get the datasets
    allocate(rad(n_set))
    if (tsk == TSK_K) then
      allocate(tjs(n_set))
      call fill_rad (rad, set, obs, xi, y=y, lrttov=.true., lprint=.false.,&
                                                                tj_set=tjs, mask_rs=mask_rs)
    else
      call fill_rad (rad, set, obs, xi, y=y, lrttov=.true., lprint=.false., mask_rs=mask_rs)
    end if
FTRACE_END("tovs_mult_prof:setup")

    !------------------------------------------------------------
    ! loop over sets of instruments interpolated to the same grid
    !------------------------------------------------------------
!NEC$ nomove
    loop_datasets: do iset = 1, n_set
      if (present(mask_rs)) then
        if (.not.mask_rs(iset)) CYCLE
      end if
      if (trace_mem_time) call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': loop_datasets ')
      r => rad(iset)
      s => set(iset)
      !--------------------------
      ! determine the array sizes
      !--------------------------
      n_prof = r%n_rec
      n_lev  = r%n_lev

      wr_opd = (any(s%instr(1) == (/ 10, 15, 16, 19, 21, 22, 34, 44, 56, 71, 73 /))) .and. &
               (btest(wr_opdep, 0) .and. po_context == POC_FG                .or.       &
                btest(wr_opdep, 1) .and. po_context == POC_ANA               .or.       &
                btest(wr_opdep, 2) .and. po_context >=  1 .and. tsk == tsk_k .or.       &
                btest(wr_opdep, 3) .and. po_context >=  1 .and. tsk == tsk_y .or.       &
                btest(wr_opdep, 4) .and. (po_context <= POC_TOP .or. po_lmon))
      if (wr_opd .and. wr_opd_disc .and. po_lmon .and. po_ilns < 0) wr_opd = .false. ! Prevent output from test_gradient

      n_prof_gt_0: if (n_prof <= 0) then
        !-----------------------------------------
        ! cycle loop in case of no data on this PE
        !-----------------------------------------
        if (trace_mem_time) then
          !------------------------------------------
          ! dummy calls if respective code is skipped
          !------------------------------------------
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': allocate arrays for rttov call')
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': allocate arrays for rttov_k')
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': Loop over blocks of profiles')
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': deallocate')
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': deallocate k')
        endif
        call write_opdep('alloc')
      else

        if (tsk == TSK_K) then
          n_max_calc = s% gopts% n_max_calc_k ; if (n_max_calc <= 0) n_max_calc = n_max_calc_k_def
          n_max_prof = s% gopts% n_max_prof_k ; if (n_max_prof <= 0) n_max_prof = n_max_prof_k_def
        else
          n_max_calc = s% gopts% n_max_calc_y ; if (n_max_calc <= 0) n_max_calc = n_max_calc_y_def
          n_max_prof = s% gopts% n_max_prof_y ; if (n_max_prof <= 0) n_max_prof = n_max_prof_y_def
        end if
        n_prof_block = min (n_max_prof, n_prof)
        m_calc       = max(maxval (s%n_ch_i(1:s%n_instr)),s%n_im_chan) * n_prof_block
        n_chan       = s%n_chan

        ! Allocate keep_rtstat%spt if required
        if (associated(keep_rtstat)) then
          krt => keep_rtstat(iset)
          if (.not.associated(krt%spt)) then
            allocate(krt%spt(n_prof))
            krt%nspot = n_prof
          end if
        end if

        !-----------------------------------
        ! allocate arrays for the rttov call
        !-----------------------------------
        if (trace_mem_time)                                                  &
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': allocate arrays for rttov call')

        allocate(lprofs (m_calc))
        allocate(lchans (m_calc))
        allocate(mask   (m_calc))
        allocate(istore (m_calc,2))
        if (l_chk_rl) then
          allocate(reg_lim(n_lev  ,     2,n_prof_block))
        else
          allocate(reg_lim(0,0,0))
        end if
        if (l_chk_god) then
          allocate(rflag(n_chan,n_prof_block)) ; rflag = 0
        else
          allocate(rflag(0,0))
        end if

        if (trace_mem_time)                                               &
          call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': allocate arrays for rttov_k')

        if (s%gopts%lev_mode >= 1) then
          nlt = 0
        else
          call rtifc_check_nlevs(n_lev, nlt, ierr)
          if (ierr /= 0) call finish(proc//' rtifc_check_nlevs', trim(rtifc_errMsg(ierr)))
        end if
        n_lev_ = n_lev + nlt
        n_lay  = n_lev_ - 1

        allocate(quality(n_chan,n_prof_block))

        if (tsk == TSK_K) then

          ncv   = 2*n_lev_ + nsv + max(0, n_pc)

          allocate (temp_k (n_lev,n_chan,n_prof_block))
          allocate (humi_k (n_lev,n_chan,n_prof_block))
          allocate (t2m_k  (      n_chan,n_prof_block))
          allocate (q2m_k  (      n_chan,n_prof_block))
          allocate (psurf_k(      n_chan,n_prof_block))
          allocate (ctp_k  (      n_chan,n_prof_block))
          allocate (cfr_k  (      n_chan,n_prof_block))
          allocate (u10m_k (      n_chan,n_prof_block))
          allocate (v10m_k (      n_chan,n_prof_block))
          allocate (stemp_k(      n_chan,n_prof_block))
          allocate (dppmvdq(n_lev                    ))
          allocate (emiss_k(      n_chan,n_prof_block))
          allocate (Hnew   (      n_chan,ncv         ))
          allocate (ivld   (      n_chan             ))
        end if

        l_l2c = l_first .and. any(s% iopts(:)% l2c_type > 0)
        if (l_l2c) then
          allocate(l2c_aux(n_chan, n_prof_block))
        end if
        if (l_first) then
          allocate(sinfl(n_chan, n_prof_block)) ; sinfl = -1._sp
        end if

        l_btcs = any(btest(s%iopts(1:s%n_instr)%rad_out, OUT_CSB)) .and.  po_context == POC_FG
        if (l_btcs) then
          allocate(bt_fg_cs(n_chan,n_prof_block))
        else
          allocate(bt_fg_cs(n_chan,0))
        end if

        l_im = l_first .and. (s% im_rttov_indx > 0)
        if (l_im) then
          allocate(im_fg(s%n_im_chan,n_prof_block), im_emis(s%n_im_chan,n_prof_block))
          n_im_fg    = 0 ;  n_im_emis       = 0 ;  n_im_emis_f  = 0
          n_im_tskin = 0 ;  n_im_tskin_styp = 0 ;  n_im_tskin_f = 0
          l_im_all = .true.
        end if

        l_addinterp = (s%gopts%lev_mode > 0)

        if (wr_opd .and. l_addinterp) call finish('process_tovs_mult', &
             'wr_opdep currently not supported with interpolation (lev_mode > 0)')
        call write_opdep('alloc')

        !-----------------------------
        ! Loop over blocks of profiles
        !-----------------------------
        if (trace_mem_time) &
             call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': Loop over blocks of profiles')
        n_prof_aux = 0
!NEC$ nomove
        loop_prof_blocks: do while (n_prof_aux < n_prof)
          i_prof_start =     n_prof_aux + 1
          i_prof_end   = min(n_prof_aux + n_max_prof, n_prof)
          n_prof_block = i_prof_end - i_prof_start + 1

          npr = 0 ; ipr = 0
          do i = 1, n_prof_block
            i_prof = i_prof_start + i - 1
            ib =  r% i_box(i_prof)
            is =  r% i_reprt(i_prof)
            spt => obs% o(ib)% spot(is)
            if (ldeb(spt)) then
              npr = npr + 1
              ipr(npr) = i
              if (npr==1) write(usd,*) dpref,'RTTOV call context',po_context,po_ilns,po_lmon,tsk
              if (npr == size(ipr)) exit
            end if
          end do

FTRACE_BEGIN("tovs_mult_prof:fill")
          !--------------------------
          ! provide profiles to rttov
          !--------------------------
          call rtifc_fill_input(                          &
               status        = ierr,                      &
               rad           = r,                         &
               iopts         = s%rttov_indx(1:s%n_instr), &
               istart        = i_prof_start,              &
               iend          = i_prof_end,                &
               pe            = dace% pe,                  &
               wr_profs      = spt_hd_debug(1:n_wr_profs),&
               wr_profs_fmt  = wr_profs_fmt)
          if (ierr /= 0) call finish(proc//' rtifc_fill_input', &
               trim(rtifc_errMsg(ierr)))
FTRACE_END("tovs_mult_prof:fill")
          !----------------------
          ! loop over instruments
          !----------------------
!NEC$ nomove
          loop_instruments: do instr = 1, s%n_instr
            ioff = s%o_ch_i(instr)
            nch  = s%n_ch_i(instr)

            n_calc = count(r%valid(ioff+1:ioff+nch,i_prof_start:i_prof_end))
            if (n_calc <= 0) cycle

            n_lev_rt = s%iopts(instr)%rt_nlevs
            call rtifc_get_opts(iopt=s%rttov_indx(instr), addinterp=l_aux, addclouds=l_cloud)
            l_addinterp = l_addinterp .or. l_aux

            rad_out = ibset(s%iopts(instr)%rad_out, OUT_ASB)
            allocate(l_vis(nch))
            l_vis = .false.
            if (associated(s%iflag)) then
              l_vis(1:nch) = (s%iflag(ioff+1:ioff+nch) == 1)
            end if
            if (any(l_vis)) then
              rad_out = ibset(rad_out, OUT_VIS)
              allocate(radrefl(nch,n_prof_block),refl(nch,n_prof_block))
              do i = 1,nch
                if (l_vis(i)) then
                  refl(i,:) = r%emiss(ioff+i,i_prof_start:i_prof_end)
                else
                  refl(i,:) = 0._wp
                end if
              end do
            else
              allocate(radrefl(0,0),refl(0,0))
            end if

            l_l2c_i = l_l2c .and. (s% iopts(instr)% l2c_type > 0)
            if (l_l2c_i) then
              allocate (overc  (n_lay,nch,n_prof_block))
              if (s%iopts(1)%l2c_use_rad) then
                allocate (r_clear(    nch,n_prof_block))
                rad_out = ibset(rad_out, OUT_CSR)
              else
                allocate (r_clear(0,0))
              end if
            else
              allocate (overc(0,0,0), r_clear(0,0))
            end if

            ! Surface influence calculations necessary?
            mode = 0
            m = 1
            if (associated(s%band)) m = max(m, maxval(s%band(ioff+1:ioff+nch)))
            m = min(m,size(s%iopts(instr)%surf_infl_mode))
            do i = 1, m
              if (s%iopts(instr)%surf_infl_mode(i) > 0) &
                   mode = ior(mode, s%iopts(instr)%surf_infl_mode(i))
            end do
            l_surf_infl = l_first .and. mode /= 0 .and. (po_context == POC_FG)
                      !any(btest(s%flag(1:n_chan), USE_SURFINFL))

            l_pr_tr = (s%bc%n_tr > 0) .and. (po_context == POC_FG)
            if (l_pr_tr) then
              l_pr_tr = any(s%bc%ib_ch(ioff+1:ioff+nch) > 0)
            end if


            ! Allocate transm/opdep if required
            ! TODO: transm not required for VIS??
            if (  l_hgt                                              .or. &
                  wr_opd                                             .or. &
                 (l_surf_infl .and. btest(mode, SINFL_TRANSM))       .or. &
                 debug_opdep_chan >= 0                               .or. &
                 (.not.l_addinterp .and. (i_opdep_tsk >= 0))         .or. &
                 (.not.l_addinterp .and. l_l2c .and. l2c_god_corr)   .or. &
                 l_pr_tr                                                  ) then
              allocate(transm (n_lev  ,nch,n_prof_block))
              if (l_cloud) allocate(transmcld(n_lev  ,nch,n_prof_block))
            else
              allocate(transm(0,0,0))
            end if
            if (.not.allocated(transmcld)) allocate(transmcld(0,0,0))
            if (  wr_opd                                             .or. &
                 (l_addinterp .and. (i_opdep_tsk >= 0))              .or. &
                 (l_addinterp .and. l_l2c .and. l2c_god_corr)             ) then
              allocate(opdep(n_lev_rt-1,nch,n_prof_block))
              call rtifc_coef_prop(s%rttov_indx(instr), sign_opdep=sign_opdep)
            else
              allocate(opdep(0,0,0))
            endif

            if (s%iopts(instr)%do_lambertian) then
              specularity => r%specularity(ioff+1:ioff+nch,i_prof_start:i_prof_end)
            else
              specularity => spec_dum
            end if


            ub1_rf = min(ubound(rflag,1),ioff+nch)
            lb1_rf = lbound(rflag,1)
            if (ub1_rf >= ioff) lb1_rf = ioff + 1

FTRACE_BEGIN("tovs_mult_prof:setup_instr")

            !--------------------------------------------
            ! Set up the arrays of profiles and channels,
            ! for that RTTOV shall be run.
            !--------------------------------------------
            do k = 1, n_prof_block
              i_prof = i_prof_start + k - 1
              mask  ((k-1)*nch+1:k*nch)   = r%valid(ioff+1:ioff+nch,i_prof)
              lprofs((k-1)*nch+1:k*nch)   = k
              lchans((k-1)*nch+1:k*nch)   = s%chidx(ioff+1:ioff+nch)
              istore((k-1)*nch+1:k*nch,1) = (/ (l, l=1,nch) /)
            end do
            npc                = n_prof_block*nch
            n_calc = count(mask(1:npc))
            if (n_calc <= 0) cycle
            lprofs(1:n_calc)   = pack(lprofs(1:npc)   , mask(1:npc))
            lchans(1:n_calc)   = pack(lchans(1:npc)   , mask(1:npc))
            istore(1:n_calc,1) = pack(istore(1:npc, 1), mask(1:npc))
            istore(1:n_calc,2) = lprofs(1:n_calc)
            !-----------
            ! call RTTOV
            !-----------
FTRACE_END("tovs_mult_prof:setup_instr")

FTRACE_BEGIN("tovs_mult_prof:calc")
            n_calc_aux = 0
!NEC$ nomove
            calc_loop: do while (n_calc_aux < n_calc)
              i_start =     n_calc_aux + 1
              i_end   = min(n_calc_aux + n_max_calc, n_calc)

              npr_ = 0 ; ipr_ = 0
              if (npr > 0 .and. instr==1) then
                do k = 1, npr
                  if (any(lprofs(i_start:i_end) == ipr(k))) then
                    npr_ = npr_ + 1
                    ipr_(npr_) = ipr(k)
                  end if
                end do
              end if

              l_dealloc = rt_alloc_mode < 1 .or. &
                          (rt_alloc_mode == 1 .and. i_end==n_calc .and. i_prof_end==n_prof)

              select case(tsk)
              case(TSK_Y)
FTRACE_BEGIN("tovs_mult_prof:direct_ifc")
                call rtifc_direct(                                                     &
                     iopt       = s%rttov_indx(instr),                                 &
                     lprofs     = lprofs(i_start:i_end),                               &
                     chans      = lchans(i_start:i_end),                               &
                     emissiv    = r%emiss(ioff+1:ioff+nch,i_prof_start:i_prof_end),    &
                     T_b        = r%bt_fg(ioff+1:ioff+nch,i_prof_start:i_prof_end),    &
                     t_b_clear  = bt_fg_cs(ioff+1:ioff+nch,:),                         &
#if (_RTTOV_VERSION >= 13)
                     specularity= specularity,                                         &
                     radrefl    = radrefl,                                             &
                     refl       = refl,                                                &
                     quality    = quality(ioff+1:ioff+nch,:),                          &
                     transmcld  = transmcld(:,:,:),                                    &
#endif
                     status     = ierr,                                                &
                     radclear   = r_clear,                                             &
                     radovercast= overc,                                               &
                     istore     = istore(i_start:i_end,:),                             &
                     reg_lim    = reg_lim(:,:,:),                                      &
                     rflag      = rflag (lb1_rf:ub1_rf,:),                             &
                     transm     = transm(:,:,:),                                       &
                     opdep      = opdep (:,:,:),                                       &
                     dealloc    = l_dealloc,                                           &
                     rad_out_flg= rad_out,                                             &
                     iprint     = ipr_(1:npr_),                                        &
                     pe         = dace% pe,                                            &
                     l_pio      = .true.)
                if (ierr /= 0) call finish(proc//' rtifc_direct', &
                     trim (rtifc_errMsg (ierr)))
FTRACE_END("tovs_mult_prof:direct_ifc")
              case(TSK_K)
FTRACE_BEGIN("tovs_mult_prof:k_ifc")

                call rtifc_k(                                                          &
                     iopt       = s%rttov_indx(instr),                                 &
                     lprofs     = lprofs (i_start:i_end),                              &
                     chans      = lchans (i_start:i_end),                              &
                     emissiv    = r%emiss(  ioff+1:ioff+nch,i_prof_start:i_prof_end),  &
                     emissiv_k  = emiss_k(  ioff+1:ioff+nch,:),                        &
                     temp_k     = temp_k (:,ioff+1:ioff+nch,:),                        &
                     humi_k     = humi_k (:,ioff+1:ioff+nch,:),                        &
                     t2m_k      = t2m_k  (  ioff+1:ioff+nch,:),                        &
                     q2m_k      = q2m_k  (  ioff+1:ioff+nch,:),                        &
                     stemp_k    = stemp_k(  ioff+1:ioff+nch,:),                        &
                     T_b        = r%bt_fg(  ioff+1:ioff+nch,i_prof_start:i_prof_end),  &
                     t_b_clear  = bt_fg_cs( ioff+1:ioff+nch,:),                        &
#if (_RTTOV_VERSION >= 13)
                     specularity= specularity,                                         &
                     radrefl    = radrefl,                                             &
                     refl       = refl,                                                &
                     quality    = quality(  ioff+1:ioff+nch,:),                        &
                     transmcld  = transmcld(:,:,:),                                    &
#endif
                     status     = ierr,                                                &
                     radclear   = r_clear,                                             &
                     radovercast= overc,                                               &
                     transm     = transm (:,:,:),                                      &
                     opdep      = opdep  (:,:,:),                                      &
                     psurf_k    = psurf_k(  ioff+1:ioff+nch,:),                        &
                     u10m_k     = u10m_k (  ioff+1:ioff+nch,:),                        &
                     v10m_k     = v10m_k (  ioff+1:ioff+nch,:),                        &
                     ctp_k      = ctp_k  (  ioff+1:ioff+nch,:),                        &
                     cfraction_k= cfr_k  (  ioff+1:ioff+nch,:),                        &
                     istore     = istore (i_start:i_end,:),                            &
                     reg_lim    = reg_lim(:,:,:),                                      &
                     rflag      = rflag (lb1_rf:ub1_rf,:),                             &
                     dealloc    = l_dealloc,                                           &
                     rad_out_flg= rad_out,                                             &
                     iprint     = ipr_(1:npr_),                                        &
                     pe         = dace% pe,                                            &
                     l_pio      = .true.)

                if (ierr /= 0) call finish(proc//' rtifc_k', trim (rtifc_errMsg (ierr)))

                if (allocated(refl) .and. any(l_vis(1:nch))) then
                  do i = 1,nch
                    if (l_vis(i)) r%emiss(ioff+i,i_prof_start:i_prof_end) = refl(i,:)
                  end do
                end if

                FTRACE_END("tovs_mult_prof:k_ifc")
              end select
              n_calc_aux = i_end
            end do calc_loop
FTRACE_END("tovs_mult_prof:calc")

            ! for cloudy computations the total transmission per layer is the product of clear
            ! and cloudy transmission
            if (size(transmcld) > 0 .and. rttov_version >= 13) then
               transm = transmcld * transm
            end if

            if (any(l_vis(1:nch))) then
              do i = 1, n_prof_block
                i_prof = i_prof_start + i - 1
                do j = 1, nch
                  if (l_vis(j) .and. r% valid(ioff+j, i_prof)) then
                    r% bt_fg(ioff+j, i_prof) = radrefl(j,i)
                  end if
                end do
              end do
            end if

            call check_reg_lims

#if (_RTTOV_VERSION >= 13)
            if (qflag_bitmask > 0) then
              do i = 1, n_prof_block
                i_prof = i_prof_start + i - 1
               do j = 1, nch
                  if ( r% valid(j+ioff, i_prof) .and. &
                       iand(quality(j,i),qflag_bitmask) /= 0) then
                    l  = r% i_body (j+ioff,i_prof)
                    ib = r% i_box  (       i_prof)
                    call decr_use(obs% o(ib)% body(l)% use, check=CHK_OPERATOR, lflag=.true.)
                  end if
                end do
              end do
            end if
#endif

            if (npr > 0) then
              do m = 1, npr
                i = ipr(m)
                i_prof = i_prof_start + i - 1
                is =  r% i_reprt(i_prof)
                spt => obs% o(ib)% spot(is)
                if (ldeb(spt)) then
                  do j = ioff+1, ioff+nch
                    if (r% valid(j, i_prof)) then
                      write(usd,*) dpref,'rttov result bt_fg',s%chan(j),j,i,i_prof,r% bt_fg(j, i_prof),&
                           r% emiss(j, i_prof)
                      if (btest(rad_out, OUT_CSB)) write(usd,*) dpref,'rttov result bt_fg_cs',&
                           bt_fg_cs(j, i)
                      if (tsk == tsk_k) then
                        do k = 1, n_lev
                          write(usd,*) dpref,'t/q_k',j,k,temp_k(k,j,i),humi_k(k,j,i)
                          if (size(transm) > 0) write(usd,*) dpref,'transm',j,k,transm (k,j-ioff,i)
                        end do
                      end if
                    end if
                  end do
                end if
              end do
            end if

            !--------------
            ! store results
            !--------------
            if (present(y)) then
              FTRACE_BEGIN("tovs_mult_prof:store_y")
!NEC$ ivdep
              do i = 1, n_prof_block
                i_prof = i_prof_start + i - 1
                ib = r% i_box  (i_prof)
                is =  r% i_reprt(i_prof)
                spt => obs% o(ib)% spot(is)
                do j = ioff+1, ioff+nch
                  if (r% valid(j, i_prof)) then
                    k  = r% i_body(j,i_prof)
                    y% s(ib)% x(k) = r% bt_fg(j, i_prof)
                    if (ldeb(spt)) write(usd,*) dpref,' Y:',j,s%chan(j),r% bt_fg(j, i_prof)
                  end if
                end do
              end do
              FTRACE_END("tovs_mult_prof:store_y")
            endif

            if_l2c: if (l_l2c_i) then
              FTRACE_BEGIN("tovs_mult_prof:l2c")
              if (s%iopts(instr)%l2c_use_rad) then
                call lev2chan(r_clear (1:nch,                      1:n_prof_block), &
                              overc (:,ioff+1:ioff+nch,            1:n_prof_block), &
                              r% valid(ioff+1:ioff+nch, i_prof_start:i_prof_end),   &
                              l2c_aux( ioff+1:ioff+nch,            1:n_prof_block), &
                              s%iopts(instr)%l2c_type,                              &
                              rel_lim=s%iopts(instr)%l2c_rel_lim,                   &
                              abs_lim=s%iopts(instr)%l2c_abs_lim, ipr=ipr(1:npr))
              else
                call lev2chan(r%bt_fg (ioff+1:ioff+nch, i_prof_start:i_prof_end),   &
                              overc (:,ioff+1:ioff+nch,            1:n_prof_block), &
                              r% valid(ioff+1:ioff+nch, i_prof_start:i_prof_end),   &
                              l2c_aux( ioff+1:ioff+nch,            1:n_prof_block), &
                              s%iopts(instr)%l2c_type,                              &
                              rel_lim=s%iopts(instr)%l2c_rel_lim,                   &
                              abs_lim=s%iopts(instr)%l2c_abs_lim, ipr=ipr(1:npr))
              end if
              if (l2c_god_corr) then
                allocate(l2c_corr(nch, n_prof_block))
                if (.not. l_addinterp) then
                  call rtifc_l2c_god(s%rttov_indx(instr),                                 &
                                     s%chidx(ioff+1:ioff+nch),                            &
                                     r% valid(ioff+1:ioff+nch, i_prof_start:i_prof_end),  &
                                     l2c_aux (ioff+1:ioff+nch, 1:n_prof_block),           &
                                     l2c_corr(1:nch, 1:n_prof_block),                     &
                                     transm=transm(:,:,:),                                &
                                     ideb=ipr(1:npr))
                else
                  if (ubound(r%p,2) >= i_prof_end) then
                    ip0 = i_prof_start ; ip1 = i_prof_end
                  else
                    ip0 = 1            ; ip1 = 1
                  end if
                  call rtifc_l2c_god(s%rttov_indx(instr),                                 &
                                     s%chidx (ioff+1:ioff+nch),                           &
                                     r% valid(ioff+1:ioff+nch, i_prof_start:i_prof_end),  &
                                     l2c_aux (ioff+1:ioff+nch, 1:n_prof_block),           &
                                     l2c_corr(1:nch, 1:n_prof_block),                     &
                                     opdep=opdep(:,:,:),                                  &
                                     p_l2c=r%p(:,ip0:ip1),                                &
                                     ideb=ipr(1:npr))
                end if
                l2c_aux (ioff+1:ioff+nch, 1:n_prof_block) = l2c_corr(1:nch, 1:n_prof_block)
                deallocate(l2c_corr)
              end if
              FTRACE_END("tovs_mult_prof:l2c")
            end if if_l2c

            call calc_hgt
            call calc_surf_infl
            call calc_transm_pred

            deallocate(overc, r_clear)
            deallocate(l_vis  )
            deallocate(radrefl,refl)

            if (instr == 1) then
              call check_opdep
              call write_opdep('store')
            end if
            deallocate(opdep, transm, transmcld)

          end do loop_instruments

          call check_god
          call check_ps

          k_stuff: if (tsk == TSK_K) then
FTRACE_BEGIN("tovs_mult_prof:store_k")
            ib_prev = -1
            loop_prof: do i = 1, n_prof_block
              i_prof = i_prof_start + i - 1
              ib =  r% i_box(i_prof)
              if (ib /= ib_prev) then
                call obs_block (ob,obs,ib)
                ib_prev = ib
              end if
              is =  r% i_reprt(i_prof)
              spt => obs% o(ib)% spot(is)
              n_chan = count(r%valid(:,i_prof))
              ncv    = spt% i% n
              ivld(1:n_chan) = pack((/(k,k=1,s%n_chan)/),mask=r%valid(:,i_prof))

              !-----------------
              ! Convert humidity
              !-----------------
              if (rttov_version <= 10) then
                if (rt_humi == 1) then
                  ! ppmv over moist air
                  r% q_fg(:,i_prof) = q_ppmv_moist(r% q_fg(:,i_prof)) ! ppmv_moist -> q
                  r% q2m   (i_prof) = q_ppmv_moist(r% q2m   (i_prof)) ! ppmv_moist -> q
                  dppmvdq(:) = dppmv_moist_dq(r% q_fg(:,i_prof))      ! d(ppmv_moist)/dq
                  dppmvdq_2m = dppmv_moist_dq(r% q2m   (i_prof))      ! d(ppmv_moist)/dq
                  do j = 1,n_chan
                    humi_k (:,ivld(j),i) = humi_k (:,ivld(j),i) * dppmvdq(:)
                    q2m_k    (ivld(j),i) = q2m_k    (ivld(j),i) * dppmvdq_2m
                  end do
               elseif (rt_humi == 2) then
                  ! ppmv over dry air
                  r% q_fg(:,i_prof) = q_ppmv_dry(r% q_fg(:,i_prof)) ! ppmv_dry -> q
                  r% q2m   (i_prof) = q_ppmv_dry(r% q2m   (i_prof)) ! ppmv_dry -> q
                  dppmvdq(:) = dppmv_dry_dq(r% q_fg(:,i_prof))      ! d(ppmv_dry)/dq
                  dppmvdq_2m = dppmv_dry_dq(r% q2m   (i_prof))      ! d(ppmv_dry)/dq
                  do j = 1,n_chan
                    humi_k (:,ivld(j),i) = humi_k (:,ivld(j),i) * dppmvdq(:)
                    q2m_k    (ivld(j),i) = q2m_k    (ivld(j),i) * dppmvdq_2m
                  end do
                else
                  humi_k (:,ivld(1:n_chan),i) = humi_k (:,ivld(1:n_chan),i) * 1.d6/RDRD
                  q2m_k  (  ivld(1:n_chan),i) = q2m_k  (  ivld(1:n_chan),i) * 1.d6/RDRD
                end if
              else
#if (_RTTOV_VERSION >= 12)
                select case(default_gas_units)
                case(gas_unit_specconc)
                  ! Nothing to do
                case(gas_unit_ppmv)
                  ! ppmv over moist air
                  r% q_fg(:,i_prof) = q_ppmv_moist(r% q_fg(:,i_prof)) ! ppmv_moist -> q
                  r% q2m   (i_prof) = q_ppmv_moist(r% q2m   (i_prof)) ! ppmv_moist -> q
                  dppmvdq(:) = dppmv_moist_dq(r% q_fg(:,i_prof))      ! d(ppmv_moist)/dq
                  dppmvdq_2m = dppmv_moist_dq(r% q2m   (i_prof))      ! d(ppmv_moist)/dq
                  do j = 1,n_chan
                     humi_k (:,ivld(j),i) = humi_k (:,ivld(j),i) * dppmvdq(:)
                     q2m_k    (ivld(j),i) = q2m_k    (ivld(j),i) * dppmvdq_2m
                  end do
                case(gas_unit_ppmvdry)
                  ! ppmv over dry air
                  r% q_fg(:,i_prof) = q_ppmv_dry(r% q_fg(:,i_prof)) ! ppmv_dry -> q
                  r% q2m   (i_prof) = q_ppmv_dry(r% q2m   (i_prof)) ! ppmv_dry -> q
                  dppmvdq(:) = dppmv_dry_dq(r% q_fg(:,i_prof))      ! d(ppmv_dry)/dq
                  dppmvdq_2m = dppmv_dry_dq(r% q2m   (i_prof))      ! d(ppmv_dry)/dq
                  do j = 1,n_chan
                    humi_k (:,ivld(j),i) = humi_k (:,ivld(j),i) * dppmvdq(:)
                    q2m_k    (ivld(j),i) = q2m_k    (ivld(j),i) * dppmvdq_2m
                  end do
                end select
#endif
              end if

              !--------------------------
              ! Store result in H matrix.
              !--------------------------
              call K2H(H   = Hnew   (1:n_chan,1:ncv),             &
                   temp_k  = temp_k (:,ivld(1:n_chan),i),         &
                   humi_k  = humi_k (:,ivld(1:n_chan),i),         &
                   t2m_k   = t2m_k  (ivld(1:n_chan),i),           &
                   q2m_k   = q2m_k  (ivld(1:n_chan),i),           &
    !              psurf_k = psurf_k(ivld(1:n_chan),i) * 0.01_wp, &
                   psurf_k = psurf_k(ivld(1:n_chan),i) / 100._wp, & !same as rttov_mult_prof=F
                   stemp_k = stemp_k(ivld(1:n_chan),i),           &
                   ctp_k   = ctp_k  (ivld(1:n_chan),i),           &
                   cfr_k   = cfr_k  (ivld(1:n_chan),i),           &
                   emis_k  = emiss_k(ivld(1:n_chan),i),           &
                   tj      = tjs(iset)%a(i_prof),                 &
                   pz_bg   = spt% pz_bg,                          &
                   n_chan  = n_chan,                              &
                   spot    = spt)

              call destruct(tjs(iset)% a(i_prof))
              call insert_full_block_2(ob% H, Hnew(1:n_chan,1:ncv), &
                   spt% o% i, spt% i% i)
              !------------------------------------
              ! Store result of forward calculation
              !------------------------------------
              ob% yi% x(spt% o% i+1:spt% o% i+spt% o% n) = &
                   pack(r% bt_fg(:,i_prof), r%valid(:,i_prof))

            end do loop_prof
FTRACE_END("tovs_mult_prof:store_k")

            ! call calc_surf_infl

          end if k_stuff

          if (l_im) call calc_im_fg('calc')

          !!!! TODO: emissivity calculation is usually done once.
          !!!! TODO: thus, it is not necessary in all cases to load/store tovs here:

          FTRACE_BEGIN("tovs_mult_prof:l2c store")
          if (l_l2c .or. l_surf_infl .or. l_im .or. po_context == POC_FG) then
            allocate(emis(s%n_chan))
            if (l_btcs     ) allocate(bt_cs(s%n_chan))
            if (l_surf_infl) allocate(sinf (s%n_chan))
            if (l_l2c      ) allocate(l2c  (s%n_chan))
            if (l_im       ) allocate(im_ch(mx_imch,s%n_im_chan))
            loop_tovs_store: do i = 1, n_prof_block
              i_prof = i_prof_start + i - 1
              ib =  r% i_box(i_prof)
              is =  r% i_reprt(i_prof)
              spt => obs% o(ib)% spot(is)
              l_debug = ldeb(spt)
              if (l_debug) then
                call load(obs% o(ib), spt, ttovs, tovs_io=TTOVS_CI, emis=emis)
              else
                call load(obs% o(ib), spt, ttovs, tovs_io=0, emis=emis)
              end if
              if (l_l2c) then
                j = 0
                do k = 1, s%n_instr
                  if (s%iopts(k)%l2c_type > 0) then
                    ioff = s%o_ch_i(k)
                    nch  = s%n_ch_i(k)
                    m = count(r%valid(ioff+1:ioff+nch, i_prof))
                    l2c(j+1:j+m) = pack(l2c_aux(ioff+1:ioff+nch,i), &
                         mask=r%valid(ioff+1:ioff+nch,i_prof))
                    if (l_debug) then
                      do l = ioff+1,ioff+nch
                        write(usd,*) dpref,'store l2c:',j,s%chan(l),l2c_aux(l,i),&
                             r%valid(l,i_prof)
                      end do
                    end if
                    j = j + m
                  end if
                end do
                call store(obs% o(ib), spt, ttovs, tovs_io=0, l2c=l2c)
              end if
              !            ! Keep the sign of ttovs%emis in order to keep the information, whether the
              !            ! emissivity was calculated with FASTEM (see comment in prep_rttov_prof)
              !             emis(1:ttovs%nchan)  = sign(pack(r%emiss(:, i_prof), r%valid(:, i_prof)),&
              !                                         emis(1:ttovs%nchan))
              emis(1:ttovs%nchan)  = pack(r%emiss(:, i_prof), r%valid(:, i_prof))
              call store(obs% o(ib), spt, ttovs, tovs_io=0, emis=emis)
              if (l_debug) then
                do j = 1, ttovs% nchan
                  write(usd,*) dpref,'store emis:',j,s%chan(ttovs%ci(j)),emis(j)
                end do
              end if
              if (l_btcs) then
                bt_cs(1:ttovs%nchan) = pack(bt_fg_cs(:, i), r%valid(:, i_prof))
                call store(obs% o(ib), spt, ttovs, tovs_io=0, bt_cs=bt_cs)
                if (l_debug) then
                  do j = 1, ttovs% nchan
                    write(usd,*) 'debug_spot ',spt%hd%id,spt%id,' store bt_cs:',j,s%chan(ttovs%ci(j)),bt_cs(j)
                  end do
                end if
              end if
              if (l_surf_infl) then
                sinf(1:ttovs%nchan) = pack(sinfl (:, i), r%valid(:, i_prof))
                ttovs% sinfl => sinf(1:ttovs%nchan)
                call store(obs% o(ib), spt, ttovs, tovs_io=TTOVS_BASE, sinfl=sinf)
                if (l_debug) then
                  do j = 1, ttovs% nchan
                    write(usd,*) dpref,'store sinfl:',j,s%chan(ttovs%ci(j)),sinf(j)
                  end do
                end if
              end if
              if (l_im) then
                im_ch = 0._sp
                call get_im_ch_ind(ttovs, IMCH_MEAN_B, i_imch_v)
                if (i_imch_v <= 0) call finish(proc,'No IMCH_MEAN_B found in t_tovs%im_ch_v!')
                im_ch(i_imch_v,1:s%n_im_chan) = real(im_fg(1:s%n_im_chan,i), kind=sp)
                call get_im_ch_ind(ttovs, IMCH_EMIS, i_imch_v)
                if (i_imch_v <= 0) call finish(proc,'No IMCH_EMIS found in t_tovs%im_ch_v!')
                im_ch(i_imch_v,1:s%n_im_chan) = real(im_emis(1:s%n_im_chan,i), kind=sp)
                call store(obs% o(ib), spt, ttovs, tovs_io=TTOVS_IMCH, im_ch=im_ch)
                if (l_debug) then
                  do j = 1, s%n_im_chan
                    write(usd,*) dpref,'store im_fg  :',j,s%im_chan(j),im_fg(1:s%n_im_chan,i)
                    write(usd,*) dpref,'store im_emis:',j,s%im_chan(j),im_emis(1:s%n_im_chan,i)
                  end do
                end if
              end if

              call destruct(ttovs)
            end do loop_tovs_store
            deallocate(emis)
            if (l_btcs     ) deallocate(bt_cs)
            if (l_l2c      ) deallocate(l2c)
            if (l_surf_infl) deallocate(sinf)
            if (l_im       ) deallocate(im_ch)
          end if
          FTRACE_END("tovs_mult_prof:l2c store")

          n_prof_aux = i_prof_end

        end do loop_prof_blocks

        !------------------
        ! deallocate memory
        !------------------
        if (trace_mem_time) call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': deallocate')
        if (l_im) deallocate(im_fg, im_emis)
        deallocate(quality)
        deallocate(lprofs )
        deallocate(lchans )
        deallocate(mask   )
        deallocate(istore )
        deallocate(reg_lim)
        deallocate(rflag  )
        deallocate(bt_fg_cs)
        if (l_l2c) deallocate(l2c_aux)
        if (allocated(sinfl)) deallocate(sinfl)
        if (trace_mem_time) call stop_time ('process_tovs_mult '//trim(descr_call)//'_'//char3(iset)//': deallocate k')
        if (tsk == TSK_K) then
          deallocate(temp_k )
          deallocate(humi_k )
          deallocate(t2m_k  )
          deallocate(q2m_k  )
          deallocate(ctp_k  )
          deallocate(cfr_k  )
          deallocate(psurf_k)
          deallocate(u10m_k )
          deallocate(v10m_k )
          deallocate(stemp_k)
          deallocate(Hnew   )
          deallocate(ivld   )
          deallocate(emiss_k)
          deallocate(dppmvdq)
          if (associated(tjs(iset)% a)) deallocate(tjs(iset)% a)
        end if

      endif n_prof_gt_0

      call write_opdep('write')

    end do loop_datasets

    call flush_buf

    call calc_im_fg('stat')

    l_first = .false.

    if (trace_mem_time) call stop_time ('process_tovs_mult '//trim(descr_call)//&
         ': end do loop_datasets: rtifc_dealloc_profiles')
    if (rt_alloc_mode <= 2) call rtifc_cleanup(lprof=.true., lcoef=.false., latlas=.false.)

    call destruct(rad(:))

    if (po_context /= -99 .and. .not.po_lmon) po_context_old = po_context
    if (po_ilns  /= -99 .and. .not.po_lmon) po_ilns_old  = po_ilns
    po_lmon_old = po_lmon

!   end if
FTRACE_END("tovs_mult_prof")

    if (trace_mem_time) call stop_time ('process_tovs_mult '//trim(descr_call)//&
         ': end subroutine process_tovs_mult')

  contains


    subroutine calc_im_fg(task)
      character(len=*),  intent(in)  :: task
      type(t_tskin_opt), pointer     :: tso
      real(wp)                       :: ts(1)
      real(wp),          allocatable :: tskin(:)
      integer                        :: styp
      integer,           pointer     :: im_tso(:)
      logical                        :: l_ts, l_found

      if (task == 'calc') then

        if (.not.l_im) return

        nch    = s%n_im_chan
        n_calc = n_prof_block*nch
        if (n_calc <= 0) return

        do k = 1, n_prof_block
          lprofs((k-1)*nch+1:k*nch)   = k
          lchans((k-1)*nch+1:k*nch)   = s%im_chan_indx(1:nch)
          istore((k-1)*nch+1:k*nch,1) = (/ (l, l=1,nch) /)
        end do
        istore(1:n_calc,2) = lprofs(1:n_calc)
        im_emis = -1._wp

#if (_RTTOV_VERSION >= 12)
        ! Emissivity atlas
        if (s% gopts% im_atlas_id >= 0) then
          do k = 1, n_prof_block
            i_prof = i_prof_start + k - 1
            if (any(r%stype(1,i_prof) == (/rts_land,rts_ice/))) then
              call rtifc_emis_atlas(s% im_rttov_indx, s% gopts% im_atlas_id, (/(k,j=1,nch)/), s%im_chan_indx(1:nch), &
                   im_emis(1:nch,k), ierr, max_dst=40._wp, ldebug=any(ipr(1:npr)==k))
              if (ierr /= 0) call finish(proc//':calc_im_fg rtifc_emis_atlas', trim (rtifc_errMsg (ierr)))
              n_im_emis = n_im_emis + nch
              where (im_emis(1:nch,k) <= 0._wp .or. im_emis(1:nch,k) >= 1._wp) im_emis(1:nch,k) = -1._wp
              n_im_emis_f = n_im_emis_f + count(im_emis(1:nch,k) == -1._wp)
            end if
          end do
        end if
        ! Tskin retrieval

        ! TODO: cloud cover threshold


        im_tso => s%gopts%im_tskin_opt
        if (all(im_tso(:) >= 0)) then
          if (.not.btest(s%gopts%opt_vars, OPTV_CLD_FRC)) call finish(proc, &
               'calc_im_fg: imager cloud fraction not available for imager T_skin retrieval.')
          allocate(tskin(n_prof_block))
          l_found = .false.
          tso => null()
          do i = 1, s%n_tskin_opt
            tso => s%tskin_opt(i)
            if (tso%inst == im_tso(1) .and. tso%n_cdyn > 0 .and. associated(tso%cdyn)) then
              if (tso%cdyn(1) == im_tso(2)) exit
            end if
            tso => null()
          end do
          tskin_instr_loop: do i=1,s%n_instr
            if (s%instr(i) == im_tso(1)) then
              do j=s%o_ch_i(i)+1,s%o_ch_i(i)+s%n_ch_i(i)
                if (s%chan(j) == im_tso(2)) then
                  do k = 1, n_prof_block
                    i_prof = i_prof_start + k - 1
                    l  = r% i_body (j,i_prof)
                    ib = r% i_box  (  i_prof)
                    tskin(k) = r%ts_fg(i_prof)
                    styp     = r%stype(1,i_prof)
                    if ( iand(2**styp,im_tso(3)) /= 0 .and. r%valid(j,i_prof) .and. &
                         r%cfraction(i_prof) <= s%gopts%im_tskin_cfrac) then
                      n_im_tskin = n_im_tskin + 1
                      if (associated(tso)) then
                        if (any(tso%styp == styp)) then
                          ! the required Tskin retrieval was done before in mo_tskin for the main instrument
                          if (any(ipr(1:npr)==k)) write(usd,*) 'imag_tskin_ret recycle result'
                          CYCLE
                        end if
                      end if
                      call rtifc_tskin_retrieve(s%rttov_indx(i),(/k/), s%chan(j:j), s%chidx(j:j), &
                           real(obs% o(ib)% body(l:l)%o,kind=wp), reshape((/1._wp/),(/1,1/)), r%emiss(j:j,i_prof), &
                           ts(:), ierr, l_ts, pe=dace%pe, ldeb=any(ipr(1:npr)==k))
                      if (ierr /= 0) call finish(proc//':calc_im_fg rtifc_tskin_retrieve', trim (rtifc_errMsg (ierr)))
                      if (any(ipr(1:npr)==k)) write(usd,*) 'imag_tskin_ret',ierr,l_ts,ts(:),obs%o(ib)%body(l)%o
                      if (l_ts) then
                        tskin(k) = ts(1)
                      else
                        n_im_tskin_f = n_im_tskin_f + 1
                      end if
                      if (styp >= 0 .and. styp <= n_styp) n_im_tskin_styp(styp) = n_im_tskin_styp(styp) + 1
                    end if
                  end do
                  l_found = .true.
                  exit tskin_instr_loop
                end if
              end do
            end if
          end do tskin_instr_loop
          if (.not.l_found) call finish(proc//':calc_im_fg','instr/chan not found.')
        else
          allocate(tskin(0))
        end if
#endif
        if (.not.allocated(tskin)) allocate(tskin(0))

        n_calc_aux = 0
!NEC$ nomove
        calc_loop: do while (n_calc_aux < n_calc)
          i_start =     n_calc_aux + 1
          i_end   = min(n_calc_aux + n_max_calc, n_calc)

          npr_ = 0 ; ipr_ = 0
          do k = 1, npr
            if (any(lprofs(i_start:i_end) == ipr(k))) then
              npr_ = npr_ + 1
              ipr_(npr_) = ipr(k)
            end if
          end do

          l_dealloc = rt_alloc_mode < 1 .or. &
               (rt_alloc_mode == 1 .and. i_end==n_calc .and. i_prof_end==n_prof)

          call rtifc_direct(iopt    = s% im_rttov_indx,              &
                            lprofs  = lprofs (i_start:i_end),        &
                            chans   = lchans (i_start:i_end),        &
                            emissiv = im_emis(1:nch,1:n_prof_block), &
#if (_RTTOV_VERSION >= 13)
                            tskin   = tskin,                         &
#endif
                            istore  = istore (i_start:i_end,:),      &
                            T_b     = im_fg  (1:nch,1:n_prof_block), &
                            status  = ierr,                          &
                            iprint  = ipr_(1:npr_),                  &
                            pe      = dace%pe,                       &
                            dealloc = l_dealloc,                     &
                            l_pio   = .true.)
          if (ierr /= 0) call finish(proc//':calc_im_fg rtifc_direct', trim (rtifc_errMsg (ierr)))

          n_calc_aux = i_end
        end do calc_loop

        n_im_fg = n_im_fg + n_calc

      elseif (task == 'stat') then

        l_im_all = p_or(l_im_all)
        if (.not.l_im_all) return

        n_im_fg         = p_sum(n_im_fg)
        n_im_emis       = p_sum(n_im_emis)
        n_im_emis_f     = p_sum(n_im_emis_f)
        n_im_tskin      = p_sum(n_im_tskin)
        n_im_tskin_f    = p_sum(n_im_tskin_f)
        n_im_tskin_styp = p_sum(n_im_tskin_styp)
        if (dace%lpio) then
          write(6,*)
          write(6,'(1x,A)') 'Colocated Imager first guess calculation statistics'
          write(6,'(3x,A,T30,":",I9)') '#RTTOV calcs.:',n_im_fg
          write(6,'(3x,A,T30,":",I9)') '#emiss. atlas calcs.',n_im_emis
          if (n_im_emis > 0) &
               write(6,'(3x,A,T30,":",I9,2x,F7.3,"%")') '#emiss. atlas calls failed',n_im_emis_f,(100.*n_im_emis_f)/n_im_emis
          i = sum(n_im_tskin_styp)
          write(6,'(3x,A,T30,":",I9)') '#Tskin retrievals required',n_im_tskin
          if (n_im_tskin > 0) then
            write(6,'(3x,A,T30,":",I9)') '#Tskin retrievals recycled',n_im_tskin - i
            write(6,'(3x,A,T30,":",I9)') '#Tskin retrievals done',i
            if (i > 0) &
                 write(6,'(3x,A,T30,":",I9,2x,F7.3,"%")') '#Tskin retrievals failed',n_im_tskin_f,(100.*n_im_tskin_f)/i
            write(6,'(T30," ",3(I9,1x,A,3x))') &
                 n_im_tskin_styp(rts_sea),'sea',n_im_tskin_styp(rts_land),'land',n_im_tskin_styp(rts_ice),'ice'
          end if
        end if
      else
        call finish(proc, 'calc_im_fg: unknown task')
      end if
    end subroutine calc_im_fg


    subroutine prep_check(opts, l_chk, act)
      type(t_rtstat_opts), intent(in),  target  :: opts
      logical,             intent(out)          :: l_chk
      integer,             intent(out), pointer :: act

      if (.not.po_lmon) then
        select case(po_context)
        case(POC_FG)
          act => opts%action_fg
        case(POC_ANA)
          act => opts%action_ana
        case(POC_TOP)
          act => opts%action_mon
        case(1:size(opts%action))
          act => opts%action(po_context)
        case default
          act => null()
        end select
      else
        act => opts%action_mon
      end if
      if (associated(act)) then
        l_chk = (iand(act,7) /= 0)
        if (po_context == POC_FG) then
          ! Perform test for fg also if tsk_k not set in opts%tsk
          l_chk = l_chk .and. (iand(tsk, opts%tsk+tsk_k) /= 0)
        else
          l_chk = l_chk .and. (iand(tsk, opts%tsk) /= 0)
        end if
      else
        l_chk = .false.
      end if
    end subroutine prep_check


    subroutine prep_checks
#if (_RTTOV_VERSION > 0)
      ! Prepare opdep check
      call prep_check(chk_opdep, l_chk_opdep, act_opdep)
      l_chk_opdep = l_chk_opdep .and. (po_context >= 1 .or. po_context == POC_TOP)
      if (l_chk_opdep) then
        if (po_context_old /= po_context .or. &
            po_ilns_old     > po_ilns    .or. &
            (po_lmon.neqv.po_lmon_old)   .or. &
            .not.associated(keep_rtstat)) then
          i_opdep_tsk = 0
          if (.not.associated(keep_rtstat)) then
            allocate(keep_rtstat(n_set))
          end if
        else
          i_opdep_tsk = 1
        end if
      else
        i_opdep_tsk = -1
      end if

      ! Prepare reg_lims check
      call prep_check(chk_reg_lims, l_chk_rl, act_rl)
      if (l_chk_rl) then
        chk_reg_lims_ifc = act_rl
        l_rl_var = btest(act_rl, 4)  ! distinguish variables
        l_rl_lev = btest(act_rl, 5)  ! distinguish levels
        l_rl_dir = l_rl_var .and. l_rl_lev .and. btest(act_rl, 6)  ! distinguish directions (upper/lower limit)
        if (po_context >= 1) then
          if (po_context /= po_context_old) then
            nout_chk_reg_lims = nout_chk_reg_lims + 1
            l_new_out = .true.
          else
            l_new_out = .false.
          end if
          ni = chk_reg_lims%nitout
        else
          ni = 0
        end if
        if (ni > 1) then
          ! We have to track the results of multiple RTTOV calls.
          if (.not.associated(keep_rtstat)) allocate(keep_rtstat(n_set))
          ! Prepare logical array to store information to be kept:
          nrl_lev = 1
          if (l_rl_lev) nrl_lev = n_lev
          nrl_var = 1
          if (l_rl_var) nrl_var = 2
          nrl_dir = 1
          if (l_rl_dir) nrl_dir = 2
          n1 = nrl_lev * nrl_var * nrl_dir
          allocate(l_rl(n1,chk_reg_lims%nitout))
        end if
        if (fg_prof_top == 0) then
          if (interp_strato /= 0) chk_plim_t = p_top
          chk_plim_q = hum_ana_top
        end if
      else
        chk_reg_lims_ifc = 0
      end if

      ! Prepare ps check
      call prep_check(chk_ps, l_chk_ps, act_ps)
      l_chk_ps = l_chk_ps .and. po_context >= 1
      if (l_chk_ps) then
        if (po_context_old /= po_context .or. &
            po_ilns_old     > po_ilns    .or. &
            (po_lmon.neqv.po_lmon_old)   .or. &
            .not.associated(keep_rtstat)) then
          i_ps_tsk = 0
          if (.not.associated(keep_rtstat)) then
            allocate(keep_rtstat(n_set))
          end if
        else
          i_ps_tsk = 1
        end if
      else
        i_ps_tsk = -1
      end if
#endif

#if (_RTTOV_VERSION >= 12)
      ! Prepare god check
      call prep_check(chk_god, l_chk_god, act_god)
      if (associated(act_god)) then
        god_thresh_ifc = god_thresh(max(po_context,1))
        chk_god_ifc = act_god                                ! btest(chk_god,1) fully done in mo_rtifc
        l_chk_god = l_chk_god .and. &                        ! btest(chk_god,1) .or. btest(chk_god,2)
                    (btest(act_god,1) .or. btest(act_god,2)) ! require work here
      else
        chk_god_ifc = 0
      end if
      if (l_chk_god) then
        if (chk_god%nitout > 1 .and. .not.associated(keep_rtstat) .and. po_context >= 1) then
          allocate(keep_rtstat(n_set))
        end if
        if (po_context >= 1) then
          if (po_context /= po_context_old) then
            nout_chk_god = nout_chk_god + 1
            l_new_out = .true.
          else
            l_new_out = .false.
          end if
          ni = chk_god%nitout
        else
          ni = 0
        end if
        if (ni > 1) then
          ! We have to track the results of multiple RTTOV calls.
          if (.not.associated(keep_rtstat)) allocate(keep_rtstat(n_set))
        end if
      end if
#endif

    end subroutine prep_checks

    subroutine check_god
      integer :: i, j, m, ni
      logical :: lpr
      logical, allocatable :: l_god(:,:)

      check_god_stuff: if (l_chk_god) then
        if (po_context >= 1) then
          ni = chk_god%nitout
        else
          ni = 0
        end if

        loop_god: do i = 1, n_prof_block
          i_prof = i_prof_start + i - 1

          if (ni > 1) then
            ! We have to track the results of multiple RTTOV calls
            allocate(l_god(s%n_chan, ni))
            call load(krt%spt(i_prof), larr4=l_god)
            if (l_new_out) then
              l_god(:,1:ni-1) = l_god(:,2:ni)
              l_god(:,ni)     = .false.
            end if
            l_god(:,ni) = l_god(:,ni) .or. (rflag(:,i) /= 0)
          else
            allocate(l_god(s%n_chan, 1))
            l_god(:,1) = (rflag(:,i) /= 0)
          end if

          ! Debug
          ib = r% i_box  (  i_prof)
          is = r% i_reprt(  i_prof)
          spt => obs% o(ib)% spot(is)
          if (ldeb(spt)) then
            write(usd,*) dpref,'chk_god_deb',rflag(:,i)
            write(usd,*) dpref,l_god(:,ubound(l_god,2))
          end if
          ! end Debug

          do j = 1, s%n_chan
            if (.not.r% valid(j, i_prof)) cycle

            if (ldeb(spt)) write(usd,*) dpref,'chk_god_deb',i,j,s%chan(j),rflag(j,i),l_god(j,:)

            if (all(l_god(j,:))) then
              l  = r% i_body (j,i_prof)
              ib = r% i_box  (  i_prof)
              is = r% i_reprt(  i_prof)
              spt => obs% o(ib)% spot(is)
              if (chk_god%op_na_bit >= 0) then
                l_set = btest(obs% o(ib)% body(l)% op_na,chk_god%op_na_bit)
              else
                l_set = .false.
              end if
              if (btest(act_god, 1)) then
                lpr = .false.
                if (po_context == POC_FG) then
                  call decr_use(obs% o(ib)% body(l)% use, check=CHK_OPERATOR, state=STAT_PAS_REJ, lflag=.true.)
                  lpr = .true.
                else
                  obs% o(ib)% body(l)% use% flags = &
                       ibset(obs% o(ib)% body(l)% use% flags, CHK_OPERATOR)
                  if (chk_god%op_na_bit >= 0) then
                    if (.not.l_set) then
                      obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na,chk_god%op_na_bit)
                      lpr = .true.
                    end if
                  end if
                end if
                if (lpr .and. btest(act_god, 0)) then
                  write(msg,*) 'Discard RAD/TOVS obs ',spt%hd%id,s%satid,s%grid,s%chan(j),&
                       ' (god) op_na=',obs% o(ib)% body(l)% op_na
                  call add_line(trim(msg))
                end if
              end if
              if (btest(act_god, 2)) then
                lpr = .false.
                if (po_context == POC_FG) then
                  call decr_tovs_use(spt, obs% o(ib), CHK_OPERATOR, state=STAT_PAS_REJ)
                  lpr = .true.
                else
                  spt%use%flags = ibset(spt%use%flags, CHK_OPERATOR)
                  if (chk_god%op_na_bit >= 0) then
                    do m = 1, s%n_chan
                      if (r% valid(m, i_prof)) then
                        l  = r% i_body (m,i_prof)
                        if (.not.l_set) then
                          obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na,chk_god%op_na_bit)
                          lpr = .true.
                        end if
                      end if
                    end do
                  end if
                end if
                if (lpr) then
                  write(msg,*) 'Discard RAD/TOVS spot ',spt%hd%id,s%satid,s%grid,' (god) op_na=',obs% o(ib)% body(l)% op_na
                  call add_line(trim(msg))
                end if
              end if
              call debug_spot(sp=spt, obs=obs%o(ib), hint='check_god')
            end if
          end do

          if (ni > 1) then
            call store(krt%spt(i_prof), larr4=l_god)
          end if
          deallocate(l_god)
        end do loop_god

      end if check_god_stuff
    end subroutine check_god

    subroutine check_reg_lims
      integer :: i, j, ilev, ivar, idir, ni
      logical :: l_dum, l_exceed

      reg_lims_stuff: if (l_chk_rl .and. s%grid == s%instr(instr)) then
        if (po_context >= 1) then
          ni = chk_reg_lims%nitout
        else
          ni = 0
        end if

        loop_reg_lims: do i = 1, n_prof_block
          i_prof = i_prof_start + i - 1
          ! Comprise information as required
          if (.not.l_rl_dir) where(reg_lim(:,:,i) /= 0) reg_lim(:,:,i) = 1
          if (.not.l_rl_var) where(reg_lim(:,1,i) /= 0 .and. reg_lim(:,2,i) /= 0) reg_lim(:,1,i) = 1
          if (.not.l_rl_lev) then
            do j = 1,nrl_var
              if (any(reg_lim(:,j,i) /= 0)) reg_lim(1,j,i) = 1
            end do
          end if
          if (ni > 1) then
            ! We have to track the results of multiple RTTOV calls
            call load(krt%spt(i_prof), larr2=l_rl)
            if (l_new_out) l_rl(:,1:ni-1) = l_rl(:,2:ni) ! move old results
            l_exceed = .false.
            j = 0
            do idir = -1, -1 + 2 * (nrl_dir-1), 2 ! -1 or -1, 1
              do ivar = 1, nrl_var
                do ilev = 1, nrl_lev
                  j = j + 1
                  l_dum = (reg_lim(ilev,ivar,i) /= 0)
                  if (l_dum .and. l_rl_dir) l_dum = (sign(reg_lim(ilev,ivar,i), 1) == idir)
                  if (.not.l_new_out) then
                    ! "aggregate" results within current outer loop
                    l_dum = l_dum .or. l_rl(j,ni)
                  end if
                  l_rl(j, ni) = l_dum
                  if (all(l_rl(j,:))) l_exceed = .true.
                end do
              end do
            end do
            call store(krt%spt(i_prof), larr2=l_rl)
            if (nout_chk_reg_lims < ni) cycle loop_reg_lims
          else
            l_exceed = any(reg_lim(:,:,i) /= 0)
          endif
          if (l_exceed) then
            ib = r% i_box  (  i_prof)
            is = r% i_reprt(  i_prof)
            spt => obs% o(ib)% spot(is)
            if (btest(act_rl, 0)) then
              ! Nothing to do, printing done in mo_rtifc
            end if
            if (btest(act_rl, 1).or.btest(act_rl, 2)) then
              if (po_context == POC_FG) then
                write(msg,*) 'Discard RAD/TOVS spot ',spt%hd%id,s%satid,s%grid,' (reg_limits)'
                call add_line(trim(msg))
                call decr_tovs_use(spt, obs% o(ib), CHK_OPERATOR, state=STAT_PAS_REJ)
              else
                l_set = .false.
                do j = 1, s%n_chan
                  if (r% valid(j, i_prof)) then
                    l  = r% i_body (j,i_prof)
                    obs% o(ib)% body(l)% use% flags = &
                         ibset(obs% o(ib)% body(l)% use% flags, CHK_OPERATOR)
                    if (chk_reg_lims%op_na_bit >= 0) then
                      if (.not.btest(obs% o(ib)% body(l)% op_na,chk_reg_lims%op_na_bit)) then
                        obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na, chk_reg_lims%op_na_bit)
                        l_set = .true.
                      end if
                    end if
                  end if
                end do
                if (l_set) then
                  write(msg,*) 'Discard RAD/TOVS spot ',spt%hd%id,s%satid,s%grid,' (reg_limits) op_na=',obs% o(ib)% body(l)% op_na
                  call add_line(trim(msg))
                end if
                if (btest(act_rl, 2)) then
                  spt%use%flags = ibset(spt%use%flags, CHK_OPERATOR)
                end if
              end if
            end if
          end if
        end do loop_reg_lims
      end if reg_lims_stuff
    end subroutine check_reg_lims


    subroutine check_opdep
      logical,  allocatable :: l_opdep(:,:), l_opdep_old(:,:), mask_(:)
      real(wp), allocatable :: od(:), trans_l(:)
      real(wp)              :: min_transm, max_transm, od_tot, tr
      integer               :: nl, nch
      integer               :: i, j, k, l, m, ii
      logical               :: l_set
      real(wp), parameter   :: small = 1.E-30

      if (instr > 1 .or. ioff > 0) then
        call finish('process_tovs_mult', 'invalid call to check_opdep')
      end if

      ! optical depth: tau  !NEGATIVE!
      ! transmission : T    ! T_i level i to space, T_(i-1) > T_i
      ! incident rad.: ri
      ! outgoing rad.: ro             ro_i = ri_(i-1)
      ! tau = ln(T) = ln(ro/ri) = - ln(ri/ro)
      ! tau_i = ln(ro_i/ri_i) = ln(ri_(i-1)/ri_i) = ln(T_i/T_(i-1))
      ! exp(tau_i) = T_i/T_(i-1) < exp(tau_threshold)
      !              T_(i-1)/T_i > exp(-tau_threshold)
      !              T_(i-1) > T_i * exp(-tau_threshold)

      if (debug_opdep_chan > 0) then
        do i = 1, n_prof_block
          i_prof = i_prof_start + i - 1
          is = r% i_reprt(  i_prof)
          ib = r% i_box  (  i_prof)
          spt => obs% o(ib)% spot(is)
          if (ldeb(spt)) then
            if (l_addinterp) then
              nl = n_lev_rt - 1
            else
              nl = n_lev-1
            end if
            do j = 1, s%n_chan
              k = (j-1)*nl
              if (btest(s%flag(j),USE_PASSIVE) .or. .not.r%valid(j,i_prof)) cycle
              if (s%chan(j) /= debug_opdep_chan) cycle
              do l = 1, nl
                if (transm(l,j,i) > small .and. transm(l+1,j,i) > small) then
                  write(msg,*) 'opdep',j,s%chan(j),l+nlt,&
                       transm(l+1,j,i)/transm(l,j,i),log(transm(l+1,j,i)/transm(l,j,i)),&
                       transm(l:l+1,j,i)
                else
                  write(msg,*) 'opdep',j,s%chan(j),l+nlt,'*','*',&
                       transm(l:l+1,j,i)
                end if
                write(usd,*) dpref,trim(msg)
              end do
            end do
          end if
        end do
      end if

      opdep_stuff: if (i_opdep_tsk >= 0) then
        min_transm = opdep_min_transm(max(po_context,1))
        max_transm = 1._wp - opdep_min_transm(max(po_context,1))

        nch = s%n_ch_i(1)
        if (l_addinterp) then
          nl = n_lev_rt - 1
        else
          nl = n_lev - 1
          allocate(od(nl), trans_l(nl))
        end if

        allocate(l_opdep(nch*nl,3))

        if (i_opdep_tsk == 0 .and. wr_opd_disc) then
          if (associated(keep_disc)) deallocate(keep_disc)
          n_keep_disc = 0
        end if

        do i = 1, n_prof_block
          i_prof = i_prof_start + i - 1

          is = r% i_reprt(  i_prof)
          ib = r% i_box  (  i_prof)
          spt => obs% o(ib)% spot(is)
          do j = 1, nch
            k = (j-1)*nl
            if (r% valid(j, i_prof)) then
              if (l_addinterp) then
                l_opdep(k+1:k+nl,1) = (sign_opdep * opdep(1:nl,j,i) > 0._wp)
                l_opdep(k+1:k+nl,2) = (sign_opdep * opdep(1:nl,j,i) >  abs(-opdep_thresh(max(po_context,1))))
                od_tot = 0._wp
                do m = 1, nl
                  tr = exp(-od_tot)
                  od_tot = od_tot + max(sign_opdep*opdep(m,j,i), 0._wp)
                  l_opdep(k+m,3) = tr > min_transm .and. tr < max_transm
                end do
              else
                !                         T_(i-1)              T_i
                l_opdep(k+1:k+nl,1) = (transm(1:nl,j,i) > transm(2:nl+1,j,i))
                l_opdep(k+1:k+nl,2) = (transm(1:nl,j,i) > transm(2:nl+1,j,i) * exp(-opdep_thresh(max(po_context,1))))
                l_opdep(k+1:k+nl,3) = transm(1:nl,j,i) > min_transm .and. transm(1:nl,j,i) < max_transm
              end if
            else
              l_opdep(k+1:k+nl,:) = .false.
            end if
          end do
          if (i_opdep_tsk == 1) then
            allocate(l_opdep_old(nch*nl,3),mask_(nl))
            call load(krt%spt(i_prof), larr1=l_opdep_old)
            do j = 1, nch
              k = (j-1)*nl
              if (ldeb(spt)) then
                do l = 1, nl
                  write(usd,*) dpref,'opdep',j,s%chan(j),l+nlt,l_opdep_old(k+l,:)
                end do
              end if
              if (r% valid(j, i_prof)) then
                mask_ = (l_opdep(k+1:k+nl,1) .neqv. l_opdep_old(k+1:k+nl,1)) .and. &
                        (l_opdep(k+1:k+nl,2) .or.   l_opdep_old(k+1:k+nl,2)) .and. &
                        (l_opdep(k+1:k+nl,3) .or.   l_opdep_old(k+1:k+nl,3))
                if (any(mask_)) then
                  l  = r% i_body (j,i_prof)
                  ib = r% i_box  (  i_prof)
                  is = r% i_reprt(  i_prof)
                  spt => obs% o(ib)% spot(is)
                  if (chk_opdep%op_na_bit >= 0) then
                    l_set = btest(obs% o(ib)% body(l)% op_na,chk_opdep%op_na_bit)
                  else
                    l_set = .false.
                  end if
                  if (btest(act_opdep, 0) .and. .not.l_set) then
                    if (l_addinterp) then
                      od_tot = 0._wp
                      do ii = 1, nl
                        tr = exp(-od_tot)
                        if (mask_(ii)) then
                          write(msg,'(" Discontinuity in opdep",I10,1x,I3.3,1x,I3.3,1x,I5.5,&
                               &" in levels ",I3.2,1x,F10.5,2(1x,E13.6),1x,F7.2,"N",F7.2,"E")') &
                               spt%hd%id,s%satid,s%grid,s%chan(j),ii,             &
                               opdep(ii,j,i),exp(-od_tot),exp(-(od_tot + abs(opdep(ii,j,i)))),&
                               spt%col%c%dlat,spt%col%c%dlon
                          call add_line(trim(msg))
                        end if
                        od_tot = od_tot + abs(opdep(ii,j,i))
                      end do
                    else
                      where(transm(1:nl,j,i) > small)
                        trans_l(:) = transm(2:nl+1,j,i)/transm(1:nl,j,i)
                      elsewhere(transm(2:nl+1,j,i) > 0._wp)
                        trans_l(:) = huge(0._wp)
                      elsewhere
                        trans_l(:) = 1._wp
                      end where
                      where(transm(1:nl,j,i) > small .and. transm(2:nl+1,j,i) > small)
                        od  (:) = log(transm(2:nl+1,j,i)/transm(1:nl,j,i))
                      elsewhere
                        od  (:) = 0._wp
                      end where
                      do ii = 1, nl
                        if (mask_(ii)) then
                          write(msg,'(" Discontinuity in opdep",I10,1x,I3.3,1x,I3.3,1x,I5.5,&
                               &" in levels ",I3.2,1x,F9.7,1x,F10.5,2(1x,E13.6),1x,F7.2,"N",F7.2,"E")')&
                               spt%hd%id,s%satid,s%grid,s%chan(j),ii+nlt,trans_l(ii),&
                               od(ii),transm(ii,j,i),transm(ii+1,j,i),&
                               spt%col%c%dlat,spt%col%c%dlon
                          call add_line(trim(msg))
                        end if
                      end do
                    end if
                  end if
                  if (btest(act_opdep, 1)) then
                    obs% o(ib)% body(l)% use% flags = &
                         ibset(obs% o(ib)% body(l)% use% flags, CHK_OPERATOR)
                    if (chk_opdep%op_na_bit >= 0) then
                      if (.not.l_set) then
                        obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na,chk_opdep%op_na_bit)
                        write(msg,*) 'Discard RAD/TOVS obs ',spt%hd%id,s%satid,s%grid,s%chan(j),&
                             ' (opdep) op_na=',obs% o(ib)% body(l)% op_na
                        call add_line(trim(msg))
                      end if
                    end if
                  end if
                  if (btest(act_opdep, 2)) then
                    spt%use%flags = ibset(spt%use%flags, CHK_OPERATOR)
                    if (chk_opdep%op_na_bit >= 0) then
                      l_set = .false.
                      do m = 1, nch
                        if (r% valid(m, i_prof)) then
                          l  = r% i_body (m,i_prof)
                          if (.not.l_set) then
                            l_set = .true.
                            obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na,chk_opdep%op_na_bit)
                          end if
                        end if
                      end do
                      if (l_set) then
                        write(msg,*) 'Discard RAD/TOVS spot ',spt%hd%id,s%satid,s%grid,&
                             ' (opdep) op_na=',obs% o(ib)% body(l)% op_na
                        call add_line(trim(msg))
                      end if
                    end if
                  end if

                  if (wr_opd_disc) then
                    ! Keep information about discontinuites in optical depths
                    do ii = 1, nl
                      if (mask_(ii)) call add_disc(iset, spt%hd%id, j, ii)
                    end do
                  end if

                end if
              end if
            end do
            deallocate(l_opdep_old,mask_)
          else
            call store(krt%spt(i_prof), larr1=l_opdep)
          end if
        end do
        deallocate(l_opdep)
      end if opdep_stuff

    end subroutine check_opdep


    subroutine add_disc(iset, id, chan, lev)
      integer, intent(in) :: iset
      integer, intent(in) :: id
      integer, intent(in) :: chan
      integer, intent(in) :: lev

      type(t_keep_disc), pointer :: kdp(:) => null()
      integer :: i, m_keep_disc

      do i = 1, n_keep_disc
        if ( iset == keep_disc(i)%iset .and. &
             id   == keep_disc(i)%id   .and. &
             chan == keep_disc(i)%chan .and. &
             lev  == keep_disc(i)%lev  ) RETURN
      end do
      m_keep_disc = 0
      if (associated(keep_disc)) m_keep_disc = size(keep_disc)
      if (n_keep_disc >= m_keep_disc) then
        if (associated(keep_disc)) then
          kdp => keep_disc
          allocate(keep_disc(m_keep_disc + n_alloc_keep_disc))
          keep_disc(1:n_keep_disc) = kdp(1:n_keep_disc)
          deallocate(kdp)
        else
          allocate(keep_disc(n_alloc_keep_disc))
        end if
      end if
      n_keep_disc = n_keep_disc + 1
      keep_disc(n_keep_disc) = t_keep_disc(iset, id, chan, lev)
    end subroutine add_disc

    elemental logical function l_wr_opd_chan(instr,chan,band) result(lwr)
      integer, intent(in) :: instr
      integer, intent(in) :: chan
      integer, intent(in) :: band

      select case(instr)
      case(10)
        lwr = any(chan == (/ 9, 10, 11 /))
      case(15)
        lwr = any(chan == (/ 3, 4, 5 /))
      case(16)
        lwr = (band == 3)
      case(19)
        lwr = (chan >= 18)
      case(21)
        lwr = any(chan == (/ 2, 3 /))
      case(22)
        lwr = any(chan == (/ 2 /))
      case(34)
        lwr = .true.
      case(44)
        lwr = any(chan == (/ 2, 3, 4 /))
      case(56)
        lwr = any(chan == (/ 2, 3, 4 /))
      case(71)
        lwr = (chan >= 10)
      case(73)
        lwr = (chan >= 10)
      case default
        lwr = .false.
      end select
    end function l_wr_opd_chan


    subroutine write_opdep(task)
      character(len=*), intent(in) :: task
#if (_RTTOV_VERSION >= 12)
      character(len=300)          :: fname
      type(t_wr_opd), allocatable :: wo_rcv(:)
      integer,        allocatable :: index_snd(:)
      integer                     :: n_snd(dace% npe), n_rec(dace% npe)
      real(wp)                    :: opdep_
      integer                     :: iunit
      integer                     :: n_wo_all    ! Total number of entries to be written
      integer                     :: n_wr        ! Number of entries written to file
      integer                     :: nf_wo       ! Number of files to be written
      integer                     :: nf_pe       ! Number of files to be written by a single pe
      integer                     :: pe_wo       ! Output pe
      integer                     :: i_prof, i, j, k, i_p
      logical                     :: l_new, l_exist

      if (.not.wr_opd) return

      select case(task)
      case('alloc')
        if (n_prof <= 0) then
          n_wo = 0
          allocate(wo_snd(n_wo))
        else
          ! Count entries to be written
          n_wo = 0
          do i_prof = 1, n_prof
            if (wr_opd_disc .and. (po_context <= POC_TOP .or. po_lmon)) then
              ib = r% i_box  (i_prof)
              is =  r% i_reprt(i_prof)
              spt => obs% o(ib)% spot(is)
            end if
            do j = 1, s%n_ch_i(1)
              if (btest(s%flag(j),USE_PASSIVE) .and. .not.btest(wr_opdep, 5)) cycle
              if (.not.l_wr_opd_chan(s%instr(1), s%chan(j), s%band(j))) cycle
              if (.not.r% valid(j, i_prof)) cycle
              do k = 20, n_lev -1
                if (wr_opd_disc .and. (po_context == POC_NN .or. po_lmon)) then
                  if (.not.any(keep_disc(1:n_keep_disc)%iset == iset      .and. &
                               keep_disc(1:n_keep_disc)%id   == spt%hd%id .and. &
                               keep_disc(1:n_keep_disc)%chan == j         .and. &
                               keep_disc(1:n_keep_disc)%lev  == k         ) ) cycle
                end if
                n_wo = n_wo + 1
              end do
            end do
          end do
          allocate(wo_snd(n_wo))
          i_wo = 0
        end if

      case('store')
        if (n_wo > 0) then
          do i = 1, n_prof_block
            i_prof = i_prof_start + i - 1
            ib = r% i_box  (i_prof)
            is =  r% i_reprt(i_prof)
            spt => obs% o(ib)% spot(is)
            i_p = min(i_prof, size(r%p, 2))
            do j = 1, nch
              if (btest(s%flag(j),USE_PASSIVE) .and. .not.btest(wr_opdep, 5)) cycle
              if (.not.l_wr_opd_chan(s%instr(1), s%chan(j), s%band(j))) cycle
              if (.not.r% valid(j, i_prof)) cycle
              do k = 20, n_lev -1
                if (wr_opd_disc .and. (po_context == POC_NN .or. po_lmon)) then
                  if (.not.any(keep_disc(1:n_keep_disc)%iset == iset      .and. &
                               keep_disc(1:n_keep_disc)%id   == spt%hd%id .and. &
                               keep_disc(1:n_keep_disc)%chan == j         .and. &
                               keep_disc(1:n_keep_disc)%lev  == k         ) ) cycle
                end if
                i_wo = i_wo + 1
                if (transm(k+1,j,i) > 0._wp .and. transm(k,j,i) > 0._wp) then
                  opdep_ = log(transm(k+1,j,i)/transm(k,j,i))
                else
                  opdep_ = 0._wp
                end if
                wo_snd(i_wo)%chan   = j
                wo_snd(i_wo)%lev    = k + nlt
                wo_snd(i_wo)%spt    = spt%hd%id
                wo_snd(i_wo)%opdep  = opdep(k,j,i)
                wo_snd(i_wo)%opdep_ = opdep_
                wo_snd(i_wo)%transm = transm(k:k+1,j,i)
                wo_snd(i_wo)%q      = q_ppmv_dry(0.5_wp * (ppmv_dry_q(r% q_fg(k,i_prof)) + ppmv_dry_q(r% q_fg(k+1,i_prof))))
                wo_snd(i_wo)%t      = 0.5_wp * (r% t_fg(k,i_prof) + r% t_fg(k+1,i_prof))
                wo_snd(i_wo)%rh     = rh_q(wo_snd(i_wo)%q, wo_snd(i_wo)%t, 0.5_wp * (r%p(k,i_p)+r%p(k+1,i_p)))
!                if (any(spt_hd_debug==spt%hd%id) .and. s%chan(j)==3661 .and. k+nlt==31) &
!                   print*,'cp_opd',dace% pe,i,j,k,i_wo,opdep(k,j,i),opdep_,transm(k:k+1,j,i)
              end do
            end do
          end do
        end if

      case('write')
        n_wo_all = p_sum(n_wo)
        if (n_wo_all > 0) then
          ! Determine number of files (nf_wo) and number of files per pe (nf_pe)
          nf_wo = 0
          do j = 1, s%n_ch_i(1)
            if (btest(s%flag(j),USE_PASSIVE) .and. .not.btest(wr_opdep, 5)) cycle
            if (.not.l_wr_opd_chan(s%instr(1), s%chan(j), s%band(j))) cycle
            do k = 20, n_lev -1
              nf_wo = nf_wo + 1
            end do
          end do
          nf_pe = (nf_wo/dace% npe) + 1
          ! Determine number of entries to be sent
          allocate(index_snd(n_wo))
          n_snd(:) = 0
          nf_wo = 0
          i_wo = 0
          do j = 1, s%n_ch_i(1)
            if (btest(s%flag(j),USE_PASSIVE) .and. .not.btest(wr_opdep, 5)) cycle
            if (.not.l_wr_opd_chan(s%instr(1), s%chan(j), s%band(j))) cycle
            do k = 20, n_lev -1
              pe_wo = nf_wo/nf_pe + 1
              nf_wo = nf_wo + 1
              do i = 1, n_wo
                if (wo_snd(i)%chan==j .and. wo_snd(i)%lev == k+nlt) then
                  n_snd(pe_wo) = n_snd(pe_wo) + 1
                  i_wo = i_wo + 1
                  index_snd(i_wo) = i
!                  if (any(spt_hd_debug == wo_snd(i)%spt) .and. &
!                       wo_snd(i)%chan==307 .and. wo_snd(i)%lev+nlt==31)              &
!                       print*,'snd_wo',dace% pe,pe_wo,i,i_wo,index_snd(i_wo),wo_snd(i)%opdep
                end if
              end do
            end do
          end do
          call p_alltoall(n_snd, n_rec)
          n_wo = sum(n_rec)
          allocate(wo_rcv(n_wo))
          call p_alltoall(wo_snd(index_snd(:)),wo_rcv(:),sendcounts=n_snd,recvcounts=n_rec)
          ! Write out optical depths stuff
          do i = 1, n_wo
            if (i == 1) then
              l_new = .true.
              iunit = get_unit_number()
            else
              l_new = (wo_rcv(i)%chan /= wo_rcv(i-1)%chan) .or. &
                      (wo_rcv(i)%lev  /= wo_rcv(i-1)%lev )
              if (l_new) then
                !write(*,'(1x,I4,1x,A,I9,A)') dace% pe,' Wrote ',n_wr,' entries to opdep file '//trim(fname)
                close(iunit)
              end if
            end if
            if (l_new) then
              write(fname,'("opdep_sat",I3.3,"_instr",I3.3,"_chan",I4.4,"_lev",I3.3)') &
                   s%satid, s%instr(1), s%chan(wo_rcv(i)%chan), wo_rcv(i)%lev
              fname = path_file(aux,fname)
              inquire(file=fname, exist=l_exist)
              open(iunit,file=fname, position='append')
              if (.not.l_exist) then
                write(iunit,'(" # ",A,T15,I3.3)') 'satid',s%satid
                write(iunit,'(" # ",A,T15,I3.3)') 'instr',s%instr(1)
                write(iunit,'(" # ",A,T15,I5.5)') 'chan ',s%chan(wo_rcv(i)%chan)
                write(iunit,'(" # ",A,T15,I3.3)') 'level',wo_rcv(i)%lev
                write(iunit,'(" # ",A,T15,A   )') 'date ',cyyyymmddhh (ana_time)
                write(iunit,'(5(1x,A22),3(1x,A10),3(1x,A22))') &
                     '#                opdep', '       effective opdep', '    eff. transmittance', &
                     '     top transmittance', '  bottom transmittance', &
                     '      spot', '     phase', '      ilns', &
                     '           temperature', '     relative humidity', '     specific humidity'
              end if
              n_wr = 0
            end if
            if (any(spt_hd_debug == wo_rcv(i)%spt)) write(usd,*) 'debug_spot wr_opd',&
                 wo_rcv(i)%spt,trim(fname),po_ilns,wo_rcv(i)%opdep, wo_rcv(i)%opdep_, &
                 exp(wo_rcv(i)%opdep_), wo_rcv(i)%transm(:), wo_rcv(i)%t, wo_rcv(i)%rh, wo_rcv(i)%q
            write(iunit,'(5(1x,e22.15),3(1x,i10),3(1x,e22.15))') &
                 wo_rcv(i)%opdep, wo_rcv(i)%opdep_, &
                 exp(wo_rcv(i)%opdep_), wo_rcv(i)%transm(:), &
                 wo_rcv(i)%spt, po_context, po_ilns, wo_rcv(i)%t, wo_rcv(i)%rh, wo_rcv(i)%q
            n_wr = n_wr + 1
            if (i == n_wo) then
              !write(*,'(1x,I4,1x,A,I9,A)') dace% pe,' Wrote ',n_wr,' entries to opdep file '//trim(fname)
              close(iunit)
              call return_unit_number(iunit)
            end if
          end do
          deallocate(wo_rcv, index_snd)
        end if
        deallocate(wo_snd)
      end select
#endif
    end subroutine write_opdep


    subroutine check_ps
      logical :: l_ps(10,1), l_ps_old(10,1)  ! 10 Bits to store level number
      integer :: i, j, m,  nl, ilev, ilev_old
      ps_stuff: if (i_ps_tsk >= 0) then
        do i = 1, n_prof_block
          i_prof = i_prof_start + i - 1
          ib = r% i_box  (  i_prof)
          is = r% i_reprt(  i_prof)
          spt => obs% o(ib)% spot(is)
          nl = size(r%p,1)
          ilev = nl + 1
          do j = 1, nl
            if (r%ps_fg(i_prof) < r% p(j,1)) then
              ilev = j
              exit
            end if
          end do
          do j = 1,size(l_ps,1)
            l_ps(j,1) = btest(ilev,j-1)
          end do

          if (i_ps_tsk == 1) then
            call load(krt%spt(i_prof), larr3=l_ps_old)

            ilev_old = 0
            do j = 1,size(l_ps_old,1)
              if (l_ps_old(j,1)) ilev_old = ibset(ilev_old,j-1)
            end do
            if (any(l_ps .neqv. l_ps_old)) then
              ib = r% i_box  (  i_prof)
              is = r% i_reprt(  i_prof)
              spt => obs% o(ib)% spot(is)
              if (btest(act_ps, 0)) then
                write(msg,*) 'Discontinuity due to ps-bug in RTTOV',spt%hd%id,s%satid,s%grid
                call add_line(trim(msg))
              end if
              if (btest(act_ps, 1) .or. btest(act_ps, 2)) then
                do j = 1, s%n_chan
                  if (r% valid(j, i_prof)) then
                    l  = r% i_body (j,i_prof)
                    obs% o(ib)% body(l)% use% flags = &
                         ibset(obs% o(ib)% body(l)% use% flags, CHK_OPERATOR)
                    if (chk_ps%op_na_bit >= 0) then
                      if (.not.btest(obs% o(ib)% body(l)% op_na,chk_ps%op_na_bit)) then
                        obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na,chk_ps%op_na_bit)
                        write(msg,*) 'Discard RAD/TOVS obs ',spt%hd%id,s%satid,s%grid,s%chan(j),&
                             ' (ps) op_na=',obs% o(ib)% body(l)% op_na
                        call add_line(trim(msg))
                      end if
                    end if
                  end if
                end do
              end if
              if (btest(act_ps, 2)) then
                spt%use%flags = ibset(spt%use%flags, CHK_OPERATOR)
                if (chk_ps%op_na_bit >= 0) then
                  l_set = .false.
                  do m = 1, s%n_chan
                    if (r% valid(m, i_prof)) then
                      l  = r% i_body (m,i_prof)
                      if (.not.btest(obs% o(ib)% body(l)% op_na,chk_ps%op_na_bit)) then
                        l_set = .true.
                        obs% o(ib)% body(l)% op_na = ibset(obs% o(ib)% body(l)% op_na,chk_ps%op_na_bit)
                      end if
                    end if
                  end do
                  if (l_set) then
                    write(msg,*) 'Discard RAD/TOVS spot ',spt%hd%id,s%satid,s%grid,' (ps) op_na=',obs% o(ib)% body(l)% op_na
                    call add_line(trim(msg))
                  end if
                end if
              end if
            end if
          else
            call store(krt%spt(i_prof), larr3=l_ps)
          end if
        end do
      end if ps_stuff
    end subroutine check_ps


    subroutine calc_surf_infl
      real(wp), pointer :: p(:)
      real(wp)          :: d_stemp, d_emiss, mx_infl
      real(wp)          :: surf_infl
      real(wp)          :: s1,s2,w1,w2
      real(wp)          :: t_sum !, q_sum
      integer           :: band
      integer           :: i,j,k,l,m,instr,itr
      logical           :: chan_mask(m_chan)

      if (.not.l_surf_infl) return

      FTRACE_BEGIN("tovs_mult_prof:surf_infl")
      do i = 1, n_prof_block
        i_prof = i_prof_start + i - 1
        ib =  r% i_box(i_prof)
        is =  r% i_reprt(i_prof)
        spt => obs% o(ib)% spot(is)
        m = min(i_prof, size(r%p, 2))
        p => r%p(:,m)
        k = 0
        chan_mask(1:n_chan) = .false.
        instr=1
        do j = ioff+1,ioff+nch
          if (j > s%o_ch_i(instr) + s%n_ch_i(instr)) instr = instr + 1
          band = max(1, s%band(j))
          mode    = s%iopts(instr)%surf_infl_mode(band)
          mx_infl = s%iopts(instr)%max_surf_infl (band)
          d_stemp = s%iopts(instr)%d_stemp       (band)
          d_emiss = s%iopts(instr)%d_emiss       (band)
          if (r%valid(j,i_prof)) then
            k = k + 1
            if (btest(mode, SINFL_ABS)) then
              if (tsk /= TSK_K) call finish('calc_surf_infl@process_tovs_mult',&
                   'surf_infl_mode only supported with TSK_K')
              surf_infl = abs(stemp_k(j,i) * d_stemp) + &
                   abs(emiss_k(j,i) * d_emiss)
              obs% o(ib)% body(spt%o%i+k)% obs_par(1) = surf_infl
              obs% o(ib)% body(spt%o%i+k)% obs_par(2) = stemp_k(j,i)
            elseif (btest(mode, SINFL_REL)) then
              if (tsk /= TSK_K) call finish('calc_surf_infl@process_tovs_mult',&
                   'surf_infl_mode only supported with TSK_K')
              t_sum = sum(abs(temp_k(:,j,i)), mask=(preshPa(:) < r%ps_fg(i_prof)))
              ! q_sum = sum(abs(humi_k(:,j,i)), mask=(preshPa(:) < r%ps_fg(i_prof)))
              if (t_sum > 0._wp .and. any(preshPa(:) < r%ps_fg(i_prof))) then
                surf_infl = abs(stemp_k(j,i) / t_sum)
              else
                surf_infl = huge(surf_infl)
              end if
              obs% o(ib)% body(spt%o%i+k)% obs_par(1) = surf_infl
              obs% o(ib)% body(spt%o%i+k)% obs_par(2) = stemp_k(j,i)
            elseif (btest(mode, SINFL_TRANSM)) then
              itr = j-ioff
              ! Determine surface RTTOV level/layer
              l = -1
              do m = 1, n_lev
                if (p(m) > r%ps_fg(i_prof)) then
                  l = m
                  exit
                end if
              end do
              if (l == -1) l = n_lev + 1
              if (l > 1 .and. l <= n_lev) then
                w2 = abs(p(l-1) - r%ps_fg(i_prof))
                s1 = log(transm(l-1,itr,i))
                w1 = abs(p(l) - r%ps_fg(i_prof))
                s2 = log(transm(l,itr,i))
              elseif (l==1) then
                ! should not happen
                w1 = 0._wp
                s1 = 0._wp
                w2 = 1._wp
                s2 = log(transm(l,itr,i))
              elseif (l == n_lev + 1) then
                w1 = 1._wp
                s1 = log(transm(l-1,itr,i))
                w2 = 0._wp
                s2 = 0._wp
              end if
              surf_infl = exp((w1 * s1 + w2 * s2) / (w1 + w2))
              obs% o(ib)% body(spt%o%i+k)% obs_par(1) = surf_infl
            end if
            sinfl(j,i) = real(surf_infl, kind=sp)
            if (ldeb(spt)) write(usd,*) dpref,'surf_infl:',j,s%chan(j),surf_infl,&
                 obs% o(ib)% body(spt%o%i+k)% obs_par(1),mx_infl,mode,&
                 btest(s%flag(j), USE_SURFINFL)
            if (btest(s%flag(j), USE_SURFINFL)) then
              if (surf_infl > mx_infl) then
                if (ldeb(spt)) write(usd,*) dpref,'surf_infl_flag:',j,s%chan(j),spt%stlsf,mode,&
                     'land',.not.(btest(mode, SINFL_LAND).and..not.btest(spt%stlsf, SUR_LAND)),&
                     btest(mode, SINFL_LAND),btest(spt%stlsf, SUR_LAND), &
                     'highland',.not.(btest(mode, SINFL_HIGHLAND).and..not.btest(spt%stlsf, SUR_HIGHLAND)),&
                     btest(mode, SINFL_HIGHLAND),btest(spt%stlsf, SUR_HIGHLAND)
                if (.not.(btest(mode, SINFL_LAND).and..not.btest(spt%stlsf, SUR_LAND)).and.&
                     .not.(btest(mode, SINFL_HIGHLAND).and..not.btest(spt%stlsf, SUR_HIGHLAND))) then
                  chan_mask(j) = .true.
                end if
                ! Equivalent to:
                ! if (btest(mode, SINFL_LAND).and..not.btest(spt%stlsf, SUR_LAND)) then
                !   continue
                ! elseif (btest(mode, SINFL_HIGHLAND).and..not.btest(spt%stlsf, SUR_HIGHLAND)) then
                !   continue
                ! else
                !   chan_mask(j) = .true.
                ! end if
              end if
            end if
          end if
        end do
        if (any(chan_mask(1:n_chan))) call decr_tovs_use(spt, obs%o(ib), CHK_SURF, set=s, &
             state=STAT_PAS_REJ, tovs_flag=TF_SURF_INFL, chan_mask=chan_mask(1:n_chan))

        call debug_spot(sp=spt, obs=obs%o(ib), hint='surfinfl')
      end do
      FTRACE_END("tovs_mult_prof:surf_infl")
    end subroutine calc_surf_infl


    subroutine calc_hgt
      real(wp),      parameter   :: small = 1.E-30
      real(kind=wp)              :: w, wz, z_
      real(kind=wp)              :: sw, swz, swz2
      real(kind=wp)              :: plev, plev_wd
      real(kind=wp)              :: t_, t_prev, p_, p_prev, p_s
      real(kind=wp)              :: p(3), t(3)
                                    ! p(1): height where transm == 1.-pl_wd_thresh
                                    ! p(2): height where transm == 0.5
                                    ! p(3): height where transm == pl_wd_thresh
      integer                    :: i, j, k, l, i_p, itr
      real(kind=wp), parameter   :: p_inv = -huge(0._wp)

      if (.not.l_hgt) return

      do i = 1, n_prof_block
        i_prof = i_prof_start + i - 1
        is = r% i_reprt(  i_prof)
        ib = r% i_box  (  i_prof)
        spt => obs% o(ib)% spot(is)

        i_p = min(i_prof, ubound(r%p,2))

        do j = ioff+1, ioff+nch
          itr = j - ioff
          if (.not.r% valid(j, i_prof)) cycle
          if (l_vis(j-ioff) .and. pl_vis > 0._wp) then
            l  = r% i_body (j,i_prof)
            obs%o(ib)%body(l)%plev  = pl_vis
            cycle
          end if
          select case(pl_method)
          case(5,6)
            sw   = 0._wp
            swz  = 0._wp
            swz2 = 0._wp
            do l = 1, n_lev-1
              if (transm(l,itr,i) > small .and. transm(l+1,itr,i) > small) then
                select case(pl_method)
                case(5)
                  w  = -transm(l,itr,i) * log(transm(l+1,itr,i)/transm(l,itr,i))
                case(6)
                  w  = transm(l,itr,i) - transm(l+1,itr,i)
                end select
                if (pl_log) then
                  z_ = 0.5_wp * (log(r%p(l,i_p)) + log(r%p(l+1,i_p)))
                else
                  z_ = 0.5_wp * (r%p(l,i_p) + r%p(l+1,i_p))
                end if
                wz = w * z_
                sw  = sw   + w
                swz = swz  + wz
                swz2 =swz2 + wz * z_
              end if
            end do
            if (sw > 0._wp .and. swz > 0._wp) then
              swz  =       swz  / sw
              swz2 = sqrt (max(swz2 / sw - swz * swz,0._wp))
              if (pl_log) then
                swz = exp (swz)
              else
                swz2 = swz2 / swz
              endif
              plev    = swz
              plev_wd = swz2
            else
              plev    = 0.000001_wp
              plev_wd = 0._wp
            end if
          case(7)
            p(:) = p_inv
            t(1) = 1._wp - pl_wd_thresh
            t(2) = 0.5_wp
            t(3) = pl_wd_thresh
            t_prev = 1._wp ! previous transmission
            p_prev = 0.0001_wp
            if (pl_log) then
              p_prev = log(p_prev)
            end if
            lev_loop: do l = 1, n_lev
              t_ = transm(l,itr,i)
              if (pl_log) then
                p_ = log(r%p(l,i_p))
              else
                p_ = r%p(l,i_p)
              end if
              do k = 1, size(p)
                if (t_ <= t(k) .and. t_prev >= t(k)) then
                  if (abs(t_ - t_prev) > 1.E-10) then
                    p(k) = (p_prev * (t_-t(k)) + p_ * (t(k)-t_prev)) / (t_ - t_prev)
                  else
                    p(k) = 0.5_wp * (p_ + p_prev)
                  end if
                  if (all(p(:) >= 0._wp)) exit lev_loop
                end if
              end do
              t_prev = t_
              p_prev = p_
            end do lev_loop
            if (any(p(:) == p_inv)) then
              ! p_s = obs%o(ib)%lev(spt%i%i+2*n_lev+2)
              p_s = log(r%p(n_lev,i_p))
              if (.not.pl_log) p_s = exp(p_s)
              where(p(:) == p_inv) p(:) = p_s
            end if
            plev    = p(2)
            plev_wd = 0.5_wp * (p(3) - p(1))
            if (pl_log) then
              plev = exp(plev)
            else
              if (plev > 1.E-10) plev_wd = plev_wd / plev
            end if
          end select

          ! also make sure plev is in Pa for feedback file output, so if pressure unit
          ! is hPa here, convert to Pa
          if(trim(r% p_unit) == 'hpa') then
             plev = plev * 100._wp
          end if

          ! for the local model a minimum plevel might be set to avoid
          !  observations in the layer where relaxation against driving model
          !  takes place
          plev = max(plev,plev_min)

          l  = r% i_body (j,i_prof)
          obs%o(ib)%body(l)%plev       = plev
          obs%o(ib)%body(l)%plev_width = plev_wd

          if (ldeb(spt)) write(usd,*) dpref,'hgt/width',j,obs%o(ib)%body(l)%plev,&
               obs%o(ib)%body(l)%plev_width
        end do
      end do

    end subroutine calc_hgt


    subroutine calc_transm_pred
      integer               :: n_tr, np
      integer               :: i,j,k,it,ic,l
      integer               :: j_tt  ! index of channel in t_tovs
      integer               :: i_p, kp(2)
      logical               :: l_deb
      real(wp)              :: tr(2), w, t_, t1, t2, p1, p2, d_tr, wsum
      real(wp)              :: tr_surf, w_surf, T_surf
      integer               :: k_surf
      logical,  allocatable :: l_fail(:)
      real(sp), allocatable :: tt_tr(:,:)

      if (.not.l_pr_tr) return

      n_tr = s%bc%n_tr
      allocate(tt_tr(n_chan, n_tr),l_fail(n_chan))

      do i = 1, n_prof_block
        i_prof = i_prof_start + i - 1
        is = r% i_reprt(  i_prof)
        ib = r% i_box  (  i_prof)
        spt => obs% o(ib)% spot(is)
        i_p = min(i_prof, ubound(r%p,2))
        l_deb = ldeb(spt)

        call load(obs% o(ib), spt, ttovs, tovs_io=TTOVS_TR, tr=tt_tr)

        if (any(s%bc%type(1:n_tr) == 4)) then
          ! Determine lowest level above surface
          k_surf = -1
          do k = n_lev, 1, -1
            if (r%p(k,i_p) <= r%ps_fg(i_prof)) then
              k_surf = k
              exit
            end if
          end do
          if (k_surf < 1 .or. k_surf > n_lev) then
            write(0,*) 'spt=',spt%hd%id,'k_surf',k_surf,' ps=',r%ps_fg(i_prof),' p=',r%p(:,i_p)
            call finish('calc_transm_pred','failed to find surface level')
          end if
          ! Calc. atm. T at surface
          if (k_surf == 1) then
            ! Extrapolate linearly to surface
            w_surf = (log(r%ps_fg(i_prof)) - log(r%p(k_surf,i_p))) / &
                     (log(r%p(k_surf,i_p)) - log(r%p(k_surf+1,i_p)))
            T_surf = r%t_fg(k_surf,i_prof) + w_surf * &
                    (r%t_fg(k_surf,i_prof) - r%t_fg(k_surf+1,i_prof))
          else
            ! Interpolate
            w_surf = (log(r%ps_fg(i_prof))   - log(r%p(k_surf,i_p))) / &
                     (log(r%p(k_surf-1,i_p)) - log(r%p(k_surf,i_p)))
            T_surf = w_surf * r%t_fg(k_surf-1,i_prof) + (1._wp-w_surf) * r%t_fg(k_surf,i_prof)
          end if
          if (l_deb) write(usd,*) dpref,' transm_pred surf',r%ps_fg(i_prof),&
               k_surf,r%p(k_surf,i_p),w_surf,T_surf
        end if


        l_fail = .false.
        do j = ioff+1, ioff+nch
          ic = j - ioff
          if (.not.r% valid(j, i_prof)) cycle
          j_tt = count(r% valid(1:j,i_prof)) ! index of channel in t_tovs
          if (l_deb) write(usd,*) dpref,' transm_pred',j,ic,j_tt
          tr_loop: do it = 1, n_tr
            if (btest(s%bc%ib_ch(j_tt),it)) then
              select case(s%bc%type(it))
              case(1,2) ! P_THICK_TR, P_THICK_TRF
                if (l_deb) write(usd,*) dpref,' p_thick_tr pred',j,it,s%bc%p_tr(it,:)*0.01_wp
                kp(:) = -1
                do k = 1, n_lev
                  where (kp < 0 .and. r%p(k,i_p) >= s%bc%p_tr(it,:)*0.01_wp) kp(:) = k
                  if (all(kp > 0)) exit
                end do
                if (l_deb) write(usd,*) dpref,' p_thick_tr kp',kp
                do l = 1, 2
                  if (kp(l) < 0) then
                    l_fail(j) = .true.
                    exit tr_loop
                  elseif (kp(l) == 1) then
                    tr(l) = transm(kp(l),ic,i)
                  else
                    w = (log(s%bc%p_tr(it,l)*0.01_wp) - log(r%p(kp(l)-1,i_p))) / &
                        (log(r%p(kp(l),i_p))          - log(r%p(kp(l)-1,i_p)))
                    tr(l) = w * transm(kp(l),ic,i) + (1._wp - w) * transm(kp(l)-1,ic,i)
                  end if
                end do
                if (l_deb) write(usd,*) dpref,' p_thick_tr tr',tr
                select case(s%bc%type(it))
                case(1)
                  tt_tr(j_tt,it) = tr(2) - tr(1)
                case(2)
                  tt_tr(j_tt,it) = tr(2)
                end select
              case(3,4)
                ! Calc. transm thresholds
                if (s%bc%type(it) == 4) then
                  if (k_surf == 1) then
                    tr_surf = transm(k_surf,ic,i) + w_surf * &
                             (transm(k_surf,ic,i) - transm(k_surf+1,ic,i))
                  else
                    tr_surf = w_surf * transm(k_surf-1,ic,i) + (1._wp-w_surf) * transm(k_surf,ic,i)
                  end if
                  tr_surf = max(min(tr_surf,1._wp),0._wp)
                  if (l_deb) write(usd,*) dpref,' p_trans_ts tr_surf',tr_surf
                  tr = tr_surf + (1._wp - tr_surf) * s%bc%p_tr(it,:)*0.01_wp
                else
                  tr = s%bc%p_tr(it,:)*0.01_wp
                end if
                if (tr(2) < tr(1)) then ; t_=tr(2) ; tr(2)=tr(1) ; tr(1)=t_ ; endif
                if (l_deb) write(usd,*) dpref,' p_trans_ts tr',tr
                ! Calc. nearest levels below thresholds
                kp(:) = n_lev + 1
                do k = 1, n_lev
                  where (kp > n_lev .and. transm(k,ic,i) <= tr(:)) kp(:) = k
                  if (all(kp <= n_lev)) exit
                end do
                if (s%bc%type(it) == 4) then
                  kp(1) = min(kp(1),k_surf+1)
                else
                  if (any(kp(:) > n_lev)) then
                    l_fail(j) = .true.
                    exit tr_loop
                  end if
                end if
                if (l_deb) write(usd,*) dpref,' p_trans_ts kp',kp
                wsum = 0._wp
                ! Lowest Layer
                if (s%bc%type(it) == 4 .and. kp(1) > k_surf) then
                  ! Layer between surface and lowest level above surface (k_surf)
                  d_tr = (transm(k_surf,ic,i)-tr_surf)
                  t1   = T_surf
                  t2   = r%t_fg(k_surf,i_prof)
                elseif (kp(1) <= n_lev) then
                  ! Layer between low threshold and first level above low threshold (kp(1)-1)
                  d_tr = (transm(kp(1)-1,ic,i) - tr(1))
                  w    = (tr(1) - transm(kp(1),ic,i)) / (transm(kp(1)-1,ic,i) - transm(kp(1),ic,i))
                  t1   = w * r%t_fg(kp(1)-1,i_prof) + (1._wp-w) * r%t_fg(kp(1),i_prof)
                  t2   = r%t_fg(kp(1)-1,i_prof)
                end if
                t_   = (t1 + t2) * d_tr
                wsum = wsum + d_tr
                if (l_deb) write(usd,*) dpref,' p_trans_ts t_ s',t_,(t1+t2),d_tr,t_/wsum*0.5_wp,wsum
                ! Full layers between thresholds
                do l = kp(1)-1,kp(2)+1,-1
                  d_tr = (transm(l-1,ic,i)-transm(l,ic,i))
                  t1   = r%t_fg(l  ,i_prof)
                  t2   = r%t_fg(l-1,i_prof)
                  t_   = t_ + (t1 + t2) * d_tr
                  wsum = wsum + d_tr
                  if (l_deb) write(usd,*) dpref,' p_trans_ts t_',l,t_,(t1+t2),d_tr,t_/wsum*0.5_wp,wsum
                end do
                ! Layer between level below upper threshold and upper threshold
                if (kp(2) > 1) then
                  d_tr = (tr(2) - transm(kp(2),ic,i))
                  w  = (tr(2) - transm(kp(2),ic,i)) / (transm(kp(2)-1,ic,i) - transm(kp(2),ic,i))
                  t1 = w * r%t_fg(kp(2)-1,i_prof) + (1._wp-w) * r%t_fg(kp(2),i_prof)
                  t2 = r%t_fg(kp(2),i_prof)
                else
                  d_tr = (tr(2) - transm(kp(2),ic,i))
                  t1   = r%t_fg(kp(2),i_prof)
                  t2   = t1
                end if
                t_ = t_ + (t1 + t2) * d_tr
                wsum = wsum + d_tr
                if (l_deb) write(usd,*) dpref,' p_trans_ts t_',kp(2),t_,(t1+t2),d_tr,t_/wsum*0.5_wp,wsum
                if (abs(wsum) > 0._wp) then
                  tt_tr(j_tt,it) = t_ / wsum * 0.5_wp
                  if (l_deb) write(usd,*) dpref,' p_trans_ts tr',tt_tr(j_tt,it),t_,wsum,tr(2)-tr(1)
                else
                  call finish('calc_transm_pred','wsum == 0.')
                end if

              case(5,6)
                np = 1
                if (s%bc%type(it) == 6) np = 2
                ! Calc. nearest levels below thresholds
                kp(1:np) = n_lev + 1
                do k = 1, n_lev
                  where (kp(1:np) > n_lev .and. transm(k,ic,i) <= s%bc%p_tr(it,1:np)*0.01_wp) kp(1:np) = k
                  if (all(kp(1:np) <= n_lev)) exit
                end do
                if (any(kp(1:np) > n_lev)) then
                  l_fail(j) = .true.
                  exit tr_loop
                end if
                do k = 1, np
                  if (kp(k) == 1) then
                    w = (s%bc%p_tr(it,k)*0.01_wp - transm(kp(k),ic,i)) / (1._wp - transm(kp(k),ic,i))
                    p2 = 0._wp
                    p1 = r%p(kp(k),i_p)
                  else
                    w = (s%bc%p_tr(it,k)*0.01_wp - transm(kp(k),ic,i)) / (transm(kp(k)-1,ic,i) - transm(kp(k),ic,i))
                    p2 = r%p(kp(k)-1,i_p)
                    p1 = r%p(kp(k)  ,i_p)
                  end if
                  tr(k) = w * p2 + (1._wp - w) * p1
                end do
                if (np == 1) then
                  tt_tr(j_tt,it) = tr(1)
                else
                  tt_tr(j_tt,it) = tr(2) - tr(1)
                end if
              end select
              if (l_deb) write(usd,*) dpref,' transm_pred tr',j,it,tt_tr(j_tt,it)
            end if
          end do tr_loop
        end do

        if (any(l_fail)) then
          if (l_deb) write(usd,*) dpref,' transm_pred l_fail',count(l_fail(1:n_chan)),l_fail(1:n_chan)
          call decr_tovs_use(spt, obs%o(ib), CHK_BIASCOR, set=s, tovs=ttovs,&
             state=STAT_REJECTED, tovs_flag=TF_BC_FAILED, chan_mask=l_fail(1:n_chan))
        end if

        call store(obs% o(ib), spt, ttovs, tovs_io=TTOVS_TR, tr=tt_tr)
        call destruct(ttovs)

      end do

    end subroutine calc_transm_pred


  end subroutine process_tovs_mult

  subroutine K2H (H, temp_k, humi_k, t2m_k, q2m_k, psurf_k, stemp_k, &
                     ctp_k, cfr_k, emis_k, tj, pz_bg, n_chan, spot)
  !-----------------------------------------------------------------
  ! derive Jacobi matrix in 3D-Var interpolation space
  ! depending on tv (virtual temperature),
  !              gh (generalised humidity) and
  !              hs (geopotential height at surface pressure level)
  ! from quantities provided by the RTTOV K-routine.
  !----------------------------------------------------------------
  real(wp),    intent(out) :: H      (:,:) ! Jacobi matrix for 3D-Var
  real(wp),    intent(in)  :: temp_k (:,:) ! d Tb / d T(p)
  real(wp),    intent(in)  :: humi_k (:,:) ! d Tb / d q(p)
  real(wp),    intent(in)  :: t2m_k  (:)   ! d Tb / d T2m
  real(wp),    intent(in)  :: q2m_k  (:)   ! d Tb / d q2m
  real(wp),    intent(in)  :: psurf_k(:)   ! d Tb / d ps
  real(wp),    intent(in)  :: stemp_k(:)   ! d Tb / d ts
  real(wp),    intent(in)  :: ctp_k  (:)   ! d Tb / d ctp
  real(wp),    intent(in)  :: cfr_k  (:)   ! d Tb / d cfraction
  real(wp),    intent(in)  :: emis_k (:)   ! d Tb / d emis
  type(t_jac), intent(in)  :: tj           ! coefficients: Tv,gh<->T,q
  real(wp),    intent(in)  :: pz_bg        !               ps(hPa)/z(gpm)
  integer,     intent(in)  :: n_chan       ! number of channels
  type(t_spot),intent(in)  :: spot

    integer :: npv         ! number of profile variables
    integer :: k,j
    integer :: ideb

    ideb = 0
    if (ldeb(spot)) ideb = 2

    npv  = 2*size(temp_k,1)

!+++ TODO: handle cloud top & fraction

!NEC$ move
    do k = 1, n_chan
      H (k,1:npv  :2) = tj%dt_tv * temp_k(:,k)  & ! tv -> t
           + tj%dq_tv * humi_k(:,k)               ! tv -> q
      H (k,2:npv  :2) = tj%dt_rh * temp_k(:,k)  & ! rh -> t
           + tj%dq_rh * humi_k(:,k)               ! rh -> q
      H (k,2*tj%il-1) =  H (k,2*tj%il-1)        &
           + tj%wl * tj%dt_tv_2m * t2m_k(k)     &
           + tj%wl * tj%dq_tv_2m * q2m_k(k)
      H (k,2*tj%iu-1) =  H (k,2*tj%iu-1)        &
           + tj%wu * tj%dt_tv_2m * t2m_k(k)     &
           + tj%wu * tj%dq_tv_2m * q2m_k(k)
      H (k,2*tj%il)   =  H (k,2*tj%il)          &
           + tj%wl * tj%dt_rh_2m * t2m_k(k)     &
           + tj%wl * tj%dq_rh_2m * q2m_k(k)
      H (k,2*tj%iu)   =  H (k,2*tj%iu)          &
           + tj%wu * tj%dt_rh_2m * t2m_k(k)     &
           + tj%wu * tj%dq_rh_2m * q2m_k(k)

      H (k, npv+2)    =            psurf_k(k)     ! ps
      if (.not. tj%ps_llev) then
        H (k, npv+2)    =  H (k, npv+2)            &
           + ( t2m_k(k) * tj%dt_tv_2m * tj%d_tv    &
             + t2m_k(k) * tj%dt_rh_2m * tj%d_rh    &
             + q2m_k(k) * tj%dq_tv_2m * tj%d_tv    &
             + q2m_k(k) * tj%dq_rh_2m * tj%d_rh ) / tj%d_lev
      endif
      H (k, npv + 1 ) = tj% sigma_var_tskin * stemp_k(k)   ! skin temperature
      H (k, npv + 2 ) =     pz_bg           * H (k, npv+2) ! surface pressure
      H (k, npv + 3 ) = tj% dclt_dum        * ctp_k(k)     ! cloud top pressure
      H (k, npv + 4 ) = tj% dclf_dum        * cfr_k(k)     ! cloud fraction
      H (k, npv + 5 ) = tj% dsnf_dum        &              ! snow  fraction
                      * tj% de_dpc (0 ,k)   * emis_k(k)
      H (k, npv + 6:) = tj% de_dpc (1:,k)   * emis_k(k)    ! emissivity PC
    end do

    select case(ideb)
    ! case(1)
    !   write(*,*) 'K2H',po_context,spot%hd%id,'H(crc)',crc(reshape(H,(/spot%i%n*spot%o%n/)))
    case(2)
      write(usd,*) dpref,' K2H',po_context,spot%hd%id,'w',tj%wl,tj%wu
      write(usd,*) dpref,' K2H',po_context,spot%hd%id,'d2m',tj%dt_tv_2m,tj%dq_tv_2m,tj%dt_rh_2m,tj%dq_rh_2m
      write(usd,*) dpref,' K2H',po_context,spot%hd%id,'dc',tj%dclt_dum,tj%dclf_dum
      do k = 1, n_chan
        do j = 1, npv/2
          write(usd,*) dpref,' K2H',po_context,spot%hd%id,'t',k,j,H(k,(j-1)*2+1),tj%dt_tv(j),temp_k(j,k),tj%dq_tv(j),humi_k(j,k)
        end do
        do j = 1, npv/2
          write(usd,*) dpref,' K2H',po_context,spot%hd%id,'q',k,j,H(k,(j-1)*2+2),tj%dt_rh(j),temp_k(j,k),tj%dq_rh(j),humi_k(j,k)
        end do
        write(usd,*) dpref,' K2H',po_context,spot%hd%id,'2m',k,t2m_k(k),q2m_k(k),psurf_k(k),stemp_k(k),emis_k(k)
        write(usd,*) dpref,' K2H',po_context,spot%hd%id,'aux',k,ctp_k(k),cfr_k(k)
      end do
    end select

  end subroutine K2H
!==============================================================================
  subroutine radv_obs (r, obs, entry0)
  type (t_radv) ,intent(in)    ,target :: r
  type (t_obs)  ,intent(inout)         :: obs     ! observation data type
  integer       ,intent(in)            :: entry0
  !-----------------------------------------------
  ! fill radiances read into observation data type
  !-----------------------------------------------
    character(len=*), parameter :: proc = 'radv_obs'
    type(t_rad_set),  pointer   :: s     => NULL()
    type(t_rad_set)             :: rs
    type(t_head)                :: hd
    type(t_spot),     pointer   :: spot  => NULL()
    type(t_spot)                :: empty       ! default initialised
    type(t_spot)                :: spt
    type(t_use),      pointer   :: u     => NULL()
    type(t_use)                 :: use
    type(t_tovs),     target    :: ttovs
    character(len=132)          :: msg
    real(wp)                    :: surf_hgt, land_frac, instr_temp, cld_frc
    integer                     :: i_entry
    integer                     :: k_instr, k_sens, k_topo
    integer                     :: i_ch_cld_cov
    integer                     :: surf_flag
    integer                     :: surf_type
    integer                     :: flag
    integer                     :: i,j,k ! indices
    integer                     :: stat
    integer                     :: n_acc ! number of accepted profiles
    integer                     :: n_ch, n_act, n_pas, n_rej ! number of accepted profiles
    integer                     :: n_chan
    integer,        pointer     :: ii(:)
    integer,        allocatable :: chan_indx(:)
    logical,        allocatable :: mask(:)
    logical,        allocatable :: tmp_mask(:) ! Workaround for xlf 14.1.0.2
    logical,        allocatable :: l_vis(:)
    logical                     :: l_deb
    logical                     :: l_imager

    if (r% n_rec <= 0) return

    s => r% i

    allocate(mask     (s% n_chan), chan_indx(s% n_chan), l_vis(s%n_chan))
    if (associated(s% var)) deallocate(s% var)
    allocate(s% var(s%n_chan))

    ! Some variables in sat_pp depend on sensor or instrument in sat_pp,
    ! but are simple scalars in DACE. For these variables we have to make
    ! a decision which sensor/instrument shall be selected for DACE. The
    ! sensor might be selected directly with the select_sensor@tovs_obs_chan
    ! namelist entry. The flag_instr@tovs_obs_chan namelist entry is used to select
    ! an instrument (or a sensor for select_sensor <= 0).
    ! flag_instr meanings:
    ! >=0: get flag from instrument with this RTTOV ID
    !  -1: get flag from grid instrument
    !  -2: worst flag of all requested instruments
    !      (watch out for meaning of 'worst'!)
    !  -3: as for -1, but flag as 'mismatch'
    !      if inconsistent for diff. requested instruments


    ! Sensor to be selected for scalars that depend on the sensor
    k_sens = 0
    if (s% select_sens > 0) then
      if (s% select_sens > s% n_sens) then
        write(0,*) 'WARNING: select_sens=',s%select_sens,' < ',s%n_sens,'=n_sens !'
        write(0,*) 'WARNING: satid=',s%satid,' grid=',s%grid
        s%select_sens = s%n_sens
      end if
      k_sens = s%select_sens
    elseif (s% flag_instr >= 0) then
      do j = 1, s%n_sens
        if (s%sensor_instr(j) == s% flag_instr) then
          k_sens = j
          exit
        end if
      end do
    elseif ((s%flag_instr == -1) .or. (s%flag_instr == -3)) then
      do j = 1, s%n_sens
        if (s%sensor_instr(j) == s% grid) then
          k_sens = j
          exit
        end if
      end do
    elseif (s%flag_instr == -2) then
      k_sens = -1
    end if

    ! Instrument to be selected for scalars that depend on the instrument
    k_instr = 0
    if (s% flag_instr >= 0) then
      do j = 1, s%n_instr
        if (s% instr(j) == s% flag_instr) then
          k_instr = j
          exit
        end if
      end do
    elseif ((s%flag_instr == -1) .or. (s%flag_instr == -3)) then
      do j = 1, s%n_instr
        if (s%instr(j) == s% grid) then
          k_instr = j
          exit
        end if
      end do
    elseif (s%flag_instr == -2) then
      k_instr = -1
    end if
    if (k_instr == 0) call finish(proc, 'invalid flag_instr value')

    ! The topography variables might depend on sensor or instrument
    k_topo = 0
    if (size(r%landfr,1) == s%n_sens) then
      k_topo = k_sens
    else
      k_topo = k_instr
    end if
    if (k_topo == 0) call finish(proc, 'Failed to determine index for topography variables')

    ! colocated imager info
    l_imager = (s%n_im_chan > 0 .and. s%n_im_cluster > 0 .and. s%im_rttov_id >= 0)

    ! select channel for cloud_cover
    if (associated(r%cld_frc)) then
      i_ch_cld_cov = 0
      if (s%chan_cld_cov /= '') then
        select case(trim(tolower(s%chan_cld_cov)))
        case('avg')
          i_ch_cld_cov = -1
        case('max')
          i_ch_cld_cov = -2
        case('min')
          i_ch_cld_cov = -3
        case default
          read(s%chan_cld_cov, *, iostat=stat) k
          if (stat == 0 .and. k > 0 .and. k <= s%n_chan) then
            i_ch_cld_cov = k
          else
            write(0,*) '*** WARNING '//trim(proc)//': failed to interpret chan_cld_cov namelist entry "'//&
                 trim(s%chan_cld_cov)//'"'
            write(0,*) '(satid=',s%satid,' ,grid=',s%grid,')'
          end if
        end select
      end if
    end if
    
    ! Prepare header
    hd% modtype  = TOVS
    hd% obstype  = OT_RAD
    hd% buf_type = 21
    hd% source   = r% file_id
    hd% satid    = s% satid
    hd% grid_id  = s% grid

    n_acc  = 0
    n_ch   = 0
    n_act  = 0
    n_pas  = 0
    n_rej  = 0
    i_entry  = entry0

    rec_loop: do i = 1, r% n_rec
      i_entry           = i_entry + 1
      hd% time          = get_time(r% date(i), r% time(i))
      hd% subset        = r% fov (i)
      if (associated(r% obsnum   )) then
        select case(hd_id_vers)
        case(:1)
          hd% id        = r% obsnum(i)
        case(2:)
          hd% id        = r% obsnum(i) * 1000 + r% file_id
        end select
      else
        hd% id          = i_entry
      end if
      if (associated(r% center   )) hd% center        = r% center   (i)
      if (associated(r% subcenter)) then
        if (abs(r% subcenter(i)) <= huge(hd%subcenter)) then
          hd% subcenter     = r% subcenter(i)
        else
          hd% subcenter = -1
        end if
      end if
      if (associated(r% date_d).and.associated(r% time_d)) &
           hd% db_time = get_time(r% date_d(i), r% time_d(i))

      call check_report_0 (use, hd, 1)
      if (use% state <= STAT_DISMISS) cycle

      ! Topography variables
      surf_flag  = 0
      surf_type  = 0
      surf_hgt   = 0._wp
      land_frac  = 0._wp
      if (k_topo > 0) then
        surf_type = r% stype (k_topo, i)
        surf_hgt  = r% shgt  (k_topo, i)
        land_frac = r% landfr(k_topo, i)
      else
        surf_type = maxval(r% stype (:, i))
        surf_hgt  = maxval(r% shgt  (:, i))
        land_frac = maxval(r% landfr(:, i))
      end if
      if (s%flag_instr == -3) then
        if (minval(r% stype(:, i)) /= maxval(r% stype(:, i))) &
             surf_flag = ibset(surf_flag, SUR_MISMATCH)
      end if

      ! Optional variables
      cld_frc    = 0._wp
      if (associated(r%cld_frc)) then
        if (s%n_instr == size(r%cld_frc,1)) then
          if (k_instr > 0) then
            cld_frc = r% cld_frc(k_instr,i)
          else
            cld_frc = maxval(r% cld_frc(:,i))
          end if
        elseif (s%n_chan == size(r%cld_frc,1)) then
          select case(i_ch_cld_cov)
          case(-3)
            cld_frc = minval(r%cld_frc(:,i), mask=r%cld_frc(:,i)/=NF90_FILL_INT)
          case(-2)
            cld_frc = maxval(r%cld_frc(:,i), mask=r%cld_frc(:,i)/=NF90_FILL_INT)
          case(-1)
            k = count(r%cld_frc(:,i)/=NF90_FILL_INT)
            if (k > 0) then
              cld_frc = sum(r%cld_frc(:,i), mask=r%cld_frc(:,i)/=NF90_FILL_INT)/real(k)
            else
              cld_frc = -999._wp
            end if
          case(0)
            ! We could do something useful here, but for backwards compatibility reasons
            ! we set to invalid here (relevant e.g. in mo_tskin)
            cld_frc = -999._wp
          case(1:)
            if (i_ch_cld_cov > s%n_chan) call finish(proc, 'invalid i_ch_cld_cov (>s%n_chan)')
            cld_frc = r%cld_frc(i_ch_cld_cov,i)
          case default
            call finish(proc, 'invalid i_ch_cld_cov (<-3)')
          end select
        end if
      end if
      instr_temp = 0._wp
      if (associated(r%instr_temp)) then
        if (k_sens > 0) then
          instr_temp = r% instr_temp(k_sens, i)
        else
          instr_temp = maxval(r% instr_temp(:, i))
        end if
      end if

      ! surface flag
      select case (surf_type)
      case (0) ! sea
        surf_flag = ibset(surf_flag, SUR_SEA)
      case (1) !coast
        surf_flag = ibset(surf_flag, SUR_LAND)
        surf_flag = ibset(surf_flag, SUR_SEA)
      case (2) ! land
        surf_flag = ibset(surf_flag, SUR_LAND)
      case default
        surf_flag = ibset(surf_flag, SUR_MISSING)     ! catch missing values
      end select
      if (surf_hgt > highland_thresh) then
        surf_flag = ibset(surf_flag, SUR_HIGHLAND)
      endif

      if (s% gopts% pp_flags_instr > 0 .and. s% gopts% pp_flags > 0) then
        j = s% gopts% pp_flags_instr
        if (j > size(r% pp_flags,1)) then
!           print*,'invalid:',j,size(r% pp_flags,1)
!           call finish('radv_obs', 'invalid pp_flags_instr')
          write(0,'(A,2(2x,I3))') '*** WARNING: invalid pp_flags_instr',j,size(r% pp_flags,1)
        else
          ! sat_pp does not know weighting functions, and is not able to do a weighting
          ! function dependent blacklisting. For this reason the evaluation of the PP_FLAGS
          ! was implmented, such that the information from sat_pp might be used to blacklist
          ! low peaking observations, while high peaking observations might be kept. Since
          ! the sat_pp physical checks usually detect "low" tropospheric deterioration, we
          ! use the 3dvar surface flag here.
          if (iand(r% pp_flags(j,i), s% gopts% pp_flags) /= 0) then
            surf_flag = ibset(surf_flag, SUR_BLK_PP)
          end if
        end if
      end if


      ! determine valid channels
      where (btest(s% flag(1:s%n_chan), USE_BCOR)  .or. &
             ((r% bt_obs(1:s%n_chan, i) < 500._sp) .and.&
              (r% bt_obs(1:s%n_chan, i) >   0._sp)))
        mask(:) = .true.
      elsewhere
        mask(:) = .false.
      end where
      n_chan = count(mask(:))
      chan_indx(1:n_chan) = pack( (/ (i, i=1, s%n_chan) /), mask(:))
      if (n_chan == 0) cycle

      spt               = empty
      spt% use          = use
      spt% hd           = hd
      spt% int_type     = ITY_ICOL
      spt% sttyp        = s% grid_wmo
      spt% ident        = s% satid
      spt% statid       = satname(s% satid)
      spt% actual_time  = hd% time
      spt% col% c% dlat = r% dlat(i)
      spt% col% c% dlon = r% dlon(i)
      spt% col% nlev    = n_chan
      spt% stzen        = r% stzen(i)
      spt% stazi        = r% stazi(i)
      if (r%fov(i) == -huge(r%fov(i))) then
        spt%phase       = -huge(spt%phase)
      else
        spt%phase       = r%fov(i)
      end if
      spt% char         = CHR_NONL+CHR_EXP
      spt% cost         = n_chan * 5
      spt% stlsf        = surf_flag
      spt% stclf        = 5            ! IR clear and MW clear
      spt% pcc          = 100
      if (associated(r% sunzen)) spt% sozen = r% sunzen(i)
      if (associated(r% sunazi)) spt% soazi = r% sunazi(i)

      l_deb = ldeb(spt)
      if (l_deb) write(usd,*) dpref,'radv_obs',s%satid,s%n_instr,s%instr(1:s%n_instr),n_chan,&
           &r% dlat(i),r% dlon(i),hd%time

      call check_report_1 (spt)
      if (iand(flg_sur, surf_flag) /= 0)  call decr_rpt_use(spt, CHK_SURF)
      if (spt% stzen > 75._wp)            call decr_rpt_use(spt, CHK_CONSIST, &
           use=STAT_DISMISS, comment='sat. zenith angle.')
      if (spt% use% state <= STAT_DISMISS) cycle

      n_acc = n_acc + 1
      call new_spot(obs, 1, set_id=.true.)
      spot     => obs% spot (obs% n_spot)
      j        =  spot% id
      spot     =  spt
      spot% id =  j
      spot% nr =  n_chan
      call set_xuv(spot)

      rs = reduced_rad_set(s, mask(:), oerr_par=.false.)

      call new_obs(obs, n_chan, spot=spot)
      obs% olev (spot%o%i+1:spot%o%i+n_chan)       = rs% chan(1:n_chan)
      obs% body (spot%o%i+1:spot%o%i+n_chan)%o     = r% bt_obs(chan_indx(1:n_chan),i)
      obs% body (spot%o%i+1:spot%o%i+n_chan)%bc    = 0._wp
      obs% body (spot%o%i+1:spot%o%i+n_chan)%op_na = 0
      obs% varno(spot%o%i+1:spot%o%i+n_chan)       = VN_RAWBT

      l_vis = .false.
      if (associated(rs%iflag)) then
        l_vis(1:n_chan) = (rs%iflag(1:n_chan) == 1)
        where(l_vis(1:n_chan)) obs% varno(spot%o%i+1:spot%o%i+n_chan) = VN_REFL
      end if

!!! Original code
!     where (pack(r% valid(1:s%n_chan, i), mask(:)))
!       obs% body (spot%o%i+1:spot%o%i+n_chan)%ch_bl = 0
!     elsewhere
!       obs% body (spot%o%i+1:spot%o%i+n_chan)%ch_bl = 1
!     end where
!!! Alternative code
!     obs% body (spot%o%i+1:spot%o%i+n_chan)%ch_bl = &
!          merge (0, 1, pack(r% valid(1:s%n_chan, i), mask(:)))
!!! Workaround for xlf 14.1.0.2 (won't finish compilation)
      allocate (tmp_mask (n_chan))
      tmp_mask(:) = pack(r% valid(1:s%n_chan, i), mask(:))
      where (tmp_mask(:))
        obs% body (spot%o%i+1:spot%o%i+n_chan)%ch_bl = 0
      elsewhere
        obs% body (spot%o%i+1:spot%o%i+n_chan)%ch_bl = 1
      end where
      deallocate (tmp_mask)

      do j = 1, rs%n_instr
        obs% body (spot%o%i+rs%o_ch_i(j)+1:&
             spot%o%i+rs%o_ch_i(j)+rs%n_ch_i(j))% lev_sig = rs% instr_wmo(j)
      end do

      ! Initialize t_tovs for this spot
      call construct(ttovs)
      ttovs% nchan = n_chan
      ttovs% nlev  = jplev
      ttovs% nl2c  = sum(rs%n_ch_i(1:rs%n_instr), &
                         mask=rs%iopts(1:rs%n_instr)%l2c_type > 0)
      ttovs% nspec = count(rs%iopts(1:rs%n_instr)%do_lambertian)
      if (any(btest(rs%iopts(1:rs%n_instr)%rad_out, OUT_CSB))) then
        ttovs% nbtcs = n_chan
      else
        ttovs% nbtcs = 0
      end if
      allocate(ttovs% ci  (ttovs% nchan))
      ttovs% ci         = chan_indx(1:n_chan)
      allocate(ttovs% emis(ttovs% nchan))
      ttovs% emis(:)    = -1._wp
      allocate(ttovs% flag(ttovs% nchan))
      ttovs% flag(:)    = 0
      ttovs% init       = TTOVS_BASE + TTOVS_CI + TTOVS_EMIS + TTOVS_FLAG
      if (ttovs% nl2c > 0) then
        allocate(ttovs% l2c (ttovs% nl2c ))
        ttovs% l2c (:)  = -1._sp
        ttovs% init     = ttovs% init + TTOVS_L2C
      end if
      if (ttovs% nspec > 0) then
        allocate(ttovs% spec (ttovs% nspec))
        ttovs% spec(:)  = 1._sp
        ttovs% init     = ttovs% init + TTOVS_SPEC
      end if
      if(ttovs % nbtcs > 0 ) then
        allocate(ttovs% bt_cs(ttovs%nchan))
        ttovs% bt_cs(:) = -1._sp
        ttovs% init     = ttovs% init + TTOVS_BTCS
      end if
      if (any(rs%iopts(1:rs%n_instr)%l_cldlev)) then
        ttovs% nband = 0
        do j = 1, rs%n_instr
          if (rs%iopts(j)%l_cldlev) then
            ttovs% nband = max(ttovs%nband, maxval(rs%band(rs%o_ch_i(j)+1:rs%o_ch_i(j)+rs%n_ch_i(j))))
          end if
        end do
        if (ttovs%nband > m_bd) call finish(proc, 'maximum number of bands exceeded')
        if (ttovs%nband > 0) then
          allocate(ttovs%cldlev(ttovs%nband))
          ttovs% cldlev(:) = -1._sp
          ttovs% init      = ttovs% init + TTOVS_CLDLEV
        end if
      end if
      ttovs% id_rs      = s%id
      ttovs% saza       = real(r% stzen(i), kind=sp)
      ttovs% boa        = real(r% stazi(i), kind=sp)
      if (associated(r% scanl     )) ttovs% scanl      =      r% scanl(i)
      if (associated(r% sunzen    )) ttovs% soza       = real(r% sunzen(i), kind=sp)
      if (associated(r% sunazi    )) ttovs% soa        = real(r% sunazi(i), kind=sp)
      if (associated(r% orb_phase )) ttovs% orb_phase  = real(r% orb_phase(i), kind=sp)
      if (associated(r% nwc_flg   )) ttovs% nwc_flg    =      r% nwc_flg(i)
      if (associated(r% instr_temp)) ttovs% instr_temp = real(instr_temp, kind=sp)
      if (associated(r% cld_frc   )) ttovs% cloud_imag = real(cld_frc, kind=sp)
      if (associated(r% sunzen) .and. associated(r% sunazi)) &
           ttovs% sasoa = real(sat_sun_azimuth_angle(r%stzen(i), r%stazi(i), &
                               r%sunzen(i), r%sunazi(i)), kind=sp)
      ttovs% land_frac  = real(land_frac, kind=sp)
      if (l_imager) then
        ttovs%n_im_cl   = s%n_im_cluster
        ttovs%n_im_cl_v = 2 * s%n_im_chan + 1
        ttovs%n_im_ch   = s%n_im_chan
        ttovs%n_im_ch_v = mx_imch
        allocate(ttovs% im_cl(ttovs%n_im_cl_v,ttovs%n_im_cl),ttovs% im_ch(ttovs%n_im_ch_v,ttovs%n_im_ch))
        ttovs% im_cl = 0._sp
        ttovs% im_ch = 0._sp
        k = 1
        ttovs% im_cl_v(k  ) = IMCL_FRAC
        ttovs% im_cl  (k,:) = r% im_cl_frac(1:s%n_im_cluster,i) * 0.01_sp
        if (l_deb) then
          write(usd,*) dpref,'radv_obs imager s%im_chan',s%im_chan(s%n_im_chan)
          write(usd,*) dpref,'radv_obs imager r%im_chan',r%im_chan(:)
          write(usd,*) dpref,'radv_obs imager frac',ttovs% im_cl  (k,:)
        end if
        do j = 1, s%n_im_chan
          k = k + 1
          ttovs% im_cl_v(k  ) = 100 * s%im_chan(j) + IMCL_MEAN
          ttovs% im_cl  (k,:) = real(pack(r% im_cl_mean(:,i), mask=(r%im_chan(:)==s%im_chan(j))), kind=sp)
          k = k + 1
          ttovs% im_cl_v(k  ) = 100 * s%im_chan(j) + IMCL_STDV
          ttovs% im_cl  (k,:) = real(pack(r% im_cl_stdv(:,i), mask=(r%im_chan(:)==s%im_chan(j))), kind=sp)
          if (l_deb) then
            write(usd,*) dpref,'radv_obs imager count(mask)',count(r%im_chan(:)==s%im_chan(j))
            write(usd,*) dpref,'radv_obs imager mean',j,k,ttovs% im_cl  (k-1,:)
            write(usd,*) dpref,'radv_obs imager stdv',j,k,ttovs% im_cl  (k  ,:)
          end if
        end do
        ttovs% init = ttovs% init + TTOVS_IMCL + TTOVS_IMCH
      end if

      ! This must be called before the decr_tovs_use call with tovs_flag option below!
      call store (obs, spot, ttovs)

      if (use_ch_blcklst) then
        call decr_tovs_use(spot, obs, CHK_DATASET, set=rs, tovs=ttovs, state=STAT_PAS_REJ, &
             chan_mask=.not.r%valid(chan_indx(1:n_chan),i), flag_rpt=.false., hint_debug='ch_blcklst')
      end if
      call decr_tovs_use(spot, obs, CHK_NOTUSED, set=rs, tovs=ttovs, state=STAT_PASSIVE, &
           bit=USE_PASSIVE, flag_rpt=.false., hint_debug='passive')


      !LB:: QC for obs with VN_REFL
      if (any(l_vis(1:n_chan))) then
        !OBS > 1.5
        if (any(l_vis(1:n_chan) .and. obs% body (spot%o%i+1:spot%o%i+n_chan)%o >= 1.5_wp)) then
          call decr_tovs_use(spot, obs, CHK_DATASET, set=rs, tovs=ttovs, state=STAT_REJECTED,&
               chan_mask=(l_vis(1:n_chan) .and. obs% body (spot%o%i+1:spot%o%i+n_chan)%o >= 1.5_wp))
        end if
        !SZA > 75
        if (ttovs%soza .ge. 75._wp) then
          call decr_tovs_use(spot, obs, CHK_OPERATOR, set=rs, tovs=ttovs, state=STAT_REJECTED, &
               chan_mask=l_vis(1:n_chan))
        end if
        !Snow
        if ( btest(ttovs%nwc_flg,2)) then
          call decr_tovs_use(spot, obs, CHK_CLOUD, set=rs, tovs=ttovs, state=STAT_REJECTED, &
               tovs_flag=TF_SURF_INFL, chan_mask=l_vis(1:n_chan))
        end if
      end if
      ! Dust and volcanic ash QC also for WV channels
      !Dust
      if ( btest(ttovs%nwc_flg,5)) then
        call decr_tovs_use(spot, obs, CHK_CLOUD, set=rs, tovs=ttovs, state=STAT_REJECTED, &
             tovs_flag=TF_DESERT_DUST )
      end if
      !Volcanic ash
      if ( btest(ttovs%nwc_flg,8)) then
        call decr_tovs_use(spot, obs, CHK_CLOUD, set=rs, tovs=ttovs, state=STAT_REJECTED,&
             tovs_flag=TF_VOLCANIC_ASH)
      end if

      call destruct(ttovs)

      do j = 1, n_chan
        u    => obs% body(spot%o%i+j)% use
        flag =  rs% flag(j)

        if (btest(flag, USE_BCOR) .and. u% state < STAT_PAS_REJ) &
             call incr_use(u,state=STAT_PAS_REJ,check=CHK_BIASCOR)
        n_ch = n_ch + 1
        if (u% state >= STAT_ACTIVE) then
          n_act = n_act + 1
        elseif (u% state >= STAT_PASSIVE) then
          n_pas = n_pas + 1
        else
          n_rej = n_rej + 1
        end if
      end do

    end do rec_loop

    write(msg,'(6x,"store : ",I7," FOVs accepted, ",I8,&
         &" ch., ",I8," act., ",I8," pass., ",I8," rej. (im=",L1,")")') &
         n_acc, n_ch, n_act, n_pas, n_rej, l_imager
    call add_line(msg)

  contains

    type(t_time) function get_time(date, time)
      integer, intent(in) :: date
      integer, intent(in) :: time

      integer              :: yyyy, mo, dd
      integer              :: hh, mi, sec

      yyyy =      date / 10000
      mo   = mod (date/   100 ,100)
      dd   = mod (date        ,100)
      hh   =      time/ 10000
      mi   = mod (time/   100 ,100)
      sec  = mod (time        ,100)
      call init_time (get_time, yyyy, mo, dd, hh, mi, sec)

    end function get_time

  end subroutine radv_obs


  subroutine read_fdbk_rad(task, obs, i_s, i_e, nh, nb)
    character(len=*), intent(in)             :: task
    type (t_obs),     intent(inout), target  :: obs(:)          ! observation data type variable
    integer,          intent(in),   optional :: i_s
    integer,          intent(in),   optional :: i_e
    integer,          intent(in),   optional :: nh
    integer,          intent(in),   optional :: nb

    character(len=*), parameter   :: proc = 'read_fdbk_rad'
    character(len=120)            :: msg  = ''
    integer,          save        :: n_fdbk  = 0
    type (t_obs),     pointer     :: o
    type(t_spot),     pointer     :: spt
    type (t_tovs)                 :: tovs       ! RAD obstype specific data
    type (t_rad_set), pointer     :: rs      => null()
    real(sp),         allocatable :: l2c(:), emiss(:), ins_tmp(:), orb_ph(:), cld_frc(:), sinfl(:)
    integer,          allocatable :: nwc_flg(:), tovs_flag(:), fov(:), scanl(:)
    integer                       :: i1, in, n
    integer                       :: i, j, k, l, m, ib, is, ic, iset, ii, ih
    integer                       :: i1_b, in_b, nbdy
    integer                       :: satid, grid, chan, instr
    integer                       :: ierr_emis, ierr_tf, ierr_fov, ierr_sl, ierr_si
    integer                       :: ierr(0:n_optv-1)
    integer                       :: ci(m_chan)

    ! task='post', i.e. unify rad_set:
    type t_rad_set_aux
      integer :: satid
      integer :: grid
      integer :: n_instr
      integer :: n_chan
      integer :: opt_vars
      integer, pointer :: instr (:) => null()
      integer, pointer :: n_ch_i(:) => null()
      integer, pointer :: o_ch_i(:) => null()
      integer, pointer :: chan  (:) => null()
      integer, pointer :: iflag (:) => null()
    end type t_rad_set_aux
    type(t_rad_set_aux), allocatable,target :: rsa(:)
    type(t_rad_set_aux), allocatable,target :: rsa_u(:)
    type(t_rad_set_aux), pointer     :: rsa_, rsa_o
    logical,             allocatable :: mask_rs(:)
    logical,             allocatable :: mask_iflag(:)
    logical                          :: l_iflag
    integer,             allocatable :: isend(:), isend2(:), irecv(:)
    integer                          :: mx_instr, mx_chan
    integer                          :: nrecv, nsend
    integer                          :: n_rs, n_rs_all
    integer                          :: irs
    integer                          :: n_iflag

    ! ! abort if SATPP files are read at the same time
    ! ! TODO: is this really necessary?
    ! if (n_filetype (FT_SATPP) > 0) &
    !   call finish ('read_fdbk_rad','n_filetype (FT_SATPP) > 0')

    select case(task)
    case('read')
      n_fdbk = n_fdbk + 1

      if (.not.present(nh).or..not.present(nb).or..not.present(i_s).or..not.present(i_e)) &
           call finish(proc,'nh/nb/i_s/i_e arguments required for task=read')

      allocate(l2c(nb),emiss(nb),tovs_flag(nb),sinfl(nb),fov(nh),scanl(nh),&
           nwc_flg(nh),orb_ph(nh),ins_tmp(nh),cld_frc(nh))

      stanc = 1
      counc = nb
      strnc = 1
      call get_var(emiss,     'emissiv',   ierr=ierr_emis     )
      call get_var(tovs_flag, 'tovs_flag', ierr=ierr_tf       )
      call get_var(sinfl,     'surf_infl', ierr=ierr_si       )
      call get_var(l2c,       'l2c',       ierr=ierr(OPTV_L2C))
      stanc = 1
      counc = nh
      strnc = 1
      call get_var(fov,     'fov',         ierr=ierr_fov          )
      call get_var(scanl,   'scanline',    ierr=ierr_sl           )
      call get_var(nwc_flg, 'nwc_flag',    ierr=ierr(OPTV_NWC_FLG))
      call get_var(orb_ph , 'orbit_phase', ierr=ierr(OPTV_ORB_PH ))
      call get_var(ins_tmp, 'instr_temp',  ierr=ierr(OPTV_INS_TMP))
      call get_var(cld_frc, 'cloud_frac',  ierr=ierr(OPTV_CLD_FRC))
      ! Set up tovs specific stuff, i.e. t_tovs and rad_set
      nbdy = 0
      do ib = 1, size(obs)
        o => obs(ib)
        spot_loop: do is = i_s,i_e
          spt => o% spot(is)
          if (spt% hd% obstype /= OT_RAD) cycle
          ! if (spt% hd% mon_file /= ifile) cycle
          ih = spt% hd% mon_rec
          n     = spt% o% n
          i1    = spt% o% i + 1
          in    = spt% o% i + n
          i1_b  = nbdy + 1
          in_b  = nbdy + n
          satid = spt% hd% satid
          grid  = rttov_instr (spt% sttyp, satid)
          if (grid < 0) then
            write(msg,'(I3.3)') spt%sttyp
            call finish(proc,'unknown grid (instype) = '//trim(msg))
          end if
          spt% hd% grid_id = grid
          spt% hd% id      = spt%hd%mon_rec * 1000 + spt%hd%mon_file
          ! Setup correspoding rad_set
          iset = set_indx(rad_set(1:n_set), satid=satid, grid=grid)
          if (iset > 0 .and. iset <= n_set) then
            rs => rad_set(iset)
            if (rs%id < 0) rs%id = iset  ! Indicated that this rad_set connected to data
          else
            n_set = n_set + 1
            if (n_set > m_rad_set) call finish(proc,'n_set > m_rad_set')
            rs => rad_set(n_set)
            rs%source   = 'fdbk'
            rs%id       = n_set
            rs%satid    = satid
            rs%grid     = grid
            rs%grid_wmo = instr_rttov(grid, satid)
            allocate(rs%chan(m_chan))
            call satid_bufr2rttov(rs% satid, rs% rttov_satid, platform=rs% platform)
          end if

          ! Loop over channels and add to rad_set if necessary
          do j = i1, in
            i = j-i1+1
            chan  = nint(o%olev(j))
            instr = rttov_instr(int(o%body(j)%lev_sig), satid)
            if (instr < 0) then
              write(msg,'(I3.3)') int(o%body(j)%lev_sig)
              call finish(proc,'unknown instrument (lev_sig) = '//trim(msg))
            end if
            k = chan_indx(instr, chan, rs, ii=ii)
            if (rs%source == 'nml') then
              if (k <= 0) then
                write(0,*) 'is=',is,' spt%hd%id',spt%hd%id,' j=',j
                write(0,*) 'satid=',satid,' grid=',grid,' instr=',instr,' chan=',chan
                call print_rad_set(rs, hint=dace%pe, header='crash', unit=stderr)
                call finish(proc, 'Failed to to find channel in rad_set')
              else
                ci(i) = k
              end if
            else
              ci(i) = -1
              if (k <= 0) then
                if (ii <= 0) then
                  ! Add instrument to rad_set
                  rs%n_instr = rs%n_instr + 1
                  if (rs%n_instr > m_instr) call finish(proc, 'n_instr > m_instr')
                  ii = rs%n_instr
                  rs%instr    (ii) = instr
                  rs%instr_wmo(ii) = instr_rttov(instr,satid)
                  if (ii > 1) then
                    rs%o_ch_i(ii) = rs%o_ch_i(ii-1) + rs%n_ch_i(ii-1)
                  else
                    rs%o_ch_i(ii) = 0
                  end if
                  rs%n_ch_i(ii) = 0
                end if
                ! Add channel to rad_set
                ic = rs%o_ch_i(ii) + rs%n_ch_i(ii) + 1
                do l = rs%o_ch_i(ii)+1, rs%o_ch_i(ii)+rs%n_ch_i(ii)
                  if (rs%chan(l) > chan) then
                    ic = l
                    exit
                  end if
                end do
                k = ic
                rs%chan (ic+1:rs%n_chan+1) = rs%chan(ic:rs%n_chan)
                rs%chan (ic)               = chan
                rs%n_chan = rs%n_chan + 1
                rs%n_ch_i(ii) = rs%n_ch_i(ii) + 1
                rs%o_ch_i(ii+1:rs%n_instr) = rs%o_ch_i(ii+1:rs%n_instr) + 1
                if (associated(rs%iflag)) then
                  rs%iflag(ic+1:rs%n_chan+1) = rs%iflag(ic:rs%n_chan)
                  rs%iflag(ic)               = 0
                end if
              end if
            end if
            if (o%varno(j) == VN_REFL) then
              if (.not.associated(rs%iflag)) then
                allocate(rs%iflag(1:m_chan))
                rs% iflag = 0
              end if
              rs% iflag(k) = 1
            end if
          end do

          ! Set up t_tovs
          call construct (tovs)
          if (rs%source == 'nml') then
            tovs%id_rs = rs%id
          else
            tovs%id_rs = -1
          end if
          tovs% nchan = n
          tovs% nlev  = jplev
          tovs% init  = TTOVS_BASE + TTOVS_CI + TTOVS_EMIS + TTOVS_FLAG
          allocate (tovs% ci(n), tovs% emis(n), tovs% flag(n))
          tovs% ci    = ci(1:n)
          tovs% nspec = count(rs%iopts(1:rs%n_instr)%do_lambertian)
          if (any(btest(rs%iopts(1:rs%n_instr)%rad_out, OUT_CSB))) then
            tovs%init = tovs%init + TTOVS_BTCS
            allocate(tovs% bt_cs(n)) ; tovs% bt_cs = -1._sp
          end if
          if (ierr_emis == 0) then
            tovs%emis(1:n) = emiss(i1_b:in_b)
          else
            tovs%emis(1:n) = -1._wp
          end if
          if (ierr_tf == 0) then
            tovs%flag(1:n) = tovs_flag(i1_b:in_b)
          else
            where (tovs%emis(1:n) > 0._wp .and. tovs%emis(1:n) < 1._wp)
              tovs%flag(1:n) = ibset(0,TF_EMIS_FASTEM)
            elsewhere
              tovs%flag(1:n) = 0
            end where
          end if
          if (ierr_si == 0) then
            allocate(tovs%sinfl(n))
            tovs%init = tovs%init + TTOVS_SINFL
            tovs%sinfl(1:n) = sinfl(i1_b:in_b)
          end if
          if (ierr_sl == 0) then
            tovs%scanl = scanl(ih)
          end if
          if (ierr_fov == 0) then
            spt%phase = fov(ih)
          end if
          ! Optional variables
          if (ierr(OPTV_L2C) == 0) then
            rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_L2C)
            allocate (tovs% l2c(n))
            tovs% init  = tovs% init + TTOVS_L2C
            ii = 1
            i  = 0
            do j = 1, n
              instr = o% body(spt%o%i+j)% lev_sig
              do while (rs%instr_wmo(ii) /= instr)
                ii = ii + 1
              end do
              if (rs% iopts(ii)% l2c_type > 0) then
                i = i + 1
                tovs%l2c(i) = l2c(nbdy+j)
              end if
            end do
            tovs% nl2c = i
          end if
          if (ierr(OPTV_NWC_FLG) == 0) then
            rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_NWC_FLG)
            tovs%nwc_flg = nwc_flg(ih)
          else
            tovs%nwc_flg = 0
          end if
          if (ierr(OPTV_ORB_PH ) == 0) then
            rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_ORB_PH )
            tovs%orb_phase = orb_ph(ih)
          else
            tovs%orb_phase = 0._sp
          end if
          if (ierr(OPTV_INS_TMP) == 0) then
            rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_INS_TMP)
            tovs%instr_temp = ins_tmp(ih)
          else
            tovs%instr_temp = 0._sp
          end if
          if (ierr(OPTV_CLD_FRC) == 0) then
            rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_CLD_FRC)
            tovs%cloud_imag = cld_frc(ih)
          else
            tovs%cloud_imag = 0._sp
          end if
          if (btest(spt%stlsf, SUR_SEA) .and. btest(spt%stlsf, SUR_LAND)) then
            tovs%land_frac = 0.5_sp
          elseif (btest(spt%stlsf, SUR_SEA)) then
            tovs%land_frac = 0.0_sp
          elseif (btest(spt%stlsf, SUR_LAND)) then
            tovs%land_frac = 1.0_sp
          else ! invalid, should not happen
            tovs%land_frac = 1.0_sp
          end if
 
          tovs% saza       = real (spt% stzen,sp)       ! sat.zenith angle
          tovs% boa        = real (spt% stazi,sp)       ! satellite azimuth
          tovs% soza       = real (spt% sozen,sp)       ! solar zenith angle
          tovs% soa        = real (spt% soazi,sp)       ! solar azimuth angle
          ! set obsolete or unused entries (hopefully)
          tovs% soa        = 0._sp                      ! solar     azimuth
          tovs% sasoa      = 0._sp                      ! sat-sol.  azimuth

          call store (o, spt, tovs)

          if (associated(rs%flag) .and. rs%source == 'nml') then
            call decr_tovs_use(spt, o, CHK_NOTUSED, set=rs, tovs=tovs, state=STAT_PASSIVE, &
                 bit=USE_PASSIVE, flag_rpt=.false., hint_debug='passive')
          endif

          call destruct (tovs)

          nbdy = nbdy + n

        end do spot_loop
      end do
    case('post')
      ! ----------------------------------
      ! --- Unify rad_set on all procs ---
      ! ----------------------------------
      ! --- rad_set created by read_fdbk_rad ---
      n_rs = p_max(n_set)
      if (n_rs <= 0) RETURN

      allocate(mask_rs(n_set))
      mask_rs(1:n_set) = (rad_set(1:n_set)%source == 'fdbk')
      n_rs = count(mask_rs(1:n_set))
      n_rs_all = p_sum(n_rs)
      if (n_rs_all > 0) then
        ! --- Broadcast rad_set into rad_set_aux on dace%pio ---
        if (dace%lpio) then
          nrecv = n_rs_all
        else
          nrecv = 0
        end if
        allocate(rsa(nrecv))
        ! Basic info
        call p_gatherv(pack(rad_set(1:n_set)%satid   , mask=mask_rs), rsa(:)%satid   , dace%pio)
        call p_gatherv(pack(rad_set(1:n_set)%grid    , mask=mask_rs), rsa(:)%grid    , dace%pio)
        call p_gatherv(pack(rad_set(1:n_set)%n_instr , mask=mask_rs), rsa(:)%n_instr , dace%pio)
        call p_gatherv(pack(rad_set(1:n_set)%n_chan  , mask=mask_rs), rsa(:)%n_chan  , dace%pio)
        call p_gatherv(pack(rad_set(1:n_set)%gopts%opt_vars, mask=mask_rs), rsa(:)%opt_vars, dace%pio)
        ! Instrument info
        mx_instr = maxval(rad_set(1:n_set)%n_instr, mask=mask_rs)
        mx_instr = p_max(mx_instr)
        do irs = 1, nrecv
          allocate(rsa(irs)%instr(mx_instr), rsa(irs)%n_ch_i(mx_instr), rsa(irs)%o_ch_i(mx_instr))
        end do
        allocate(irecv(nrecv))
        do i = 1, mx_instr
          call p_gatherv(pack(rad_set(1:n_set)%instr(i), mask=mask_rs), irecv, dace%pio)
          do irs = 1, nrecv ; rsa(irs)%instr(i) = irecv(irs) ; enddo
          call p_gatherv(pack(rad_set(1:n_set)%n_ch_i(i), mask=mask_rs), irecv, dace%pio)
          do irs = 1, nrecv ; rsa(irs)%n_ch_i(i) = irecv(irs) ; enddo
          call p_gatherv(pack(rad_set(1:n_set)%o_ch_i(i), mask=mask_rs), irecv, dace%pio)
          do irs = 1, nrecv ; rsa(irs)%o_ch_i(i) = irecv(irs) ; enddo
        end do
        ! Channel info
        nsend = 0
        do irs = 1, n_set
          if (mask_rs(irs)) then
            nsend = nsend + rad_set(irs)%n_chan
          end if
        end do
        allocate(isend(nsend))
        nsend = 0
        do irs = 1, n_set
          if (mask_rs(irs)) then
            n = rad_set(irs)%n_chan
            isend(nsend+1:nsend+n) = rad_set(irs)%chan(1:n)
            nsend = nsend + n
          end if
        end do
        nrecv = 0
        if (dace%lpio) then
          do irs = 1, n_rs_all
            nrecv = nrecv + rsa(irs)%n_chan
          end do
        end if
        deallocate(irecv)
        allocate(irecv(nrecv))
        call p_gatherv(isend, irecv, dace%pio)
        if (dace%lpio) then
          nrecv = 0
          do irs = 1, n_rs_all
            n = rsa(irs)%n_chan
            allocate(rsa(irs)%chan(n))
            rsa(irs)%chan(1:n) = irecv(nrecv+1:nrecv+n)
            nrecv = nrecv + n
          end do
        end if
        ! iflag info
        allocate(mask_iflag(n_set)) ; mask_iflag = .false.
        do irs = 1, n_set
          mask_iflag(irs) = mask_rs(irs) .and. associated(rad_set(irs)%iflag)
        end do
        n_iflag = count(mask_iflag(1:n_set))
        l_iflag = (p_sum(n_iflag) > 0)
        if (l_iflag) then
          nsend = 0
          do irs = 1, n_set
            if (mask_rs(irs)) then
              n = rad_set(irs)%n_chan
              if (mask_iflag(irs)) then
                isend(nsend+1:nsend+n) = rad_set(irs)%iflag(1:n)
              else
                isend(nsend+1:nsend+n) = 0
              end if
              nsend = nsend + n
            end if
          end do
          call p_gatherv(isend, irecv, dace%pio)
          if (dace%lpio) then
            nrecv = 0
            do irs = 1, n_rs_all
              n = rsa(irs)%n_chan
              allocate(rsa(irs)%iflag(n))
              rsa(irs)%iflag(1:n) = irecv(nrecv+1:nrecv+n)
              nrecv = nrecv + n
            end do
          end if
        end if

        ! Unify rsa on dace%pio
        if (dace%lpio) then
          allocate(rsa_u(n_rs_all))
          n_rs = 0
          rs_all_loop: do i = 1, n_rs_all
            rsa_o => rsa(i)
            irs = -1
            do j = 1, n_rs
              if (rsa_u(j)%satid == rsa_o%satid .and. rsa_u(j)%grid == rsa_o%grid) then
                irs = j
                exit
              end if
            end do
            if (irs <= 0) then
              ! Copy rsa_o to unified rsa_
              n_rs = n_rs + 1
              rsa_ => rsa_u(n_rs)
              allocate(rsa_%instr(m_instr),rsa_%o_ch_i(m_instr),rsa_%n_ch_i(m_instr),rsa_%chan(m_chan))
              if (l_iflag) allocate(rsa_%iflag(m_chan))
              n = rsa_o%n_instr
              rsa_%satid                  = rsa_o%satid
              rsa_%grid                   = rsa_o%grid
              rsa_%opt_vars               = rsa_o%opt_vars
              rsa_%n_instr                = rsa_o%n_instr
              rsa_%n_chan                 = rsa_o%n_chan
              rsa_%instr  (1:n)           = rsa_o%instr  (1:n)
              rsa_%o_ch_i (1:n)           = rsa_o%o_ch_i (1:n)
              rsa_%n_ch_i (1:n)           = rsa_o%n_ch_i (1:n)
              rsa_%chan   (1:rsa_%n_chan) = rsa_o%chan   (1:rsa_%n_chan)
              if (l_iflag) rsa_%iflag  (1:rsa_%n_chan) = rsa_o%iflag  (1:rsa_%n_chan)
            else
              rsa_ => rsa_u(irs)
              rsa_%opt_vars = ior(rsa_%opt_vars, rsa_o%opt_vars)
              do j = 1, rsa_o%n_instr
                ! Unify instruments
                ii = -1
                do k = 1, rsa_%n_instr
                  if (rsa_%instr(k) == rsa_o%instr(j)) then
                    ii = k
                    exit
                  end if
                end do
                if (ii <= 0) then
                  ! Instrument not contained so far -> copy it
                  rsa_%n_instr = rsa_%n_instr + 1
                  ii = rsa_%n_instr
                  rsa_%instr(ii)  = rsa_o%instr(j)
                  rsa_%n_ch_i(ii) = rsa_o%n_ch_i(j)
                  if (ii > 1) rsa_%o_ch_i(ii) = rsa_o%o_ch_i(ii-1) + rsa_o%n_ch_i(ii-1)
                  rsa_%chan(rsa_%o_ch_i(ii)+1:rsa_%o_ch_i(ii)+rsa_%n_ch_i(ii)) = &
                       rsa_o%chan(rsa_o%o_ch_i(ii)+1:rsa_o%o_ch_i(ii)+rsa_o%n_ch_i(ii))
                  if (l_iflag) then
                    rsa_%iflag(rsa_%o_ch_i(ii)+1:rsa_%o_ch_i(ii)+rsa_%n_ch_i(ii)) = &
                         rsa_o%iflag(rsa_o%o_ch_i(ii)+1:rsa_o%o_ch_i(ii)+rsa_o%n_ch_i(ii))
                  end if
                else
                  ! Instrument present already -> unify channels
                  do k = rsa_o%o_ch_i(j)+1, rsa_o%o_ch_i(j)+rsa_o%n_ch_i(j)
                    ic = rsa_%o_ch_i(ii)+rsa_%n_ch_i(ii)+1
                    do m = rsa_%o_ch_i(ii)+1, rsa_%o_ch_i(ii)+rsa_%n_ch_i(ii)
                      if (rsa_o%chan(k) == rsa_%chan(m)) then
                        ic = 0
                        exit
                      elseif  (rsa_o%chan(k) < rsa_%chan(m)) then
                        ic = m
                      end if
                    end do
                    if (ic > 0) then
                      ! Channel not contained so far -> copy it
                      rsa_%chan(ic+1:rsa_%n_chan+1) = rsa_%chan(ic:rsa_%n_chan)
                      rsa_%chan(ic)   = rsa_o%chan(k)
                      rsa_%n_ch_i(ii) = rsa_%n_ch_i(ii) + 1
                      rsa_%n_chan     = rsa_%n_chan + 1
                      rsa_%o_ch_i(ii+1:rsa_%n_instr) = rsa_%o_ch_i(ii+1:rsa_%n_instr) + 1
                      if (l_iflag) then
                        rsa_%iflag(ic+1:rsa_%n_chan+1) = rsa_%iflag(ic:rsa_%n_chan)
                        rsa_%iflag(ic)   = rsa_o%iflag(k)
                      end if
                    end if
                  end do
                end if
              end do
            end if
          end do rs_all_loop
        end if

        if (allocated(rsa)) then
          call destruct_rsa(rsa)
          deallocate(rsa)
        end if

        call p_bcast(n_rs, dace%pio)
        if (.not.dace%lpio) allocate(rsa_u(n_rs))
        do irs = 1, n_rs
          call bcast_rad_set_aux(rsa_u(irs), dace%pio)
          rsa_ => rsa_u(irs)
          if (.not.dace%lpio) then
            if (associated(rsa_%instr )) allocate(rsa_%instr (rsa_%n_instr))
            if (associated(rsa_%n_ch_i)) allocate(rsa_%n_ch_i(rsa_%n_instr))
            if (associated(rsa_%o_ch_i)) allocate(rsa_%o_ch_i(rsa_%n_instr))
            if (associated(rsa_%chan  )) allocate(rsa_%chan  (rsa_%n_chan ))
            if (associated(rsa_%iflag )) allocate(rsa_%iflag (rsa_%n_chan ))
          end if
          call p_bcast(rsa_%instr (1:rsa_%n_instr),dace%pio)
          call p_bcast(rsa_%n_ch_i(1:rsa_%n_instr),dace%pio)
          call p_bcast(rsa_%o_ch_i(1:rsa_%n_instr),dace%pio)
          call p_bcast(rsa_%chan  (1:rsa_%n_chan ),dace%pio)
          if (l_iflag) &
               call p_bcast(rsa_%iflag(1:rsa_%n_chan ),dace%pio)

        end do

        ! Destruct old (non-unified) rad_set entries (set up by read_fdbk_rad)
        do i = 1, n_set
          if (mask_rs(i)) call destruct(rad_set(i))
        end do
        n_set = n_set - count(mask_rs)

        ! Build new (unified) rad_set entries
        do irs = 1, n_rs
          rsa_ => rsa_u(irs)
          n_set = n_set + 1
          rs   => rad_set(n_set)
          rs%id             = n_set
          rs%source         = 'fdbk'
          rs%satid          = rsa_%satid
          rs%gopts%opt_vars = rsa_%opt_vars
          rs%grid           = rsa_%grid
          rs%grid_wmo       = instr_rttov(rs%grid, rs%satid)
          rs%n_instr        = rsa_%n_instr
          rs%n_chan         = rsa_%n_chan
          n = rs%n_instr
          rs%n_sens      = n
          rs%instr (1:n) = rsa_%instr (1:n)
          rs%n_ch_i(1:n) = rsa_%n_ch_i(1:n)
          rs%o_ch_i(1:n) = rsa_%o_ch_i(1:n)
          do i = 1, n
            rs%instr_wmo(i) = instr_rttov(rs%instr(i), rs%satid)
          end do
          n = rs%n_chan
          allocate(rs%chan(n))
          rs%chan  (1:n) = rsa_%chan  (1:n)
          if (l_iflag) then
            if (any(rsa_%iflag(1:n) /= 0)) then
              allocate(rs%iflag(n))
              rs%iflag(1:n) = rsa_%iflag(1:n)
            end if
          end if

          if (dace%lpio) call print_rad_set(rs, hint=n_set, header='rad_set set up by '//proc)
        end do

        if (allocated(rsa_u)) then
          call destruct_rsa(rsa_u)
          deallocate(rsa_u)
        end if

        ! Link observations to new rad_set entries
        do ib = 1, size(obs)
          o => obs(ib)
          do is = 1, o% n_spot
            spt => o% spot(is)
            if (spt% hd% obstype /= OT_RAD) cycle
            call link_tovs_rs(spt, o)
          end do
        end do
      end if

      ! --- rad_set read from namelist ---
      n_fdbk = p_max(n_fdbk)
      if (n_fdbk == 0) RETURN
      ! Consistency check
      i = p_max(n_set)
      j = p_min(n_set)
      if (i /= j) then
        write(0,*) i,j,n_set
        call finish(proc,'n_set not equal on all PEs')
      end if
      ! communicate values that might have changed here
      do irs = 1, n_set
        rs => rad_set(irs)
        if (rs%source == 'nml') then
          rs%id             = p_max(rs%id            )
          rs%gopts%opt_vars = p_ior(rs%gopts%opt_vars)
          l_iflag = p_or(associated(rs%iflag))
          if (l_iflag) then
            n = rs%n_chan
            if (.not.associated(rs%iflag)) then
              allocate(rs%iflag(n)) ; rs%iflag = 0
            end if
            rs%iflag(1:n) = p_max(rs%iflag(1:n))
          end if
        end if
      end do


    case default
      call finish(proc,'invalid task='//trim(task))
    end select

  contains

    elemental subroutine destruct_rsa(rsa)
      type(t_rad_set_aux), intent(inout) :: rsa

      if (associated(rsa%instr )) deallocate(rsa%instr )
      if (associated(rsa%n_ch_i)) deallocate(rsa%n_ch_i)
      if (associated(rsa%o_ch_i)) deallocate(rsa%o_ch_i)
      if (associated(rsa%chan  )) deallocate(rsa%chan  )

    end subroutine destruct_rsa

#undef  VECTOR
#undef  DERIVED
#define DERIVED type(t_rad_set_aux)
#define p_bcast_DERIVED bcast_rad_set_aux
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_bcast_DERIVED


  end subroutine read_fdbk_rad




!------------------------------------------------------------------------------
  subroutine calc_tovs_quality(tsk, qmode, e_bg, e_o, instrs, satids, tl, obs, oi, si, qual)
  integer        ,intent(in)                      :: tsk     ! 1:preparation, 2:calculation, 3:cleanup
  integer        ,intent(in)    ,optional         :: qmode   ! quality mode
  type(t_vector) ,intent(in)    ,optional         :: e_bg    ! background error
  type(t_vector) ,intent(in)    ,optional         :: e_o     ! observation error
  integer        ,intent(in)    ,optional         :: instrs     (:) ! target instrs to be prepared
  integer        ,intent(in)    ,optional         :: satids     (:) ! satids to be prepared
  type(t_ilev)   ,intent(in)    ,optional ,target :: tl         (:) ! channels to be included
  type(t_obs)    ,intent(inout) ,optional ,target :: obs (:) ! observations
  type(t_obs)    ,intent(in)    ,optional ,target :: oi
  type(t_spot)   ,intent(in)    ,optional ,target :: si
  real(wp)       ,intent(out)   ,optional         :: qual

    integer         ,save        :: n_qu
    real(wp)        ,save        :: qu_max   (m_rad_set) = 0
    integer         ,save        :: qu_instrs(m_rad_set) = 0
    integer         ,save        :: qu_mod   (m_rad_set) = 0

    character(len=180)           :: str
    integer                      :: ib
    integer                      :: is
    type(t_spot)    ,pointer     :: s   => null()     ! pointer to report
    type(t_rad_set) ,pointer     :: r   => null()
    type(t_obs)     ,pointer     :: o   => null()
    type(t_tovs)                 :: tovs
    integer         ,allocatable :: ind   (:)
    real(wp)                     :: val
    integer                      :: qm
    integer                      :: ii, iset
    integer                      :: igrid
    integer                      :: n_lev
    logical                      :: luse

    if (btest(tsk, 0)) then
      call prep_tovs_quality
    end if

    if (btest(tsk, 1)) then
      if (present(obs)) then
        do ib = 1, size(obs)
          o => obs(ib)
          if (o% pe /= dace% pe) cycle
          do is = 1, o% n_spot
            s => o% spot (is)
            if (s% hd% obstype /= OT_RAD) cycle   ! only treat radiances
            call tovs_quality(val)
            s% pcc = val
            o% body(s%o%i+1:s%o%i+s%o%n)% pcc = val
          end do
        end do
      elseif (present(oi) .and. present(si)) then
        if (.not.present(qual)) then
          call finish('calc_tovs_quality','qual parameter missing')
        end if
        s => si
        o => oi
        call tovs_quality(qual)
      else
        call finish('calc_tovs_quality','no obs. given for calculation')
      end if
    end if

  contains

    subroutine prep_tovs_quality
    integer              :: qu_iatms(n_set)
    logical              :: rs_mask(n_set)
    integer              :: i, j, iset, igrid, ii, k, ib, is, satid, iatms

      ! Preparations
      rs_mask(:) = .false.
      n_qu = 0
      do iset = 1, n_set
        ! Determine datasets, that shall be included
        r => rad_set(iset)
        do igrid = 1, r%n_instr
          if (r%instr(igrid) == r%grid) exit
        end do
        call l_use(rs_mask(iset), instr=r%grid_wmo, sat=r%satid,             &
             chns=r%chan(r%o_ch_i(igrid)+1:r%o_ch_i(igrid)+r%n_ch_i(igrid)), &
             iatms=iatms)
        if (.not.rs_mask(iset)) cycle
        !if (dace% lpio) print*,'thin rs',r%satid, r%grid

        ! Set default quality_mode if necessary
        if (r% gopts% quality_mode < 0) r% gopts% quality_mode = quality_mode
        if (r% gopts% quality_mode == 6 .and. .not.hss_instr(r% instr(1))) &
             call finish('calc_tovs_quality','quality_mode=6 only allowed for hyperspectral sounders')

        ! Determine instruments, that require a normalization
        qm = r% gopts% quality_mode
        if (present(qmode)) then
          if (qmode >= 0) qm = qmode
        end if
        if (qm == 1 .or. qm == 2) then
          igrid = r% grid                    ! do not consider mapped instruments
          ii = -1
          do k = 1, n_qu
            if (qu_instrs(k) == igrid .and. qu_mod(k) == qm) then
              ii = k
              exit
            end if
          end do
          if (ii <= 0) then
            n_qu = n_qu + 1
            qu_instrs(n_qu) = igrid
            qu_mod   (n_qu) = qm
            qu_max   (n_qu) = 0._wp
            qu_iatms (n_qu) = iatms
            if (qm == 2 .and. .not.(present(e_bg).and.present(e_o))) &
                 call finish('calc_tovs_quality','quality_mode=2 without e_bg and e_o is not valid.')
          end if
        end if
      end do

      ! Calculate maximum values for normalization
      if (n_qu > 0) then
        do ib = 1, size(obs)
          if (obs(ib)% pe /= dace% pe) cycle
          do is = 1, obs(ib)% n_spot
            s => obs(ib)% spot (is)
            if (s% hd% obstype /= OT_RAD) cycle   ! only treat radiances
            igrid = int(s% hd% grid_id)
            satid = int(s% hd% satid)
            iset  = set_indx(rad_set(1:n_set), satid=satid, grid=igrid)
            if (.not.rs_mask(iset)) cycle
            r => rad_set(iset)
            qm = r%gopts%quality_mode
            if (present(qmode)) then
              if (qmode >= 0) qm = qmode
            end if
            allocate(ind(s%o%n))
            call l_use(luse, si=s, oi=obs(ib), n_lev=n_lev, ind=ind)
            if (luse) then
              do k = 1, n_qu
                if (qu_instrs(k) == igrid .and. qu_mod(k) == qm) then
                  val = 0._wp
                  do i = 1, n_lev
                    j = ind(i)
                    if (obs(ib)% body(j)% use% state > STAT_REJECTED) then
                      select case(qu_mod(k))
                      case(1)
                        val = val + 1
                      case(2)
                        val = val + (e_bg% s(ib)% x(j) / e_o% s(ib)% x(j))**2
                      end select
                    end if
                  end do
                  qu_max(k) = max(qu_max(k), val)
                end if
              end do
            end if
            deallocate(ind)
          end do
        end do
        if (dace% lpio) then
          write(*,*) 'Maximum TOVS qualities:'
          write(*,'(T3,A,T10,A)') 'instr','max(quality)'
        end if
        do i = 1, n_qu
          qu_max(i) = p_max(qu_max(i))
          if (dace% lpio) write(*,'(T6,I3.3,T10,F12.7)') qu_instrs(i),qu_max(i)
          do j = 1, n_qu
            if ( qu_iatms(i) == 1 .and. qu_instrs(j)==3  .or. &
                 qu_iatms(j) == 1 .and. qu_instrs(i)==3  .or. &
                 qu_iatms(i) == 2 .and. qu_instrs(j)==15 .or. &
                 qu_iatms(j) == 2 .and. qu_instrs(i)==15) then
              val       = qu_max(i)
              qu_max(i) = max(qu_max(i), qu_max(j))
              qu_max(i) = p_max(qu_max(i))
              if (dace% lpio .and. val /= qu_max(i)) write(*,'(T7,"-> ",F12.7)') qu_max(i)
            end if
          end do
        end do
      endif
    end subroutine prep_tovs_quality

    subroutine tovs_quality(quality)
      real(wp), intent(out) :: quality
      integer :: i, j, k

      ! s, o must be set
      quality = 0._wp
      igrid =  int(s% hd% grid_id)
      iset  =  set_indx(rad_set(1:n_set), satid=int(s% hd% satid), grid=igrid)
      r     => rad_set(iset)
      qm = r% gopts% quality_mode
      if (present(qmode)) then
        if (qmode >= 0) qm = qmode
      end if
      allocate(ind(s%o%n))
      call l_use(luse, si=s, oi=o, n_lev=n_lev, ind=ind)
      if (.not.luse) then
        deallocate(ind)
        return
      end if

      select case (qm)
      case(0)
        quality = 0
      case(1)
        val = count(o% body(ind(1:n_lev))% use% state > STAT_REJECTED)
      case(2)
        if (.not.(present(e_bg).and.present(e_o))) &
             call finish('calc_tovs_quality','quality_mode=2 without e_bg and e_o is not valid.')
        val = 0._wp
        do i = 1, n_lev
          j = ind(i)
          if (o% body(j)% use% state > STAT_REJECTED) then
            val = val + (e_bg% s(o%ibox)% x(j) / e_o% s(o%ibox)% x(j))**2
          end if
        end do
      case(3)
        if ( any(btest(o%body(ind(1:n_lev))%use%flags, CHK_CLOUD)) .or. &
             btest(s% stlsf, SUR_LAND)                             .or. &
             btest(s% stlsf, SUR_MWSURF) ) then
          quality = 0
        else
          quality = 100
        end if
      case(4)
        if ( any(btest(o%body(ind(1:n_lev))%use%flags, CHK_CLOUD)) .or. &
             any(btest(o%body(ind(1:n_lev))%use%flags, CHK_SURF )) .or. &
             btest(s% stlsf, SUR_LAND)                             .or. &
             btest(s% stlsf, SUR_MWSURF)) then
          quality = 0
        else
          quality = 100
        end if
      case(5)
        call load(o, s, tovs, tovs_io=TTOVS_CI)
        if ( any(btest(o%body(ind(1:n_lev))%use%flags, CHK_CLOUD))       .or. &
             any(btest(o%body(ind(1:n_lev))%use%flags, CHK_SURF ))       .or. &
             any(btest(o%body(ind(1:n_lev))%use%flags, CHK_DATASET) .and.     &
                 .not.btest(r%flag(tovs%ci(1:k)),USE_PASSIVE) )             .or. &
             btest(s% stlsf, SUR_LAND)                                   .or. &
             btest(s% stlsf, SUR_MWSURF)) then
          quality = 0
        else
          quality = 100
        end if
        call destruct(tovs)
      case(6)
        ! Information on cloud check results is kept in s%ppc for IASI
        if (quality == 100._wp) then
          if ( any(btest(o%body(ind(1:n_lev))%use%flags, CHK_SURF )) .or. &
               btest(s% stlsf, SUR_LAND)                                       .or. &
               btest(s% stlsf, SUR_MWSURF)) quality = 0._wp
        else
          quality = 0._wp
        end if
      case(7)
        ! Nothing to do. Quality was calculated in decr_quality
        quality = s%pcc
      case(8)
        call load(o, s, tovs, tovs_io=0)
        quality = 100._wp * cos(abs(tovs%saza) * d2r)
        call destruct(tovs)
      case(9) ! hybrid 4 (major) + 8 (minor)
        if ( any(btest(o%body(ind(1:n_lev))%use%flags, CHK_CLOUD)) .or. &
             any(btest(o%body(ind(1:n_lev))%use%flags, CHK_SURF )) .or. &
             btest(s% stlsf, SUR_LAND)                             .or. &
             btest(s% stlsf, SUR_MWSURF)) then
          quality = 0._wp
        else
          quality = 50._wp
        end if
        call load(o, s, tovs, tovs_io=0)
        quality = quality + 50._wp * cos(abs(tovs%saza) * d2r)
        call destruct(tovs)
      case default
        call finish('calc_tovs_quality','unknown quality_mode')
      end select
      ! Normalization
      if (qm == 1 .or. qm == 2) then
        ii = -1
        do k = 1, n_qu
          if (qu_instrs(k) == igrid .and. qu_mod(k) == qm) then
            ii = k
            exit
          end if
        end do
        if (ii < 0) then
          write(str,'(I3.3,1x,I2,"(",100(I3.3,1x,I2,2x),")")') igrid, qu_instrs(1:n_qu)
          str = trim(str)//')'
          call finish('calc_tovs_quality','invalid instrument/quality_mode combination '//trim(str))
        end if
        if (qu_max(ii) > 0._wp) then
          quality = int(100. * val/qu_max(ii))
        else
          quality = 0
        end if
      end if
      deallocate(ind)

!       s%pcc = quality
!       call debug_spot(sp=s, obs=o, hint='calc_tovs_quality')
!       !!!!!!!!!!!!!!!!!!!!!!
!       if (present(tl)) then
!         o%body(ind(1:n_lev))% pcc = quality
!       end if


    end subroutine tovs_quality

    subroutine l_use(luse, instr, sat, chns, si, oi, n_lev, ind, iatms)
    logical      ,intent(out)                   :: luse
    integer      ,intent(in)  ,optional ,target :: instr
    integer      ,intent(in)  ,optional ,target :: sat
    integer      ,intent(in)  ,optional ,target :: chns(:)
    type(t_spot) ,intent(in)  ,optional ,target :: si
    type(t_obs)  ,intent(in)  ,optional         :: oi
    integer      ,intent(out) ,optional         :: n_lev
    integer      ,intent(out) ,optional         :: ind(:)
    integer      ,intent(out) ,optional         :: iatms ! 1 for ATMS temp chans, 2 for ATMS hum. chans, 0 otherwise

      type(t_ilev)  ,target  :: tl_dum(0)
      type(t_ilev)  ,pointer :: tp(:) => null()
      integer       ,pointer :: ii    => null()
      integer       ,pointer :: satid => null()
      integer                :: ni, id

      ! Preparations
      luse = .true.

      if (present(instr)) then
        ii => instr
      elseif (present(si)) then
        ii => si% sttyp
      else
        call finish('calc_tovs_quality/l_use','Either instr or si parameter required.')
      end if
      if (present(sat)) then
        satid => sat
      elseif (present(si)) then
        satid => si% ident
      else
        call finish('calc_tovs_quality/l_use','Either satid or si parameter required.')
      end if
      id = rttov_instr(ii, satid)

      ! Check, whether instr/sat/channels match to requirements
      if (present(instrs)) then
        if (instrs(1) >= 0) then
          if (.not.any(instrs(:) == id)) then
            luse = .false.
          end if
        end if
      end if
      if (present(satids)) then
        if (satids(1) >= 0) then
          if (.not.any(satids(:) == satid)) then
            luse = .false.
          end if
        end if
      end if
      if (.not.luse) then
        if (present(n_lev)) n_lev  = 0
        if (present(ind  )) ind(:) = 0
        if (present(iatms)) iatms  = 0
        return
      end if

      if (present(ind) .or. present(tl) .or. present(n_lev) .or. present(iatms)) then
        if (present(tl)) then
          tp => tl
        else
          tp => tl_dum
        end if
        call match_ilev(tp, ni, instr=instr, levs=chns, si=si, oi=oi, ind=ind, iatms=iatms)
        if (present(n_lev)) n_lev = ni
        luse = (ni > 0)
      end if

    end subroutine l_use

  end subroutine calc_tovs_quality


!------------------------------------------------------------------------------
  real(wp) function get_l2c_max(lat, l2c_max_trop, l2c_max_midlat, l2c_max_polar)
    real(wp), intent(in) :: lat
    real(wp), intent(in) :: l2c_max_trop
    real(wp), intent(in) :: l2c_max_midlat
    real(wp), intent(in) :: l2c_max_polar
    real(wp), parameter  :: pi2 = acos(-1._wp)*0.5_wp
    real(wp)             :: alat, l2c_max_p, tr2ml_1, tr2ml_2, ml2pl_1, ml2pl_2
    integer              :: i

    tr2ml_1 = l2c_max_tr2ml - 0.5_wp * l2c_max_tr2ml_width
    tr2ml_2 = l2c_max_tr2ml + 0.5_wp * l2c_max_tr2ml_width
    ml2pl_1 = l2c_max_ml2pl - 0.5_wp * l2c_max_ml2pl_width
    ml2pl_2 = l2c_max_ml2pl + 0.5_wp * l2c_max_ml2pl_width

    alat = abs(lat)
    if (alat > ml2pl_2) then
      l2c_max_p = l2c_max_polar
    elseif (alat > ml2pl_1) then
      l2c_max_p = l2c_max_polar + &
           sin((ml2pl_2-alat)/l2c_max_ml2pl_width*pi2)**2 * (l2c_max_midlat - l2c_max_polar)
    elseif (alat > tr2ml_2) then
      l2c_max_p = l2c_max_midlat
    elseif (alat > tr2ml_1) then
      l2c_max_p = l2c_max_midlat + &
           sin((tr2ml_2-alat)/l2c_max_tr2ml_width*pi2)**2 * (l2c_max_trop   - l2c_max_midlat)
    else
      l2c_max_p = l2c_max_trop
    end if
    l2c_max_p = l2c_max_p

    ! convert l2c_max in rttov level
    do i = 1, size(preshPa)
      if (preshPa(i) > l2c_max_p) exit
    end do
    if (i < 2 .or. i > size(preshPa)) call finish('get_l2c_max','invalid l2c_max')
    get_l2c_max = ( (i-1) * abs(l2c_max_p - preshPa(i))  + &
                    i     * abs(l2c_max_p - preshPa(i-1))) / abs(preshPa(i) - preshPa(i-1))

  end function get_l2c_max

!==============================================================================
  subroutine read_tovs_nml
  !-----------------------------------------------
  ! read namelists /TOVS_OBS/ and /TOVS_OBS_CHAN/
  !-----------------------------------------------
    character(len=80) :: header
    integer           :: ierr_pos
    integer           :: ierr_read
    integer           :: ierr, i, j, i_nml, i_read
    logical           :: first   ! first occurence of /TOVS_OBS_CHAN/ in file
    integer           :: nf      ! actual number of feedback files
    type(t_trg_clim_trend), pointer :: tct
    !---------------------------
    ! 1. read namelist TOVS_OBS
    !---------------------------
    ! 1.a set defaults
    !------------------
    ! namelist TOVS_OBS
    netcdf_path           = obsinput
    data_path             = ''
    feedbk_files(:)       = ''
    read1pe               = .true.   ! Read coeffs.only on I/O PE
                                     ! (Only effective with -D_RTIFC_DISTRIBCOEF)
    rttov_version         = -1       ! rttov version to use
    rttov_mult_prof       = .true.
    e_bg_ts_sea           = 1._wp    ! surf.temp. bg.error sea
    e_bg_ts_land          = 3._wp    ! surf.temp. bg.error land
    e_bg_drts_land        = 1.5_wp   ! retrieved surf. temp. bg. error land
    e_bg_ts_ice           = 3._wp    ! surf.temp. bg.error ice
    e_bg_t_top            = 3._wp    ! bg.error temperature top level
    e_bg_lnq_top          = 1._wp    ! bg.error ln(q) top level
    flg_sur               =  0       ! 1dvar surface flag
    lhum_dum_ana          = .true.   ! dummy humidity above analysis levels
    p_top                 = -1._wp   ! model top level pressure
    monitor_prof          = 0        ! monitor profiles 1:ascii 2:netcdf 4:psas-ana
    mon_ana_prof          = 0        !                  8:H    16:Bii
    monitor_prof_split    = 0
    mon_ana_prof_split    = 0
    monitor_prof_sat      = -1
    mon_ana_prof_sat      = -1
    monitor_prof_grid     = -1
    mon_ana_prof_grid     = -1
    monitor_prof_wrbxf    = 1._wp
    mon_ana_prof_wrbxf    = 1._wp

    lBii                  = .false.  ! bit 16 set in mon_ana_prof
    trace_mem_time        = .false.  ! print memory & time usage
!   param_fil_amsub       = 'amsub_ch18_clear_thresh.asc' ! AMSU-B clchk parameter file name
    param_fil_amsub       = ''                            ! AMSU-B clchk parameter file name
    clchk_amsub           = 0                             ! AMSU-B cloud check
    amsub_ec_bnd          = param_amsub% ec_bnd           ! ECMWF cloud check bounds
    amsub_delta_ch20_ch18 = param_amsub% delta_ch20_ch18  ! AMSU-B ch(20)-ch(18) threshold
    hd_id_vers            = 2
    use_ch_blcklst        = .false. ! use preprocessing channel blacklist
    req_nmlst             = .true.  ! finish 3dvar if sat input file without corresp. namelist
    filter_ch_stat        = -1
    surf_class_vers       =  0
    snw_frc_mode          = -1
    snw_frc_par           =  0._wp
    no_ps_dep             = .false. ! no dependence on ps (fix RTTOV10 adjoint)
!    no_t_dep              = 0       ! no dependence on top t-levels
    rt_min_od             = 1.e-5_wp! parameter min_od in RTTOV10
    rtt10_min_od          = 1.e-5_wp! parameter min_od in RTTOV10
    cld_top_def           = 500._wp ! first guess cloud top
    cld_frc_def           = 0._wp   ! first guess cloud fraction
    cld_top_e             = 0._wp   ! cloud top background error
    cld_frc_e             = 0._wp   ! cloud fraction background error
    bound_constr          = 1       ! cloud top/frc boundary constraint flag
    quality_mode          = quality_mode_def
    l2c_type              = l2c_type_def
    l2c_god_corr          = .false.
    l2c_rel_lim           = l2c_rel_lim_def
    l2c_abs_lim           = l2c_abs_lim_def
    l2c_use_rad           = l2c_use_rad_def
    l2c_max               = l2c_max_def
    l2c_max_trop          = l2c_max_trop_def
    l2c_max_midlat        = l2c_max_midlat_def
    l2c_max_polar         = l2c_max_polar_def
    surf_infl_mode        = surf_infl_mode_def
    max_surf_infl         = max_surf_infl_def
    d_stemp               = d_stemp_def
    d_emiss               = d_emiss_def
    l_max_surf_hgt        = l_max_surf_hgt_def
    qflag_bitmask         = 0
    pl_log                = .true. ! use log(p) for mean, stdev
    pl_method             = 2      ! method for pressure level estimation
    pl_e_bg_t             = -1._wp ! temp.    background error for pl_method
    pl_e_bg_rh            = -1._wp ! rel.hum. background error for pl_method
    pl_wd_thresh          = 0.1_wp
    pl_vis                = -1._wp
    plev_min              = 0.0_wp ! minimum value for plevel
    use_reff              = .false.
    scale_qi              = .false.
    min_dia_qc            = 5.0_wp
    min_dia_qi            = 10.0_wp
#if (_RTTOV_VERSION >= 12)
    rt_ir_emis_mod        = 2      ! IREMIS
#else
    rt_ir_emis_mod        = 1      ! ISEM
#endif
    rt_mw_emis_mod        = 5      ! FASTEM version
    ! force_fastem: please see the note in prep_rttov_prof
    select case(method)
    case('MEC', 'GMESTAT', 'VERI_ENS', 'FC_SENS')
      ! For mec the emissivity read from the feedback file should be used.
      force_fastem        = .false.
    case default
      ! Force to call FASTEM for each RTTOV call.
      force_fastem        = .true.
    end select
#if (_RTTOV_VERSION >= 12)
    rt_humi               = gas_unit_specconc
#else
    rt_humi               = 0      ! rttov10 humidity
#endif
    rt_use_q2m            = .false.
    rt_do_lambertian      = .false.
    rt_salinity           = 0._wp
    rt_levs_exact         = .true.
    rt_hard_limits        = 2
    app_reg_lims          = .true.
    fix_hgpl              = 0
    chk_reg_lims          = empty_rtstat_opts
    chk_opdep             = empty_rtstat_opts
    chk_ps                = empty_rtstat_opts
    opdep_thresh          = 1._wp
    opdep_min_transm      = 0._wp
    use_hum_top           = 0
    warn_gen_hum          = .false.
    chk_god               = empty_rtstat_opts
    god_thresh            = 99._wp
    god_par_file          = ''
    wr_god                = .false.
    wr_opdep              = 0
    ! tracegases
    trg_file              = ''
    trg_clim_file         = ''
    trg_clim_trend        = empty_trg_clim_trend
    trg_hist_file         = 'bias_RAD_TRGHIST._YYYYMMDDHHMM_'
    trg_hist_inidate      = '2099010100'
    ! obsolete
    max_scan              = 1000
    flg_prc               = 15       ! 1=data,2=min,4=sur,8=cld
    flg_prc_amsua         =  0       ! 1=data,2=min,4=sur,8=cld
    flg_prc_amsub         =  0       ! 1=data,2=min,4=sur,8=cld
    flg_prc_hirs          =  0       ! 1=data,2=min,3=sur,8=cld

    if (dace% lpio) then
      call position_nml ('TOVS_OBS', status=ierr)
    end if
    call p_bcast (ierr, dace% pio)
    if (ierr /= POSITIONED) then
      if (dace% lpio) then
        write (6,'(a)')   repeat('-',79)
        write (6,'(a)') ' Namelist /TOVS_OBS/ not present'
      end if
      return
    end if
    if (dace% lpio) then
      if (ierr==POSITIONED) then
        read (nnml ,nml=TOVS_OBS)
      end if
      !--------------------------
      ! 1.b consistency checks
      !--------------------------
      if (max_scan > 1000) max_scan = 1000

      ! backwards compatibility of *_min_od
      if (rtt10_min_od /= 1.e-5_wp) then
         if (rt_min_od == 1.e-5_wp) rt_min_od = rtt10_min_od
      end if
      !---------------------------------------------------------------
      ! temperature vs. humidity weights for pressure level estimation
      !---------------------------------------------------------------
      if (pl_method == 4) then
        if (pl_e_bg_t  < 0._wp) pl_e_bg_t  =  1._wp  ! relative weight
        if (pl_e_bg_rh < 0._wp) pl_e_bg_rh =  1._wp
      else
        if (pl_e_bg_t  < 0._wp) pl_e_bg_t  = 0.5_wp  ! temperature error
        if (pl_e_bg_rh < 0._wp) pl_e_bg_rh = 0.1_wp  ! humidity error
      endif
      !------------------------
      ! monitoring restrictions
      !------------------------
      if (iand (mon_ana_prof, MP_H) == 0) &
        mon_ana_prof = mon_ana_prof - iand (mon_ana_prof, MP_B)
      monitor_prof   = monitor_prof - iand (monitor_prof, MP_B)
      monitor_prof   = monitor_prof - iand (monitor_prof, MP_HBH) ! only use this option for active obs
      !-------------
      ! RTTOV checks
      !-------------
      j = 1
      do i = 2, size(opdep_thresh)
        if (opdep_thresh(i) > 0._wp) then
          if (j > 0) opdep_thresh(i) = opdep_thresh(j)
        else
          j = i
        end if
      end do
      j = 1
      do i = 2, size(god_thresh)
        if (god_thresh(i) < 0._wp .or. god_thresh(i) >= 1._wp) then
          if (j > 0) god_thresh(i) = god_thresh(j)
        else
          j = i
        end if
      end do
      !-------------------
      ! derived quantities
      !-------------------
      lBii = iand (mon_ana_prof, MP_B) /= 0
      !------------------------
      ! god_par_file
      !------------------------
      if (god_par_file /= '') then
        god_par_file = adjustl(god_par_file)
        if (scan(god_par_file, '/') == 0) god_par_file = path_file (data, god_par_file)
      end if
      n_filter_ch                   = count(filter_ch_stat(:) >= 0)
      filter_ch_stat(1:n_filter_ch) = pack(filter_ch_stat(:), mask = filter_ch_stat(:) >= 0)

      ! Translate use_hum_top (which is deprecated) into other options
      if (use_hum_top > 0) then
        fg_prof_top = 1
        if (use_hum_top <= 1) e_bg_lnq_top = 0._wp
      end if
      p_top_ext = max(p_top_ext, 0.5_wp)

      dynret_w_pref = min(dynret_w_pref, 1._wp)

      !---------
      ! printout
      !---------
      nf = count (feedbk_files /= '')
      write (6,'(a)')         repeat('-',79)
      write (6,'(a)')       ' Namelist /TOVS_OBS/ read:'
      write (6,'( )')
      write (6,'(a)')       ' ! Input'
      write (6,'(a,a)')     ' netcdf_path             = ', trim (netcdf_path)
      write (6,'(a,a)')     ' data_path               = ', trim (data_path)
      write (6,'(a,a)')     ' feedbk_files            = ', trim (feedbk_files(1))
      do i   = 2, nf
        write (6,'(a,a)')   '                           ', trim (feedbk_files(i))
      end do
      write (6,'(a,l1)')    ' read1pe                 = ', read1pe
      write (6,'( )')
      write (6,'(a,i6)')    ' flg_sur                 =' , flg_sur
      write (6,'(a,l1)')    ' use_ch_blcklst          = ', use_ch_blcklst
      write (6,'(a,l1)')    ' req_nmlst               = ', req_nmlst
      write (6,'(a,20(1x,i2))') ' filter_ch_stat          =' , filter_ch_stat(1:n_filter_ch)
      write (6,'( )')
      write (6,'(a)')       ' ! Cloud check'
      write (6,'(a,a)')     ' param_fil_amsub         = ', trim (param_fil_amsub)
      write (6,'(a,i6)')    ' clchk_amsub             = ', clchk_amsub
      write (6,'(a,2f9.2)') ' amsub_ec_bnd            = ', amsub_ec_bnd
      write (6,'(a,f9.2)')  ' amsub_delta_ch20_ch18   = ', amsub_delta_ch20_ch18
      write (6,'( )')
      write (6,'(a)')       ' ! Cloud analysis'
      write (6,'(a,f9.2)')  ' cld_top_def             =' , cld_top_def
      write (6,'(a,f9.2)')  ' cld_frc_def             =' , cld_frc_def
      write (6,'(a,f9.2)')  ' cld_top_e               =' , cld_top_e
      write (6,'(a,f9.2)')  ' cld_frc_e               =' , cld_frc_e
      write (6,'(a,i6)')    ' bound_constr            =' , bound_constr
      write (6,'(a,f9.2)')  ' cld_top_def_v           =' , cld_top_def_v
      write (6,'(a,f9.2)')  ' cld_frc_def_v           =' , cld_frc_def_v
      write (6,'(a,f9.2)')  ' cld_top_e_v             =' , cld_top_e_v
      write (6,'(a,f9.2)')  ' cld_frc_e_v             =' , cld_frc_e_v
      write (6,'( )')
      write (6,'(a)')       ' ! Surface classification'
      write (6,'(a,i2)')    ' surf_class_vers         =' , surf_class_vers
      write (6,'(a,f9.2)')  ' e_bg_ts_sea             =' , e_bg_ts_sea
      write (6,'(a,f9.2)')  ' e_bg_ts_land            =' , e_bg_ts_land
      write (6,'(a,f9.2)')  ' e_bg_drts_land          =' , e_bg_drts_land
      write (6,'(a,f9.2)')  ' e_bg_ts_ice             =' , e_bg_ts_ice
      write (6,'(a,i2)')    ' snw_frc_mode            =' , snw_frc_mode
      write (6,'(a,5(1x,G11.4))') ' snw_frc_par             =' , snw_frc_par
      write (6,'(a,f9.2)')  ' dynret_w_pref           =' , dynret_w_pref
      write (6,'(a,i2)')    ' dynret_avg              =' , dynret_avg
      write (6,'(a,l1)')    ' check_recalc_stype      = ', check_recalc_stype
      write (6,'( )')
      write (6,'(a)')       ' ! Top of atmosphere'
      write (6,'(a,f9.2)')  ' e_bg_t_top              =' , e_bg_t_top
      write (6,'(a,f9.2)')  ' e_bg_lnq_top            =' , e_bg_lnq_top
      write (6,'(a,i1)')    ' fg_prof_top             = ', fg_prof_top
      write (6,'(a,i1)')    ' use_hum_top(deprecated) = ', use_hum_top
      write (6,'(a,f9.2)')  ' p_top                   =' , p_top
      write (6,'(a,f9.2)')  ' p_top_ext               =' , p_top_ext
      write (6,'(a,f9.2)')  ' p_blend_ext             =' , p_blend_ext
      write (6,'(a,i2)')    ' nlev_ext                = ', nlev_ext
      write (6,'( )')
      write (6,'(a)')       ' ! Humidity options'
      write (6,'(a,l1)')    ' lhum_dum_ana            = ', lhum_dum_ana
      write (6,'(a,l1)')    ' warn_gen_hum            = ', warn_gen_hum
      write (6,'( )')
      write (6,'(a)')       ' ! Channel height (l2c) checks (do not confuse with l2c for McNally-Watts)'
      write (6,'(a,f9.2)')  ' l2c_max_tr2ml           = ', l2c_max_tr2ml
      write (6,'(a,f9.2)')  ' l2c_max_tr2ml_width     = ', l2c_max_tr2ml_width
      write (6,'(a,f9.2)')  ' l2c_max_ml2pl           = ', l2c_max_ml2pl
      write (6,'(a,f9.2)')  ' l2c_max_ml2pl_width     = ', l2c_max_ml2pl_width
      write (6,'( )')
      write (6,'(a)')       ' ! t_rad_gopt and t_rad_iopt default values (see TOVS_OBS_CHAN namelists)'
      write (6,'(a,i6)')    ' flag_instr              = ', flag_instr
      write (6,'(a,i7)')    ' n_max_calc_k_def        =' , n_max_calc_k_def
      write (6,'(a,i7)')    ' n_max_prof_k_def        =' , n_max_prof_k_def
      write (6,'(a,i7)')    ' n_max_calc_y_def        =' , n_max_calc_y_def
      write (6,'(a,i7)')    ' n_max_prof_y_def        =' , n_max_prof_y_def
      write (6,'(a,i1)')    ' quality_mode            = ', quality_mode
      write (6,'(a,i1)')    ' l2c_type                = ', l2c_type
      write (6,'(a,l1)')    ' l2c_god_corr            = ', l2c_god_corr
      write (6,'(a,g13.6)') ' l2c_rel_lim             = ', l2c_rel_lim
      write (6,'(a,g13.6)') ' l2c_abs_lim             = ', l2c_abs_lim
      write (6,'(a,l1)')    ' l2c_use_rad             = ', l2c_use_rad
      write (6,'(a,f9.2)')  ' l2c_max                 = ', l2c_max
      write (6,'(a,f9.2)')  ' l2c_max_trop            = ', l2c_max_trop
      write (6,'(a,f9.2)')  ' l2c_max_midlat          = ', l2c_max_midlat
      write (6,'(a,f9.2)')  ' l2c_max_polar           = ', l2c_max_polar
      write (6,'(a,i1)')    ' surf_infl_mode          = ', surf_infl_mode
      write (6,'(a,f9.2)')  ' max_surf_infl           = ', max_surf_infl
      write (6,'(a,f9.2)')  ' d_stemp                 =' , d_stemp
      write (6,'(a,f9.2)')  ' d_emiss                 =' , d_emiss
      write (6,'(a,l1)')    ' l_max_surf_hgt          = ', l_max_surf_hgt
      write (6,'( )')
      write (6,'(a)')       ' ! RTTOV options'
      write (6,'(a,i6)')    ' rttov_version           =' , rttov_version
      write (6,'(a,i6)')    ' rttov_levels            =' , rttov_levels
      write (6,'(a,l1)')    ' rt13_clip_gas_opdep     = ', rt13_clip_gas_opdep
      write (6,'(a,l1)')    ' rt_levs_exact           = ', rt_levs_exact
      write (6,'(a,l1)')    ' rttov_mult_prof         = ', rttov_mult_prof
      write (6,'(a,i1)')    ' rt_ir_emis_mod          = ', rt_ir_emis_mod
      write (6,'(a,i1)')    ' rt_mw_emis_mod          = ', rt_mw_emis_mod
      write (6,'(a,l1)')    ' force_fastem            = ', force_fastem
      write (6,'(a,i1)')    ' rt_humi                 = ', rt_humi
      write (6,'(a,l1)')    ' rt_use_q2m              = ', rt_use_q2m
      write (6,'(a,l1)')    ' rt_do_lambertian        = ', rt_do_lambertian
      write (6,'(a,f9.2)')  ' rt_salinity             = ', rt_salinity
      write (6,'(a,e9.1)')  ' rt_min_od               =' , rt_min_od
      write (6,'(a,i1)')    ' rt_hard_limits          = ', rt_hard_limits
      write (6,'(a,l1)')    ' app_reg_lims            = ', app_reg_lims
      write (6,'(a,i1)')    ' fix_hgpl                = ', fix_hgpl
      write (6,'(a,a)')     ' god_par_file            = ', trim(god_par_file)
      write (6,'(a,l1)')    ' wr_god                  =' , wr_god
      write (6,'(a,l1)')    ' atlas_single_inst       = ', atlas_single_inst
      write (6,'(a,i1)')    ' rt_alloc_mode           = ', rt_alloc_mode
      write (6,'( )')
      write (6,'(a)')       ' ! Checks on RTTOV applicability'
      write (6,'(a,i1)')    ' chk_reg_lims%action_fg  = ', chk_reg_lims%action_fg
      write (6,'(a,i1)')    ' chk_reg_lims%action_ana = ', chk_reg_lims%action_ana
      write (6,'(a,*(25(1x,i1,:),/,27x))') ' chk_reg_lims%action     = ', chk_reg_lims%action
      write (6,'(a,i1)')    ' chk_reg_lims%action_mon = ', chk_reg_lims%action_mon
      write (6,'(a,i1)')    ' chk_reg_lims%op_na_bit  = ', chk_reg_lims%op_na_bit
      write (6,'(a,i6)')    ' chk_reg_lims%tsk        = ', chk_reg_lims%tsk
      write (6,'(a,i2)')    ' chk_reg_lims%nitout     = ', chk_reg_lims%nitout
      write (6,'(a,*(25(1x,i1,:),/,27x))') ' chk_opdep%action        = ', chk_opdep%action
      write (6,'(a,i1)')    ' chk_opdep%action_mon    = ', chk_opdep%action_mon
      write (6,'(a,i1)')    ' chk_opdep%op_na_bit     = ', chk_opdep%op_na_bit
      write (6,'(a,i6)')    ' chk_opdep%tsk           = ', chk_opdep%tsk
      write (6,'(a,*(5(e9.2,1x,:),/,27x))') ' opdep_thresh            = ', opdep_thresh
      write (6,'(a,*(5(f9.5,1x,:),/,27x))') ' opdep_min_transm        = ', opdep_min_transm
      write (6,'(a,*(25(1x,i1,:),/,27x))') ' chk_ps%action           = ', chk_ps%action
      write (6,'(a,i1)')    ' chk_ps%action_mon       = ', chk_ps%action_mon
      write (6,'(a,i1)')    ' chk_ps%op_na_bit        = ', chk_ps%op_na_bit
      write (6,'(a,i6)')    ' chk_ps%tsk              = ', chk_ps%tsk
#if (_RTTOV_VERSION >= 12)
      write (6,'(a,i1)')    ' chk_god%action_fg       = ', chk_god%action_fg
      write (6,'(a,i1)')    ' chk_god%action_ana      = ', chk_god%action_ana
      write (6,'(a,*(25(1x,i1,:),/,27x))') ' chk_god%action          = ', chk_god%action
      write (6,'(a,i1)')    ' chk_god%action_mon      = ', chk_god%action_mon
      write (6,'(a,i1)')    ' chk_god%op_na_bit       = ', chk_god%op_na_bit
      write (6,'(a,i6)')    ' chk_god%tsk             = ', chk_god%tsk
      write (6,'(a,*(5(e9.2,1x,:),/,27x))') ' god_thresh              = ', god_thresh
      write (6,'(a,i3)')    ' wr_opdep                = ' , wr_opdep
#endif
      write (6,'(a,i10)')   ' qflag_bitmask           = ' , qflag_bitmask
      write (6,'( )')
      write (6,'(a)')       ' ! Tracegases'
      write (6,'(a,a)')     ' trg_file                = ', trim(trg_file)
      write (6,'(a,a)')     ' trg_clim_file           = ', trim(trg_clim_file)
      do i = 1, size(trg_clim_trend)
        tct => trg_clim_trend(i)
        if (trim(tct%gas_name) /= '') then
          write (6,'(a,I1,a,a10,2x,a14,2x,e13.6,2x,e13.6)')     ' trg_clim_trend(',i,')       = ', &
               trim(tct%gas_name), trim(tct%date), tct%val, tct%trend
        end if
      end do
      write (6,'(a,a)')     ' trg_hist_file           = ', trim(trg_hist_file)
      write (6,'(a,a)')     ' trg_hist_inidate        = ', trim(trg_hist_inidate)
      write (6,'( )')
      write (6,'(a)')       ' ! Clouds in RTTOV'
      write (6,'(a,l1)')    ' use_reff                =' , use_reff
      write (6,'(a,l1)')    ' scale_qi                =' , scale_qi
      write (6,'(a,f9.2)')  ' min_dia_qc              =' , min_dia_qc
      write (6,'(a,f9.2)')  ' min_dia_qi              =' , min_dia_qi
      write (6,'( )')
      write (6,'(a)')       ' ! Minimization'
      write (6,'(a,l1)')    ' no_ps_dep               = ', no_ps_dep
!      write (6,'(a,i6)')    ' no_t_dep                =' , no_t_dep
      write (6,'( )')
      write (6,'(a)')       ' ! Feedback file options'
      write (6,'(a,i6)')    ' monitor_prof            =' , monitor_prof
      write (6,'(a,i6)')    ' mon_ana_prof            =' , mon_ana_prof
      write (6,'(a,i6)')    ' monitor_prof_split      =' , monitor_prof_split
      write (6,'(a,i6)')    ' mon_ana_prof_split      =' , mon_ana_prof_split
      write (6,'(a,30(i4))')' monitor_prof_sat        =' , pack(monitor_prof_sat,  mask=monitor_prof_sat  >=0)
      write (6,'(a,30(i4))')' mon_ana_prof_sat        =' , pack(mon_ana_prof_sat,  mask=mon_ana_prof_sat  >=0)
      write (6,'(a,30(i4))')' monitor_prof_grid       =' , pack(monitor_prof_grid, mask=monitor_prof_grid >=0)
      write (6,'(a,30(i4))')' mon_ana_prof_grid       =' , pack(mon_ana_prof_grid, mask=mon_ana_prof_grid >=0)
      write (6,'(a,f9.2)')  ' monitor_prof_wrbxf      =' , monitor_prof_wrbxf
      write (6,'(a,f9.2)')  ' mon_ana_prof_wrbxf      =' , mon_ana_prof_wrbxf
      write (6,'(a,l1)')    ' lBii                    = ', lBii
      write (6,'(a,i1)')    ' hd_id_vers              = ', hd_id_vers
      write (6,'( )')
      write (6,'(a)')       ' ! Pressure level estimation'
      write (6,'(a,f9.2)')  ' pl_e_bg_t               = ', pl_e_bg_t
      write (6,'(a,f9.2)')  ' pl_e_bg_rh              = ', pl_e_bg_rh
      write (6,'(a,i1)')    ' pl_method               = ', pl_method
      write (6,'(a,l1)')    ' pl_log                  = ', pl_log
      write (6,'(a,f9.2)')  ' pl_wd_thresh            = ', pl_wd_thresh
      write (6,'(a,f9.2)')  ' pl_vis                  = ', pl_vis
      write (6,'(a,f9.2)')  ' plev_min                = ', plev_min
      write (6,'( )')
      write (6,'(a)')       ' ! Diagnostics'
      if (any(spt_debug > 0)) &
           write (6,'(a,10(1x,i9))')    ' spt_debug(obsolete)   = ', spt_debug
      if (any(spt_hd_debug > 0)) &
           write (6,'(a,10(1x,i9))')    ' spt_hd_debug(obsolete)= ', spt_hd_debug
      write (6,'(a,i1)')    ' usd(obsolete)           = ', usd
      write (6,'(a,i1)')    ' spt_debug_prof          = ', spt_debug_prof
      write (6,'(a,l1)')    ' trace_mem_time          = ', trace_mem_time
      write (6,'(a,i5)')    ' debug_opdep_chan        = ', debug_opdep_chan
      ! write (6,'(a,l1)')    ' ltest                   = ', ltest
      write (6,'( )')
      write (6,'(a)')       ' ! Obsolete'
      write (6,'(a,i6)')    ' max_scan                = ', max_scan
      write (6,'(a,i6)')    ' flg_prc                 =' , flg_prc
      write (6,'(a,i6)')    ' flg_prc_amsua           =' , flg_prc_amsua
      write (6,'(a,i6)')    ' flg_prc_amsub           =' , flg_prc_amsub
      write (6,'(a,i6)')    ' flg_prc_hirs            =' , flg_prc_hirs
      write (6,'( )')

      param_file_amsub = param_fil_amsub

    endif
    !--------------------------
    ! 1.c broadcast namelist
    !--------------------------
    ! Input
    call p_bcast (netcdf_path,          dace% pio)
    call p_bcast (data_path,            dace% pio)
    call p_bcast (feedbk_files,         dace% pio)
    call p_bcast (read1pe,              dace% pio)
    call p_bcast (flg_sur,              dace% pio)
    call p_bcast (use_ch_blcklst,       dace% pio)
    call p_bcast (req_nmlst,            dace% pio)
    call p_bcast (filter_ch_stat,       dace% pio)
    call p_bcast (n_filter_ch,          dace% pio)
    ! Cloud check
    call p_bcast (param_file_amsub,     dace% pio)
    call p_bcast (clchk_amsub,          dace% pio)
    call p_bcast (amsub_ec_bnd,         dace% pio)
    call p_bcast (amsub_delta_ch20_ch18,dace% pio)
    ! Cloud analysis
    call p_bcast (cld_top_def,          dace% pio)
    call p_bcast (cld_frc_def,          dace% pio)
    call p_bcast (cld_top_e,            dace% pio)
    call p_bcast (cld_frc_e,            dace% pio)
    call p_bcast (bound_constr,         dace% pio)
    call p_bcast (cld_top_def_v,        dace% pio)
    call p_bcast (cld_frc_def_v,        dace% pio)
    call p_bcast (cld_top_e_v,          dace% pio)
    call p_bcast (cld_frc_e_v,          dace% pio)
    ! Surface_classification
    call p_bcast (surf_class_vers,      dace% pio)
    call p_bcast (snw_frc_mode,         dace% pio)
    call p_bcast (snw_frc_par,          dace% pio)
    call p_bcast (dynret_w_pref,        dace% pio)
    call p_bcast (dynret_avg,           dace% pio)
    call p_bcast (check_recalc_stype,   dace% pio)
    ! Error specification
    call p_bcast (e_bg_ts_sea,          dace% pio)
    call p_bcast (e_bg_ts_land,         dace% pio)
    call p_bcast (e_bg_drts_land,       dace% pio)
    call p_bcast (e_bg_ts_ice,          dace% pio)
    call p_bcast (e_bg_t_top,           dace% pio)
    call p_bcast (e_bg_lnq_top,         dace% pio)
    ! Humidity options
    call p_bcast (lhum_dum_ana,         dace% pio)
    call p_bcast (fg_prof_top,          dace% pio)
    call p_bcast (use_hum_top,          dace% pio)
    call p_bcast (warn_gen_hum,         dace% pio)
    ! Blending IFS
    call p_bcast (p_top,                dace% pio)
    call p_bcast (p_top_ext,            dace% pio)
    call p_bcast (p_blend_ext,          dace% pio)
    call p_bcast (nlev_ext,             dace% pio)
    ! Channel height (l2c) checks
    call p_bcast (l2c_max_tr2ml,        dace% pio)
    call p_bcast (l2c_max_tr2ml_width,  dace% pio)
    call p_bcast (l2c_max_ml2pl,        dace% pio)
    call p_bcast (l2c_max_ml2pl_width,  dace% pio)
    ! t_rad_gopt and t_rad_iopt default values (see TOVS_OBS_CHAN namelists)
    call p_bcast (flag_instr,           dace% pio)
    call p_bcast (n_max_calc_k_def,     dace% pio)
    call p_bcast (n_max_prof_k_def,     dace% pio)
    call p_bcast (n_max_calc_y_def,     dace% pio)
    call p_bcast (n_max_prof_y_def,     dace% pio)
    call p_bcast (quality_mode,         dace% pio)
    call p_bcast (l2c_type,             dace% pio)
    call p_bcast (l2c_god_corr,         dace% pio)
    call p_bcast (l2c_rel_lim,          dace% pio)
    call p_bcast (l2c_abs_lim,          dace% pio)
    call p_bcast (l2c_use_rad,          dace% pio)
    call p_bcast (l2c_max,              dace% pio)
    call p_bcast (l2c_max_trop,         dace% pio)
    call p_bcast (l2c_max_midlat,       dace% pio)
    call p_bcast (l2c_max_polar,        dace% pio)
    call p_bcast (surf_infl_mode,       dace% pio)
    call p_bcast (max_surf_infl,        dace% pio)
    call p_bcast (d_stemp,              dace% pio)
    call p_bcast (d_emiss,              dace% pio)
    call p_bcast (l_max_surf_hgt,       dace% pio)
    ! RTTOV options
    call p_bcast (rttov_version,        dace% pio)
    call p_bcast (rttov_levels,         dace% pio)
    call p_bcast (rt13_clip_gas_opdep,  dace% pio)
    call p_bcast (rt_levs_exact,        dace% pio)
    call p_bcast (rttov_mult_prof,      dace% pio)
    call p_bcast (rt_ir_emis_mod,       dace% pio)
    call p_bcast (rt_mw_emis_mod,       dace% pio)
    call p_bcast (force_fastem,         dace% pio)
    call p_bcast (rt_humi,              dace% pio)
    call p_bcast (rt_use_q2m,           dace% pio)
    call p_bcast (rt_do_lambertian,     dace% pio)
    call p_bcast (rt_salinity,          dace% pio)
    call p_bcast (rt_min_od,            dace% pio)
    call p_bcast (rt_hard_limits,       dace% pio)
    call p_bcast (app_reg_lims,         dace% pio)
    call p_bcast (fix_hgpl,             dace% pio)
    call p_bcast (god_par_file,         dace% pio)
    call p_bcast (wr_god,               dace% pio)
    call p_bcast (atlas_single_inst,    dace% pio)
    call p_bcast (rt_alloc_mode,        dace% pio)
    ! Checks on RTTOV applicability
    call p_bcast (chk_reg_lims,         dace% pio)
    call p_bcast (chk_opdep,            dace% pio)
    call p_bcast (opdep_thresh,         dace% pio)
    call p_bcast (opdep_min_transm,     dace% pio)
    call p_bcast (chk_ps,               dace% pio)
    call p_bcast (chk_god,              dace% pio)
    call p_bcast (god_thresh,           dace% pio)
    call p_bcast (wr_opdep,             dace% pio)
    call p_bcast (qflag_bitmask,        dace% pio)
    ! Tracegases
    call p_bcast (trg_file,             dace% pio)
    call p_bcast (trg_clim_file,        dace% pio)
    do j = 1, size(trg_clim_trend)
      call p_bcast (trg_clim_trend(j),  dace% pio)
    end do
    call p_bcast (trg_hist_file,        dace% pio)
    call p_bcast (trg_hist_inidate,     dace% pio)
    ! Clouds in RTTOV
    call p_bcast (use_reff,             dace% pio)
    call p_bcast (scale_qi,             dace% pio)
    call p_bcast (min_dia_qc,           dace% pio)
    call p_bcast (min_dia_qi,           dace% pio)
    ! Minimization
    call p_bcast (no_ps_dep,            dace% pio)
!    call p_bcast (no_t_dep,             dace% pio)
    ! Feedback file options
    call p_bcast (monitor_prof,         dace% pio)
    call p_bcast (mon_ana_prof,         dace% pio)
    call p_bcast (monitor_prof_split,   dace% pio)
    call p_bcast (mon_ana_prof_split,   dace% pio)
    call p_bcast (monitor_prof_sat,     dace% pio)
    call p_bcast (mon_ana_prof_sat,     dace% pio)
    call p_bcast (monitor_prof_grid,    dace% pio)
    call p_bcast (mon_ana_prof_grid,    dace% pio)
    call p_bcast (monitor_prof_wrbxf,   dace% pio)
    call p_bcast (mon_ana_prof_wrbxf,   dace% pio)
    call p_bcast (lBii,                 dace% pio)
    call p_bcast (hd_id_vers,           dace% pio)
    ! Pressure level estimation
    call p_bcast (pl_e_bg_t,            dace% pio)
    call p_bcast (pl_e_bg_rh,           dace% pio)
    call p_bcast (pl_method,            dace% pio)
    call p_bcast (pl_log,               dace% pio)
    call p_bcast (pl_wd_thresh,         dace% pio)
    call p_bcast (pl_vis,               dace% pio)
    call p_bcast (plev_min,             dace% pio)
    ! Diagnostics
    call p_bcast (spt_debug,            dace% pio)
    call p_bcast (spt_hd_debug,         dace% pio)
    call p_bcast (usd,                  dace% pio)
    call p_bcast (spt_debug_prof,       dace% pio)
    call p_bcast (trace_mem_time,       dace% pio)
    call p_bcast (debug_opdep_chan,     dace% pio)
    ! call p_bcast (ltest,                dace% pio)
    ! Obsolete
    call p_bcast (max_scan,             dace% pio)
    call p_bcast (flg_prc,              dace% pio)
    call p_bcast (flg_prc_amsua,        dace% pio)
    call p_bcast (flg_prc_amsub,        dace% pio)
    call p_bcast (flg_prc_hirs,         dace% pio)

    !-------------------------------------------
    ! set defaults for t_rad_gopt and t_rad_iopt
    !-------------------------------------------
    quality_mode_def   = quality_mode
    l2c_type_def       = l2c_type
    l2c_rel_lim_def    = l2c_rel_lim
    l2c_abs_lim_def    = l2c_abs_lim
    l2c_use_rad_def    = l2c_use_rad
    l2c_max_def        = l2c_max
    l2c_max_trop_def   = l2c_max_trop
    l2c_max_midlat_def = l2c_max_midlat
    l2c_max_polar_def  = l2c_max_polar
    surf_infl_mode_def = surf_infl_mode
    max_surf_infl_def  = max_surf_infl
    d_stemp_def        = d_stemp
    d_emiss_def        = d_emiss
    l_max_surf_hgt_def = l_max_surf_hgt
    ! prepare defaults for l2c check
    if (l2c_max_trop_def   < 0._wp) l2c_max_trop_def   = l2c_max_def
    if (l2c_max_midlat_def < 0._wp) l2c_max_midlat_def = l2c_max_def
    if (l2c_max_polar_def  < 0._wp) l2c_max_polar_def  = l2c_max_def

    call set_rttov_vers (rttov_version, rttov_levels)
    if (rt_alloc_mode >= 2) rtifc_alloc_mode = 1

    trg_file      = path_file (input, trg_file     )
    trg_clim_file = path_file (input, trg_clim_file)

    usd_rad = usd
    ldeb_spot = any(spt_hd_debug > 0 .or. spt_debug > 0)

    first  = .true.
    i_read = 0
    i_nml  = 1
    if (dace% lpio) then
      loop_nml: do
        call position_nml ('TOVS_OBS_CHAN' ,lrewind=first ,status=ierr_pos)
        if (ierr_pos == POSITIONED) then
          i_read = i_read + 1
          call read_tovs_obs_chan_nml(nnml, rad_set(i_nml), status=ierr_read, i_read=i_read)
          where(rad_set(i_nml)%iopts(:)%rt_mw_emis_mod < 0) rad_set(i_nml)%iopts(:)%rt_mw_emis_mod = rt_mw_emis_mod
          if (ierr_read == 0) then
            if (first) print*,repeat(' ',79)
            write(header, '("Read TOVS_OBS_CHAN namelist number ",I2," (",I2,")")') &
                 &i_nml, i_read
            print*
            call print_rad_set(rad_set(i_nml), header=header)
            i_nml = i_nml + 1
            if (i_nml > size(rad_set)) call  finish('read_tovs_nml', &
                 &'Too many TOVS_OBS_CHAN namelists')
          else
            call finish('read_tovs_nml', 'Failed to read TOVS_OBS_CHAN namelist')
          end if
        else
          exit loop_nml
        end if
        first=.false.
      end do loop_nml

      ! Read emissivity namelists
      call read_nml_emis

      ! Read tskin namelists
      call read_nml_tskin

    endif
    n_set = i_nml - 1
    call p_bcast(n_set, dace% pio)
    do i=1, n_set
      call p_bcast(rad_set(i), dace% pio)
    end do

    ! Analyze trace gas options
    glob_use_o3  = 0
    glob_use_co2 = 0
    do i = 1, n_set
      do j = 1, rad_set(i)%n_instr
        glob_use_o3  = ior(glob_use_o3 , rad_set(i)%iopts(j)%use_o3 )
        glob_use_co2 = ior(glob_use_co2, rad_set(i)%iopts(j)%use_co2)
      end do
    end do

  end subroutine read_tovs_nml
!==============================================================================
  subroutine scan_satpp_feedback ()
    !------------------------------------------------------------
    ! 4. read ATOVS data from feedback file, setup file inventory
    !------------------------------------------------------------
    type (t_radv) :: rad             ! radiances data type
    integer       :: i, ierr, n_rec, n_chan
    logical       :: valid
    integer       :: n_valid        ! number of valid input files

    n_valid = 0
    do i = 1, mf
      if (feedbk_files(i)/='') then
        n_rec  = 0
        n_chan = 0
        if (dace% lpio) then
          call read_satpp_feedbk (path_file (netcdf_path,feedbk_files(i)), &
               rad, valid=valid, lread=.false., status=ierr)
          !----------------------------------
          ! write warning if file not present
          !----------------------------------
          if (ierr/=0) then
            write (6,'(a)') repeat('*',79)
            write (0,'(a)') repeat('*',79)
            write (6,*) 'WARNING, cannot open: '//trim(feedbk_files(i))
            write (6,*) trim(err_msg)
            write (0,*) 'WARNING, cannot open: '//trim(feedbk_files(i))
            write (0,*) trim(err_msg)
            write (6,'(a)') repeat('*',79)
            write (0,'(a)') repeat('*',79)
          endif
          if (.not.valid) then
            write (6,'(a)') repeat('*',79)
            write (0,'(a)') repeat('*',79)
            write (6,*) 'WARNING, invalid satpp-file :',trim(feedbk_files(i))
            write (6,*) trim(err_msg)
            write (0,*) 'WARNING, invalid satpp-file :',trim(feedbk_files(i))
            write (6,*) trim(err_msg)
            write (6,'(a)') repeat('*',79)
            write (0,'(a)') repeat('*',79)
            ierr = 1
          endif
          if (ierr == 0) then
            n_rec  = rad%   n_rec
            n_chan = rad%i% n_chan
          endif
          call destruct(rad)
        endif
        !---------------------
        ! store file inventory
        !---------------------
        call p_bcast (format_fdbk, dace% pio)
        call p_bcast (ierr,        dace% pio)
        call p_bcast (n_rec,       dace% pio)
        call p_bcast (ierr,        dace% pio)
        if (ierr==0) then
          n_valid = n_valid + 1
          call add_source (netcdf_path, feedbk_files(i), &
                           filetype    = FT_SATPP,       &
                           obstype     = OT_RAD,         &
                           entries     = n_rec,          &
                           nobs        = n_rec * n_chan, &
                           file_format = format_fdbk(i)  )
          if (dace% lpio) then
            if(n_source>1) bufr_inv% subseto (n_source) = bufr_inv% subseto (n_source-1)
            bufr_inv   (OT_RAD)% file    (n_source)   = .true.
            bufr_inv   (OT_RAD)% nrec                 = bufr_inv (OT_RAD)% nrec + n_rec
            bufr_inv   (OT_RAD)% nsubset              = bufr_inv (OT_RAD)% nsubset  + n_rec
            bufr_inv   (OT_RAD)% subseto(n_source)    = bufr_inv (OT_RAD)% subseto(n_source) + n_rec
          endif
          call p_bcast (bufr_inv, dace% pio)
        endif
      endif
    end do
    !----------------------------------------------------
    ! check if merging with IFS stratosphere is necessary
    ! (but do not force for ICON)
    !----------------------------------------------------
    if (interp_strato == 0 .and. n_valid > 0 .and. model /= "ICON") then
      interp_strato = 1
      if (dace% lpio) then
        write (6,*)
        write (6,*) 'scan_satpp_feedback: interp_strato forced to 1'
        write (6,*)
      endif
    endif

  end subroutine scan_satpp_feedback
!==============================================================================
  subroutine decr_tovs_use(sp, obs, check, set, tovs, state, bit, not_bit, instr, instr_ind, &
       chan_mask, flag_rpt, hint_debug, tovs_flag, tovs_flag_val)
    type(t_spot),     intent(inout)                  :: sp           ! spot to be flagged
    type(t_obs),      intent(inout)                  :: obs          ! obs to be flagged
    integer,          intent(in)                     :: check        ! check
    type(t_rad_set),  intent(in),   optional, target :: set          ! rad_set
    type(t_tovs),     intent(in),   optional, target :: tovs         ! t_tovs
    integer,          intent(in),   optional         :: state        ! state
    integer,          intent(in),   optional         :: bit          ! flag only, if bit is set in set%flag
    integer,          intent(in),   optional         :: not_bit      ! flag only, if not_bit is not set in set%flag
    integer,          intent(in),   optional         :: instr        ! flag only instrument "instr"
    integer,          intent(in),   optional         :: instr_ind    ! flag only instrument "instr(instr_ind)"
    logical,          intent(in),   optional         :: chan_mask(:) ! flag only channels with chan_mask=.true.
    logical,          intent(in),   optional         :: flag_rpt     ! whether to flag the report
    character(len=*), intent(in),   optional         :: hint_debug   ! hint for debug output
    integer,          intent(in),   optional         :: tovs_flag    ! flag to be set in t_tovs%flag(:)
    integer,          intent(in),   optional         :: tovs_flag_val! value "added" to tovs_flag (with ior)

    character(len=300)           :: prefix
    type(t_tovs),    target      :: tovs_loc
    type(t_tovs),    pointer     :: tt => null()
    type(t_rad_set), pointer     :: rs => null()
    type(t_tovs_instr)           :: ti
    integer                      :: new_state
    integer                      :: i, n, n_chan, l2c_bit, tovs_io
    logical                      :: flag_rpt_aux
    logical                      :: l_tovs, l_l2c, l_mask, l_debug, l_tf
    logical,         allocatable :: mask(:)

    real(wp)                     :: l2c_max

    l_tovs = .false.
    l_l2c  = (check == CHK_CLOUD .or. check == CHK_SURF)
    l_tf   = present(tovs_flag) .or. present(tovs_flag_val)

    if (present(instr) .and. present(instr_ind)) &
         call finish('decr_tovs_use', 'instrument ID and instrument index not supported.')

    if (present(set)) then
     if (associated(set%flag)) then
      l_l2c = l_l2c .and. (any(btest(set%flag(1:set%n_chan), USE_L2C_SURF)) .or. &
                           any(btest(set%flag(1:set%n_chan), USE_L2C_CLD )))
     end if
    end if
    if (present(tovs)) then
      do i = 1, n_set
        rs => rad_set(i)
        if (tovs%id_rs == rs%id .and. associated(rs%flag)) then
          l_l2c = l_l2c .and. (any(btest(rs%flag(1:rs%n_chan), USE_L2C_SURF)) .or. &
                               any(btest(rs%flag(1:rs%n_chan), USE_L2C_CLD )))
          exit
        end if
      end do
    end if

    l_mask = present(bit)       .or. present(not_bit) .or. present(instr) .or. &
             present(chan_mask) .or. l_l2c

    l_debug = ldeb(sp)

    tt => null()
    if (present(set) .and. .not. (l_l2c .or. l_mask .or. l_debug .or. l_tf)) then
      ! No t_tovs required
      rs => set
    else
      tovs_io = 0
      if (l_mask .or. l_debug) tovs_io = tovs_io + TTOVS_CI
      if (l_tf               ) tovs_io = tovs_io + TTOVS_FLAG
      if (l_l2c              ) tovs_io = tovs_io + TTOVS_L2C
      if (present(tovs)) then
        if (iand(tovs%init, tovs_io) == tovs_io) then
          tt => tovs
          call get_tovs_rs(tt, rs=rs, ti=ti)
        end if
      end if
      if (.not.associated(tt)) then
        if (sp%p%n > 0) then
          call load(obs, sp, tovs_loc, rs=rs, ti=ti, tovs_io=tovs_io)
          l_tovs = .true.
          tt => tovs_loc
        else
          call finish('decr_tovs_use','t_tovs not initialized')
        end if
      end if
    end if
    if (l_mask .and. associated(rs%flag)) then
      l_l2c = l_l2c                                                      .and. &
              (any(btest(rs%flag(tt%ci(1:tt%nchan)), USE_L2C_SURF)).or.        &
               any(btest(rs%flag(tt%ci(1:tt%nchan)), USE_L2C_CLD ))    )
      ! for backwards compatibility
      if (present(not_bit)) then
        if (not_bit==USE_MINTSURF) then
          if (.not.any(btest(rs%flag(tt%ci(1:tt%nchan)), USE_MINTSURF) .or. &
                       btest(rs%flag(tt%ci(1:tt%nchan)), USE_L2C_SURF))) l_mask = .false.
        end if
      end if
    end if

    flag_rpt_aux = .true.
    if (present(instr))    flag_rpt_aux = (rs%grid == instr)
    if (present(flag_rpt)) flag_rpt_aux = flag_rpt

    if (l_debug) then
      write(usd,*) dpref,'decr_tovs_use check',check
      if (present(hint_debug))    write(usd,*) dpref,' '//trim(hint_debug)
      if (present(set))           write(usd,*) dpref,' set',set%id
      if (present(state))         write(usd,*) dpref,' state',state
      if (present(bit))           write(usd,*) dpref,' bit',bit
      if (present(not_bit))       write(usd,*) dpref,' not_bit',not_bit
      if (present(instr))         write(usd,*) dpref,' instr',instr
      if (present(instr_ind))     write(usd,*) dpref,' instr_ind',instr_ind
      if (present(chan_mask))     write(usd,*) dpref,' chan_mask',chan_mask
      if (present(flag_rpt))      write(usd,*) dpref,' flag_rpt',flag_rpt
      if (present(tovs_flag))     write(usd,*) dpref,' tovs_flag',tovs_flag
      if (present(tovs_flag_val)) write(usd,*) dpref,' tovs_flag_val',tovs_flag_val
      prefix = 'decr_tovs_use'
      if (present(hint_debug)) prefix=trim(prefix)//' '//trim(hint_debug)
      if (associated(rs%flag)) then
        call debug_spot(sp=sp, obs=obs, flags=rs%flag(tt%ci(1:tt%nchan)), hint=trim(prefix)//' start')
      else
        call debug_spot(sp=sp, obs=obs, hint=trim(prefix)//' start')
      end if
    end if

    if (present(state)) then
      new_state = state
    else
      new_state = rept_use(sp% hd% obstype)% use (check)
    end if

    if (l_mask) then
      n_chan = tt%nchan
      allocate(mask(n_chan))

      if (present(instr)) then
         mask(:) = .false.
         do i = 1, ti%n_instr
            if (rs%instr(ti%ii(i)) == instr) then
               if (present(chan_mask)) then
                  n = min(size(chan_mask), ti%n_ch_i(i))
                  mask(ti%o_ch_i(i)+1:ti%o_ch_i(i)+n) = chan_mask(1:n)
               else
                  mask(ti%o_ch_i(i)+1:ti%o_ch_i(i)+ti%n_ch_i(i)) = .true.
               end if
               exit
            end if
         end do
      else if (present(instr_ind)) then
         mask(:) = .false.
         if (present(chan_mask)) then
            n = min(size(chan_mask), ti%n_ch_i(instr_ind))
            mask(ti%o_ch_i(instr_ind)+1:ti%o_ch_i(instr_ind)+n) = chan_mask(1:n)
         else
            mask(ti%o_ch_i(instr_ind)+1:ti%o_ch_i(instr_ind)+ti%n_ch_i(instr_ind)) = .true.
         end if
      else if (present(chan_mask)) then
        n = min(size(chan_mask), n_chan)
        mask(1:n) = chan_mask(1:n)
        mask(n+1:n_chan) = .false.
      else
        mask(:) = .true.
      end if

      if (associated(rs%flag)) then
        if (present(bit)) then
          mask(1:n_chan) = mask(1:n_chan) .and. btest(rs%flag(tt%ci(1:n_chan)), bit)
        end if
        if (present(not_bit)) then
          mask(1:n_chan) = mask(1:n_chan) .and. .not.btest(rs%flag(tt%ci(1:n_chan)), not_bit)
        end if
      end if
      if (l_l2c) then
        if (rs%gopts%lev_mode > 0) call finish('decr_tovs_use', 'l2c_max not valid with lev_mode > 0')
        l2c_max = get_l2c_max(sp% col% c% dlat, rs%iopts(1)%l2c_max_trop, &
                              rs%iopts(1)%l2c_max_midlat, rs%iopts(1)%l2c_max_polar)
        if (check == CHK_CLOUD) then
          l2c_bit = USE_L2C_CLD
        else
          l2c_bit = USE_L2C_SURF
        end if
        if (rs% gopts% quality_mode == 7) then
          ! Check whether this spot is only surviving because of the l2c rules
          if (new_state <= STAT_REJECTED .and. sp%pcc > 0) then
            if ( all(obs% body(sp%o%i+1:sp%o%i+n_chan)% use% state <= STAT_REJECTED .or. &
                 mask(1:n_chan)) ) sp%pcc = 0
          end if
          if (l_debug) write(usd,*) dpref,'qm=7 pcc=',sp%pcc,check
        end if
        if (associated(rs%flag)) then
          where (btest(rs%flag(tt%ci(1:n_chan)), l2c_bit))
            mask(1:n_chan) = mask(1:n_chan) .and. (tt% l2c(1:n_chan) > l2c_max)
          end where
        end if
      end if

      if (l_debug) write(usd,*) dpref,'mask=',mask

      do i = 1, n_chan
        if (mask(i)) then
          call decr_use(obs% body(sp%o%i+i)% use, check=check, state=new_state, lflag=.true.)
          if (l_tf) then
            if (present(tovs_flag    )) tt%flag(i) = ibset(tt%flag(i), tovs_flag    )
            if (present(tovs_flag_val)) tt%flag(i) = ior  (tt%flag(i), tovs_flag_val)
          end if
        end if
      end do

!      flag_rpt_aux = flag_rpt_aux .and. any(mask(:))
    else
      call decr_use(obs% body(sp%o%i+1:sp%o%i+sp%o%n)% use, check=check, lflag=.true.)
      if (l_tf) then
        if (present(tovs_flag    )) tt%flag(i) = ibset(tt%flag(i), tovs_flag    )
        if (present(tovs_flag_val)) tt%flag(i) = ior  (tt%flag(i), tovs_flag_val)
      end if
    end if

    if (flag_rpt_aux) then
      ! Set flag, but do not decrease status
      call decr_use(sp% use, state=int(sp% use% state), check=check, lflag=.true.)
    end if

    if (l_debug) then
      if (associated(tt) .and. associated(rs%flag)) then
        call debug_spot(sp=sp, obs=obs, flags=rs%flag(tt%ci(1:tt%nchan)),hint=trim(prefix)//' end')
      else
        call debug_spot(sp=sp, obs=obs, hint=trim(prefix)//' end')
      end if
    end if

    if (l_tf) then
      call store(obs, sp, tt, tovs_io=TTOVS_FLAG)
    end if
    if (l_tovs) call destruct(tovs_loc)

  end subroutine decr_tovs_use

!------------------------------------------------------------------------------

  subroutine write_rttov_prof (obs, bg, ana, e_bg, lana, Pb_aprx, &
                                zi, an_bg_ps, dJo, time, HBHR, Roo, rhs )
  type(t_obs_set)  ,intent(in)             :: obs     ! observations
  type(t_vector)   ,intent(in)             :: bg      ! background
  type(t_vector)   ,intent(in)             :: ana     ! analysis
  type(t_vector)   ,intent(in)             :: e_bg    ! background error
  logical          ,intent(in)             :: lana    ! .true.for analysis scan
  integer          ,intent(in)             :: Pb_aprx ! approximation of Pb
  type(t_vector)   ,intent(in)   ,optional :: zi      ! 'z' in intp.space
  type(t_vector)   ,intent(inout),optional :: an_bg_ps! analysis in PSAS space
  type(t_vector)   ,intent(in)   ,optional :: dJo     ! gradient of J_obs
  type(t_time)     ,intent(in)   ,optional :: time    ! analysis time
  type(t_matrix)   ,intent(in)   ,optional :: HBHR    ! HBH+R
  type(t_matrix)   ,intent(in)   ,optional :: Roo     ! R
  type(t_vector)   ,intent(in)   ,optional :: rhs     ! a-fg=BHt(HBHt+R)^-1*rhs

    type t_profile
      real(wp), pointer :: t_fg(:,:)      => null()
      real(wp), pointer :: q_fg(:,:)      => null()
      real(wp), pointer :: ts_fg(:)       => null()
      real(wp), pointer :: ps_fg(:)       => null()
      real(wp), pointer :: cld_top(:)     => null()
      real(wp), pointer :: cld_frc(:,:)   => null()
      real(wp), pointer :: snw_frc(:)     => null()
    end type t_profile

    character(len=20)          :: proc = 'write_rttov_prof'
    character(len=300)         :: fname
    character(len=12)          :: rtime, atime
    type(t_obs) ,      pointer :: o    => null()
    type(t_spot),      pointer :: si   => null()
    type(t_radv),      pointer :: r(:) => null()
    type(t_profile),   pointer :: p(:) => null()
    type(t_matrix)             :: Pbi
    type(t_vector)             :: ana_
    type(t_vector)             :: afg_psas
    type(t_nc_attr)            :: attrs(6)
    integer                    :: m_prof
    integer                    :: ib, is, i, j, i1, in   ! indices
    integer                    :: status                 ! exit status
    integer                    :: n_prof                 ! number of profiles
    integer                    :: n_profb (size(obs%o))  ! no. profiles in box
    integer                    :: m_chan                 ! max. no. channels
    integer                    :: n_par
    logical                    :: lsort, lH, lK, lB, lpsas, lHBH
    ! split into separate files
    integer                    :: m_prof_split
    integer                    :: m_prof_sat(m_rad_set)
    integer                    :: m_prof_grid(m_rad_set)
    logical                    :: l_sat, l_grid
    integer                    :: nf                     ! number of files
    integer                    :: if_set(m_rad_set)      ! file for rad_set
    integer                    :: i_aux(m_rad_set)       ! auxiliary arry
    character(len=300)         :: fnames(m_rad_set)
    integer                    :: nset_f                 ! max. number of sets in file
    integer                    :: i_f, iset
    logical                    :: l_partial
    ! sqauential writing of boxes
    real(wp)                   :: m_prof_wrbxf
    integer                    :: nbox, mbox
    integer                    :: ioff, ioff_new         ! offset
    integer                    :: ibox_s, ibox_e


    if (lana) then
      m_prof       = mon_ana_prof
      m_prof_split = mon_ana_prof_split
      m_prof_sat   = mon_ana_prof_sat
      m_prof_grid  = mon_ana_prof_grid
      m_prof_wrbxf = mon_ana_prof_wrbxf
    else
      m_prof       = monitor_prof
      m_prof_split = monitor_prof_split
      m_prof_sat   = monitor_prof_sat
      m_prof_grid  = monitor_prof_grid
      m_prof_wrbxf = monitor_prof_wrbxf
    endif

    if (n_set <= 0) return
    if (iand(m_prof, MP_NETCDF) == 0) return
    if (iand(m_prof, MP_OLD)    /= 0) call finish('write_rttov_prof','old version of &
         &write_rttov_prof not supported anymore.')

    FTRACE_BEGIN('write_rttov_prof')

    lsort = (iand(m_prof, MP_NOSORT ) == 0)
    lH    = (iand(m_prof, MP_H      ) /= 0)
    lK    = lH .and. present(HBHR) .and. present(zi) .and. present(rhs)
    lpsas = (iand(m_prof,MP_PSAS_ANA) /= 0 .and. present(zi) .and. present(time))
    lB    = (iand(m_prof, MP_B      ) /= 0)
    lHBH  = (iand(m_prof, MP_HBH    ) /= 0)


    !-----------------------
    ! Prepare analysis stuff
    !-----------------------
    if (lana) then
      call construct   (ana_,    bg% info, 'ana' )
      ana_ = ana
      if (present(zi)) then
        if (lpsas) then
          call construct (afg_psas, bg% info, 'a-fg_psas')
          call construct (Pbi, bg% info, 'Pbi')
          call set_Pb (Pbi, obs, rearth, time, i_ensb=0, HHt=.false.)
          afg_psas = Pbi * zi
          call destruct (Pbi)
        end if
        do ib = 1,size(obs% o)
          o => obs% o(ib)
          if (dace% pe == o% pe) then
            do is=1,o% n_spot
              si => o% spot(is)
              if (si% hd% obstype /= OT_RAD)       cycle
              if (si% use% state  <= STAT_DISMISS) cycle
              i1 = si%i%i+1
              in = si%i%i+si%i%n
              do i = i1, in
                if (o% t_int (i) == OBS_DUM) then
                  ana_% s(ib)% x(i) = ana_% s(ib)% x(i) + zi%   s(ib)% x(i)    &
                       * e_bg% s(ib)% x(i) ** 2
                endif
              end do
              if (lpsas) afg_psas% s(ib)% x(i1:in) = afg_psas% s(ib)% x(i1:in) &
                   + bg% s(ib)% x(i1:in)
            end do
          end if
        end do
      end if
    end if
    if (lK) call uncompress_cov

    !-------------------------------------------------------
    ! Determine number of files and satid/grid for each file
    !-------------------------------------------------------
    if_set(1:n_set) = 0
    if (any(m_prof_sat(:) >= 0)) then
      do i = 1, n_set
        if (.not.any(m_prof_sat(:) == rad_set(i)% satid)) if_set(i) = -1
      end do
    end if
    if (any(m_prof_grid(:) >= 0)) then
      do i = 1, n_set
        if (.not.any(m_prof_grid(:) == rad_set(i)% grid)) if_set(i) = -1
      end do
    end if
    l_partial = .true.
    l_sat  = iand(m_prof_split, 1) > 0
    l_grid = iand(m_prof_split, 2) > 0
    if (.not.(l_sat.or.l_grid)) then
      nf = 1
      where(if_set(1:n_set) >= 0) if_set(1:n_set) = nf
      fnames(nf) = ''
      l_partial = .false.
    else if (l_sat .and. l_grid) then
      nf = n_set
      do i = 1, n_set
        if (if_set(i) < 0) cycle
        if_set(i) = i
        write(fnames(i), '(I3.3,I3.3)') rad_set(i)% satid, rad_set(i)% grid
      end do
    else if (l_sat) then
      nf = 0
      set_loop1: do i = 1, n_set
        if (if_set(i) < 0) cycle
        do j = 1, nf
          if (i_aux(j) == rad_set(i)% satid) then
            if_set(i) = j
            cycle set_loop1
          end if
        end do
        nf = nf + 1
        i_aux(nf) = rad_set(i)% satid
        if_set(i) = nf
        write(fnames(nf), '(I3.3)') rad_set(i)% satid
      end do set_loop1
    else if (l_grid) then
      nf = 0
      set_loop2: do i = 1, n_set
        if (if_set(i) < 0) cycle
        do j = 1, nf
          if (i_aux(j) == rad_set(i)% grid) then
            if_set(i) = j
            cycle set_loop2
          end if
        end do
        nf = nf + 1
        i_aux(nf) = rad_set(i)% grid
        if_set(i) = nf
        write(fnames(nf), '(I3.3)') rad_set(i)% grid
      end do set_loop2
    end if

    do i_f = 1, nf
      !---------------------------------------
      ! determine rad_set entries in this file
      !---------------------------------------
      nset_f = 0
      do i = 1, n_set
        if (if_set(i) == i_f) then
          nset_f = nset_f + 1
          i_aux(nset_f) = i
        end if
      end do

      ! Check consistent lev_mode
      if (any(rad_set(i_aux(2:nset_f))%gopts%lev_mode /= rad_set(i_aux(1))%gopts%lev_mode)) then
        call finish('write_rttov_prof', 'inconsistent lev_mode within *RTOVP* file')
      end if

      !------------------------------
      ! count no.profiles on all PE's
      !------------------------------
      n_prof  = 0
      n_profb = 0
      m_chan  = 0
      do ib = 1,size(obs% o)
        o => obs% o(ib)
        if (dace% pe /= o%pe) cycle
        do is = 1, o% n_spot
          si => o% spot(is)
          if (si% hd% obstype /= OT_RAD)       cycle
          if (si% use% state  <= STAT_DISMISS) cycle
          iset = set_indx(rad_set(i_aux(1:nset_f)), satid=int(si%hd%satid), grid=int(si%hd%grid_id))
          if (iset <= 0) cycle
          n_profb(ib) = n_profb(ib) + 1
          m_chan      = max (m_chan, si% o% n)
        end do
      end do
      n_profb = p_sum  (n_profb)
      n_prof  = sum    (n_profb)
      m_chan  = p_max  (m_chan)
      if (n_prof==0) cycle

      allocate(r(nset_f))
      if (lpsas) allocate(p(nset_f))
      if (dace% lpio) then
        do i = 1, nset_f              ! No fill_rad call for "open" call of write_rtovp
          r(i)%i = rad_set(i_aux(i))  ! This allows access to the options within write_rtovp
        end do                        ! nevertheless
        if (lana) then
          fname = path_file(aux,'cofRTOVP')
          n_par = 3
          if (lpsas) n_par = 4
        else
          fname = path_file(aux,'monRTOVP')
          n_par = 2
        endif
        fname = trim(fname)//trim(fnames(i_f))//'.nc'
        atime = cyyyymmddhhmm (ana_time)
        rtime = cyyyymmddhhmm (run_time)
        attrs(1) = t_nc_attr('title', '3dvar ATOVS monitoring file: '//trim(fname), NF90_FILL_INT)
        attrs(2) = t_nc_attr('experiment', '', nex)
        attrs(3) = t_nc_attr('analysis_time', atime, NF90_FILL_INT)
        attrs(4) = t_nc_attr('run_time', rtime, NF90_FILL_INT)
        attrs(5) = t_nc_attr('monitor_prof', '', m_prof)
        attrs(6) = t_nc_attr('b_aprx', '', Pb_aprx)
        call write_rtovp(r, m_prof, mode=WRM_OPEN,  status=status, name=fname, n_prof=n_prof, &
             n_par=n_par, n_chan=m_chan, n_pc=n_pc, lH=lH, lK=lK, lB=lB, lHBH=lHBH, attrs=attrs)
        if (status /= 0) call finish(proc, 'open '//trim(fname))
        ! if (dace% lpio) print*,'write file: '//trim(fname)
        ioff = 0
        call destruct(r(:))
      end if

      mbox = size(obs% o)
      nbox = nint(mbox*m_prof_wrbxf)
      nbox = max(1, min(nbox, mbox))
      if (nbox < mbox .and. lsort) then
        lsort = .false.
        if (dace% lpio) call message('write_rttov_prof', 'warning: profile sorting and &
             &sequential writing of boxes is not implemented as yet -> We do not sort the profiles!!')
      end if

      box_loop: do ibox_s = 1, mbox, nbox
        ibox_e = min(ibox_s + nbox - 1, mbox)
        ! if (dace% lpio) print*,'write boxes',ibox_s,ibox_e

        call fill_rad(r, rad_set(i_aux(1:nset_f)), obs, bg, e_bg=e_bg,          &
          min_stat=STAT_DISMISS, pe_gather=dace% pio, ibox_s=ibox_s, ibox_e=ibox_e,    &
          lmon=.true., lH=lH, lB=lB, lHBH=lHBH, HBHR=HBHR, Roo=Roo, zi=zi, rhs=rhs, lpartial=l_partial)
        if (dace% lpio) then
          call write_rtovp(r, m_prof, mode=WRM_WRITE, status=status, name=fname, n_chan=m_chan, &
               pe=dace% pe, lsort=lsort, lH=lH, lK=lK, lB=lB, lHBH=lHBH, keep_sorting=.true., offset=ioff)
          if (status /= 0) call finish(proc, 'write '//trim(fname))
          ioff_new = ioff + sum(r(:)%n_rec)

          ! Prepare writing of background error and keep fg profile for psas
          do i = 1, nset_f
            if (r(i)% n_rec > 0) then
              if (lpsas) then
                p(i)%t_fg    => r(i)%t_fg
                p(i)%q_fg    => r(i)%q_fg
                p(i)%ts_fg   => r(i)%ts_fg
                p(i)%ps_fg   => r(i)%ps_fg
                p(i)%cld_top => r(i)%cld_top
                p(i)%cld_frc => r(i)%cld_frc
                p(i)%snw_frc => r(i)%snw_frc
              else
                deallocate(r(i)%t_fg)
                deallocate(r(i)%q_fg)
                deallocate(r(i)%ts_fg)
                deallocate(r(i)%tsm_fg)
                deallocate(r(i)%ps_fg)
                if (associated(r(i)%cld_top)) deallocate(r(i)%cld_top)
                if (associated(r(i)%cld_frc)) deallocate(r(i)%cld_frc)
                if (associated(r(i)%snw_frc)) deallocate(r(i)%snw_frc)
              end if
              deallocate(r(i)%v10_abs_fg)
              if (associated(r(i)%emis_pc)) deallocate(r(i)%emis_pc)
              if (associated(r(i)%emiss  )) deallocate(r(i)%emiss  )
              r(i)%t_fg    => r(i)%t_eb       ; r(i)%t_eb       => null()
              r(i)%q_fg    => r(i)%q_eb       ; r(i)%q_eb       => null()
              r(i)%ts_fg   => r(i)%ts_eb      ; r(i)%ts_eb      => null()
              r(i)%ps_fg   => r(i)%ps_eb      ; r(i)%ps_eb      => null()
              r(i)%cld_top => r(i)%cld_top_eb ; r(i)%cld_top_eb => null()
              r(i)%cld_frc => r(i)%cld_frc_eb ; r(i)%cld_frc_eb => null()
              r(i)%snw_frc => r(i)%snw_frc_eb ; r(i)%snw_frc_eb => null()
              if (associated(r(i)%emis_pc_eb)) then
                r(i)%emis_pc => r(i)%emis_pc_eb ; r(i)%emis_pc_eb => null()
              end if
            end if
          end do
          call write_rtovp(r, m_prof, mode=WRM_WRITE, status=status, name=fname, n_chan=m_chan, &
               par='fg_error', pe=dace% pe, lsort=lsort, keep_sorting=.true., offset=ioff)
          if (status /= 0) call finish(proc, 'write '//trim(fname))
          call destruct(r(:))
        end if
        ! Analysis
        if (lana) then
          call fill_rad(r, rad_set(i_aux(1:nset_f)), obs, ana_, min_stat=STAT_DISMISS, &
               pe_gather=dace% pio, ibox_s=ibox_s, ibox_e=ibox_e, lmon=.true., lmeta=.false., lpartial=l_partial)
          if (dace% lpio) then
            call write_rtovp(r, m_prof, mode=WRM_WRITE, status=status, name=fname, n_chan=m_chan, &
                 par='ana', pe=dace% pe, lsort=lsort, keep_sorting=.true., offset=ioff)
            if (status /= 0) call finish(proc, 'write '//trim(fname))
            call destruct(r(:))
          end if

          if (lpsas) then
            call fill_rad(r, rad_set(i_aux(1:nset_f)), obs, afg_psas, min_stat=STAT_DISMISS, &
                 pe_gather=dace% pio, ibox_s=ibox_s, ibox_e=ibox_e, lmon=.true., lmeta=.false., lpartial=l_partial)
            if (dace% lpio) then
              do i = 1, nset_f
                if (r(i)% n_rec > 0) then
                  r(i)%t_fg    = r(i)%t_fg    - p(i)%t_fg     ; deallocate(p(i)%t_fg   )
                  r(i)%q_fg    = r(i)%q_fg    - p(i)%q_fg     ; deallocate(p(i)%q_fg   )
                  r(i)%ts_fg   = r(i)%ts_fg   - p(i)%ts_fg    ; deallocate(p(i)%ts_fg  )
                  r(i)%ps_fg   = r(i)%ps_fg   - p(i)%ps_fg    ; deallocate(p(i)%ps_fg  )
                  r(i)%cld_top = r(i)%cld_top - p(i)%cld_top  ; deallocate(p(i)%cld_top)
                  r(i)%cld_frc = r(i)%cld_frc - p(i)%cld_frc  ; deallocate(p(i)%cld_frc)
                  r(i)%snw_frc = r(i)%snw_frc - p(i)%snw_frc  ; deallocate(p(i)%snw_frc)
                  if (associated(r(i)%emis_pc)) deallocate(r(i)%emis_pc)
                  if (associated(r(i)%emiss  )) deallocate(r(i)%emiss  )
                end if
              end do
              call write_rtovp(r, m_prof, mode=WRM_WRITE, status=status, name=fname, n_chan=m_chan, &
                   par='a-fg_aprox_psas', pe=dace% pe, lsort=lsort, keep_sorting=.true., offset=ioff)
              if (status /= 0) call finish(proc, 'write '//trim(fname))
              call destruct(r(:))
            end if
          end if
        end if
        if (dace% lpio) ioff = ioff_new
      end do box_loop

      if (dace% lpio) then
        call write_rtovp(r, m_prof, mode=WRM_CLOSE, status=status, name=fname)
        if (status /= 0) call finish(proc, 'close '//trim(fname))
      end if
      deallocate(r)

    end do

    if (lana) then
      call destruct(ana_)
      if (lpsas) call destruct(afg_psas)
    end if
    if (lK)   call compress_cov

    FTRACE_END('write_rttov_prof')

  end subroutine write_rttov_prof


  !> Superobbing of TOVS observations
  subroutine thin_superob_tovs(obs, y, ldism)
    type(t_obs_set), intent(inout)          :: obs    ! observation data (3dvar type)
    type(t_vector),  intent(inout)          :: y      ! background observation space
    logical,         intent(out),  optional :: ldism  ! whether obs were dismissed

    type(t_radv), target      :: rad(1)      ! "full" satellite data set
    type(t_radv), pointer     :: r => null() !
    integer                   :: nrec        ! Number of obs. in t_radv

    ! Comment on variable names: All variables named *_pio are only meaningful on dace%pio.
    !                            Usually they contain information to be distributed to other PEs.

    ! Information on superobbing Boxes
    integer,      allocatable :: ibase     (:) ! Full info on so-boxes for ALL obs (from all PEs)
    integer,      allocatable :: ibase_pio (:) ! ibase reorganized for alltoall communication
    integer,      allocatable :: ibase_    (:) ! Info on so-boxes on this PE
    integer,      allocatable :: ns_pio    (:) ! Number of spots in superobbing on all PEs
    integer                   :: ns_(1), ns    ! Number of spots in superobbing on this PE
    integer                   :: prev_line     ! Previous line with "base" obs. (on PIO)
    integer                   :: nbox          ! Total number of boxes (on PIO)
    integer                   :: nobs          ! Number of obs. in a box
    integer                   :: nl            ! Number of Lines in so-box
    integer                   :: nf            ! Number of FOVs in so_box

    ! Auxiliary information to identify obs. to be superobbed
    integer,      allocatable :: id_pio    (:) ! obs-ID of obs. to be superobbed (from all PEs)
    integer,      allocatable :: id        (:) ! obs-ID of obs. to be superobbed on this PE
    integer,      allocatable :: ind_pio   (:) ! Index of obs (in original t_radv on origin PE) from all PEs
    integer,      allocatable :: ind       (:) ! Index of obs in t_radv on this PE

    ! MPI communication: send (alltoall) obs required on other PEs
    integer,      allocatable :: n_snd_pio (:,:) ! Number of obs. to be sent between all PEs
    integer,      allocatable :: n_snd     (:)   ! Number of obs. to be send from this PE to other PEs
    integer,      allocatable :: n_rcv     (:)   ! Number of obs. to be received on this PE
    integer,      allocatable :: n_snd_    (:)   ! Number of entries in id_snd/pe_snd to be sent to PEs (from PIO)
    integer                   :: mx_snd          ! Max. number of obs. to be sent between PEs
    integer,      allocatable :: id_snd_pio(:,:) ! obs-IDs to be sent to other PEs (on PIO)
    integer,      allocatable :: id_snd    (:)   ! obs-IDs to send to other PEs (from this PE)
    integer,      allocatable :: pe_snd_pio(:,:) ! PEs to which the obs shall be sent (on PIO)
    integer,      allocatable :: pe_snd    (:)   ! PEs to which the obs shall be sent (from this PE)
    integer                   :: pe_rcv          ! Receiving PE


    ! Auxiliary array for MPI communication
    integer,      allocatable :: ibuf      (:) ! Array to be sent to other PEs
    integer,      allocatable :: displ     (:) ! displacements (for preparation of alltoall comm.)
    logical,      allocatable :: mask      (:) ! Mask of obs. to be kept/dismissed

    ! Sorting of obs (on PIO)
    integer,      allocatable :: ii        (:) ! Indices for sorted access to obs
    integer(i8),  allocatable :: ival      (:) ! Values to be sorted ascendingly

    integer                   :: iset, i, j, i_, j_ ! Loop indices
    integer                   :: n, stat

    ! Variables for additional check
    real(wp)                  :: distkm

    ! Shrink
    type(t_obs_block)         :: ob
    integer                   :: ib

    ! Access to satellite data is greatly simplified by transforming the data
    ! from the general dace-observation types to the "t_radv" type. We do the
    ! superobbing in the t_radv type and not in the general types. In this type
    ! each variable is stored in one array (for all obs). The last index in these
    ! arrays defines the FOV/spot. In the following, "index" usually means this
    ! index of a spot in the t_radv type.

    !   Grid of satellite data:
    !
    !   Scanline
    !      4  |   +(25) +(26) +(27) +(28) +(29) +(30) +(31) +(32)
    !         |
    !         |
    !      3  |   +(17) +(18) +(19) +(20) +(21) +(22) +(23) +(24)
    !         |
    !         |
    !      2  |   +(9)  +(10) +(11) +(12) +(13) +(14) +(15) +(16)
    !         |
    !         |
    !      1  |   +(1)  +(2)  +(3)  +(4)  +(5)  +(6)  +(7)  +(8)
    !         |
    !         !
    !         ------------------------------------------------------------ FOV
    !             1     2     3     4     5     6     7     8
    !
    ! We assume that we have sorted access to the data just as sketched by the
    ! (.)-bracketed numbers above. For the moment it is useful to assume that the
    ! bracketed numbers are the index of the obs.
    ! The boxes of obs to be superobbed ("so-boxes" in the following) are defined
    ! by an array "ibase", that contains the index of the "base" obs. of each box:
    !
    !   Scanline
    !      4  |   +(25) +(26) +(27) +(28) +(29) +(30) +(31) +(32)
    !         |    <17>  <17>  <19>  <19>  <21>  <21>  <23>  <23>
    !         |
    !      3  |   +(17) +(18) +(19) +(20) +(21) +(22) +(23) +(24)
    !         |    <17>  <17>  <19>  <19>  <21>  <21>  <23>  <23>
    !         |
    !      2  |   +(9)  +(10) +(11) +(12) +(13) +(14) +(15) +(16)
    !         |    <1>   <1>   <3>   <3>   <5>   <5>   <7>   <7>
    !         |
    !      1  |   +(1)  +(2)  +(3)  +(4)  +(5)  +(6)  +(7)  +(8)
    !         |    <1>   <1>   <3>   <3>   <5>   <5>   <7>   <7>
    !         !
    !         ------------------------------------------------------------ FOV
    !             1     2     3     4     5     6     7     8
    !
    ! where the <.>-bracketed numbers are the "ibase" array.
    ! The base obs (spots/FOVs) are kept in the superobbing, the other obs are dismissed.
    !
    ! Usually the obs are stored on different PEs. In order to chose the so-boxes
    ! we we have to know all obs. Thus the so-boxes are selected on dace%pio ("PIO"
    ! in the following). The information on the so-boxes has to be distributed to all
    ! the other PEs then. Let's introduce a boundary between PEs in the above example:
    !
    !   Scanline  <-----PE1-------> <-------------PE2------------>
    !      4  |   +(25) +(26) +(27)|+(28) +(29) +(30) +(31) +(32)
    !         |    <17>  <17>  <19>| <19>  <21>  <21>  <23>  <23>
    !         |                    -------------
    !      3  |   +(17) +(18) +(19) +(20) +(21)|+(22) +(23) +(24)
    !         |    <17>  <17>  <19>  <19>  <21>| <21>  <23>  <23>
    !         |                                |
    !      2  |   +(9)  +(10) +(11) +(12) +(13)|+(14) +(15) +(16)
    !         |    <1>   <1>   <3>   <3>   <5> | <5>   <7>   <7>
    !         |                          -------
    !      1  |   +(1)  +(2)  +(3)  +(4) |+(5)  +(6)  +(7)  +(8)
    !         |    <1>   <1>   <3>   <3> | <5>   <5>   <7>   <7>
    !         !
    !         ------------------------------------------------------------ FOV
    !             1     2     3     4     5     6     7     8
    !
    ! Let's assume, that the left area is PE1 and the right area is PE2.
    ! In this case for all so-boxes, that are crossed by the PE-boundary, those
    ! obs that are not on the PE of the "base" obs, must be sent to the PE
    ! of the "base" FOV, on which all the obs information will be merged.
    ! In the example above the following obs have to be sent between PE1 and PE2
    ! PE1 -> PE2: (13)  (5 is the "base" FOV of the so-box, that contains (13))
    ! PE2 -> PE1: (22), (29), (30) (all have the "base" FOV 21)
    !             (28)   ("base" FOV 19)
    !
    ! This is complicated by the fact, that the indices during the calculation of the
    ! so-boxes (on PIO, where all obs. from all PEs are contained in the t_radv arrays)
    ! are not the same as those that have to be used later, when the obs are merged on
    ! their respective PE (where only the obs known to the PE are contained in the arrays).
    ! During the calculation of the so-boxes (on PIO) the indices
    ! on the origin PE are stored in r% ind(:). These indices as well as "local" versions
    ! of the "ibase" array are distributed to the responsible PEs. For those obs, that
    ! have to be sent from their origin PE to another PE this index is not known during
    ! the calculation of the so-boxes. Thus, we also distribute the obs-IDs (t_radv%obsnum/
    ! t_spot%hd%id) to the responsible PEs. These are used to find the correct obs
    ! received from other PEs (and to do a safeguard consistency check of the indices).

    if (present(ldism)) ldism = .false.

    do iset = 1, n_set
      if ( rad_set(iset)%gopts%thin_superob_box_lines > 1 .or. &
           rad_set(iset)%gopts%thin_superob_box_fovs  > 1) then

        ! 1. Calculate superobbing boxes
        ! Get meta information of ALL obs in t_radv on PIO:
        call fill_rad(rad(:), rad_set(iset:iset), obs, y=y, pe_gather=dace%pio, &
                      lmon=.true., lpartial=.true.)
        r => rad(1)
        n = r%n_rec
        call p_bcast(n, dace%pio)
        if (n <= 0) cycle
        if (associated(r%pe)) r%pe = r%pe + 1 ! pe from 0..npe-1, but here we count 1..npe

        if (dace% lpio) then
          ! Calc. boxes on pio, since meta data of all obs. must be known to calc. boxes
          nrec = r%n_rec
          nl = r%i%gopts%thin_superob_box_lines
          nf = r%i%gopts%thin_superob_box_fovs

          write(*,*)
          write(*,'("Superobbing/Thinning of instrument ",I3.3," on sat. ",I3.3)') &
               rad_set(iset)%grid,rad_set(iset)%satid
          write(*,'(2x,"Boxes: ",I2," lines, ",I2," FOVs")') nl,nf
          select case(r%i%gopts%thin_superob_mode)
          case(1)
            write(*,'(2x,"Mode : 1 (Thinning)")')
          case(2)
            write(*,'(2x,"Mode : 2 (Superobbing)")')
          case default
            call finish('thin_superob_tovs','Unknown mode!')
          end select

          ! Set up index array for sorted access to r. Sort by (1.) scanline, (2.) FOV
          allocate(ii(nrec),ival(nrec))
          j = maxval(r%fov(1:nrec), mask=(r%fov(1:nrec)/=rinvalid))
          if (j > 0) then
            j = ceiling(log10(j*1._wp))
          else
            j = 0
          endif
          j = 10**j
          ival(1:nrec) = r%scanl(1:nrec) * j + r%fov(1:nrec)
          ii(1:nrec) = (/ (i,i=1,nrec) /)
          call sortrx(nrec, ival, ii)  !now ival(1) is smallest number in ival, ival(nrec) largest
          deallocate(ival)

          ! Determine superobbing boxes
          allocate(ibase(nrec))
          ibase(:) = 0
          nbox = 0
          prev_line = -1
          do i_ = 1, nrec
            i = ii(i_)
            if (ibase(i) > 0) cycle    ! this pixel has already a base pixel
            ! Constraint on difference between scanline numbers of "base" FOVs of so-boxes.
            ! This is not necessary, but gives a more homogeneous horizontal coverage:
            if (prev_line > 0) then
              if (mod(r%scanl(i)-prev_line, nl) /= 0) cycle
            end if

            nbox = nbox + 1
            ibase(i) = i
            nobs = 1
            ! Search forward for FOVs, that are within box
            ! obs i is the "lower, right" edge of box (scanl is y-axis, fov is x_axis)
            do j_ = i_, nrec
              j = ii(j_)
              if (r%scanl(j) - r%scanl(i) > nl-1) exit   !outside superob radius for lines
              if (r%fov(j) >= r%fov(i) .and. r%fov(j) - r%fov(i) <= nf-1 .and. &
                   ibase(j) <= 0) then
                 ! Check distance in time and space between pixels (currently set to 50km and 15 minutes)
                 distkm = distance(r%dlat(j),r%dlon(j),r%dlat(i),r%dlon(i))
                 if (distkm > r%i%gopts%max_size_sobox ) then
                    write(*, '(2x, "Attention, distance between 2 superobbing pixels too large: ", F6.2)') &
                                       distkm
                 else
                    if (abs(r%time(j) -r%time(i)) > r%i%gopts%max_timediff_sobox) then
                       write(*, '(2x, "Attention, time of superobbing pixels differs too much: ", F6.2)') &
                                       abs(r%time(j) -r%time(i))
                    else
                      nobs = nobs + 1
                      ibase(j) = i
                    end if
                 end if
              end if
            end do  !search loop for boxes in so-box
            if (nobs /= nl * nf) then
              ! Delete boxes, that do not have the desired extent?
              where(ibase(ii(i_:)) == i) ibase(ii(i_:)) = 0   ! reset ibase to zero
              nbox = nbox - 1   ! if not desired extent don't count this box
              if (nobs > nl*nf) write(0,*) 'thin_superob_tovs: too large box:',nobs,'>',nl*nf,'i=',i
            else
              prev_line = r%scanl(i)
            end if
          end do
          deallocate(ii)

          write(*,'(2x,"Number of so-boxes: ",I9)') nbox
          write(*,'(2x,"Number of obs in so-boxes: ",I10,"  (",F6.2,"%)")') &
               nbox * nl * nf, (nbox * nl * nf * 100._wp) / real(r%n_rec)

        end if

        ! 2. Distribute info on so-boxes to the responsible PEs
        allocate(ns_pio(dace%npe))
        ns_pio = 0
        if (dace%lpio) then
          do i = 1, nrec
            if (ibase(i) > 0) then
              ns_pio(r%pe(ibase(i))) = ns_pio(r%pe(ibase(i))) + 1
            end if
          end do
          n = sum(ns_pio)
          allocate(displ(dace%npe), id_pio(n), ibase_pio(n), ind_pio(n))
          displ(1) = 0
          do i = 2, dace%npe
            displ(i) = displ(i-1) + ns_pio(i-1)
          end do
          ns_pio = 0
          do i = 1, nrec
            if (ibase(i) > 0) then
              pe_rcv = r%pe(ibase(i))
              ns_pio(pe_rcv) = ns_pio(pe_rcv) + 1
              id_pio   (displ(pe_rcv)+ns_pio(pe_rcv)) = r% obsnum(i)
              ind_pio  (displ(pe_rcv)+ns_pio(pe_rcv)) = r% ind(i)        ! index on origin PE
              ibase_pio(displ(pe_rcv)+ns_pio(pe_rcv)) = r% ind(ibase(i)) ! index of "base" FOV on origin PE
              if (r%pe(i) /= r%pe(ibase(i))) then
                ! Spot and "base" spot are NOT on the same PE. The spot will be appended to
                ! t_radv array. Thus, the index (on the receiving PE!) is not known at the
                ! moment.
                ind_pio(displ(pe_rcv)+ns_pio(pe_rcv)) = -1
              end if
            end if
          end do
          deallocate(displ)
        else
          allocate(id_pio(0), ibase_pio(0), ind_pio(0))
        end if
        call p_scatterv(ns_pio, (/(1,i=1,dace%npe)/), ns_(1:1), dace%pio)
        ns = ns_(1)
        allocate(ibase_(ns), id(ns), ind(ns))
        call p_scatterv(ind_pio  (:), ns_pio, ind   (:), dace%pio)
        call p_scatterv(id_pio   (:), ns_pio, id    (:), dace%pio)
        call p_scatterv(ibase_pio(:), ns_pio, ibase_(:), dace%pio)
        deallocate(ibase_pio, id_pio, ind_pio, ns_pio)

        ! 3. Prepare alltoall communication of obs, that are required from other PEs
        ! determine and distribute numbers of spots to be sent/received
        if (dace%lpio) then
          allocate(n_snd_pio(dace%npe,dace%npe))
          n_snd_pio = 0

          do i = 1, nrec
            if (ibase(i) > 0) then
              if (r%pe(i) /= r%pe(ibase(i))) then
                ! Base spot is on other PE than current spot
                n_snd_pio(r%pe(i), r%pe(ibase(i))) = n_snd_pio(r%pe(i), r%pe(ibase(i))) + 1
              end if
            end if
          end do
          write(*,'(2x,"Number of FOVs to be distributed: ",I9)') sum(n_snd_pio)
        else
          allocate(n_snd_pio(0,0))
        end if
        allocate(n_rcv(dace%npe), n_snd(dace%npe))
        call p_scatterv(reshape(n_snd_pio,(/size(n_snd_pio)/)), &
             (/(dace%npe,i=1,dace%npe)/), n_rcv(:), dace%pio)
        call p_scatterv(reshape(transpose(n_snd_pio),(/size(n_snd_pio)/)), &
             (/(dace%npe,i=1,dace%npe)/), n_snd(:), dace%pio)
        ! Set up spot lists
        if (dace%lpio) then
          ! max. number of spots to be sent by any PE
          mx_snd = 0
          do i = 1, dace%npe
            mx_snd = max(mx_snd, sum(n_snd_pio(i,:)))
          end do
          allocate(n_snd_(dace%npe), id_snd_pio(dace%npe,mx_snd), pe_snd_pio(dace%npe,mx_snd))
          ! for each PE: set up list of spots to be sent and list of target PEs
          n_snd_ = 0
          do i = 1, nrec
            if (ibase(i) > 0) then
              if (r%pe(i) /= r%pe(ibase(i))) then
                n_snd_(r%pe(i)) = n_snd_(r%pe(i)) + 1
                id_snd_pio(r%pe(i), n_snd_(r%pe(i))) = r% obsnum(i)
                pe_snd_pio(r%pe(i), n_snd_(r%pe(i))) = r% pe(ibase(i))
              end if
            end if
          end do
          deallocate(ibase)
        else
          allocate(n_snd_(dace%npe))
          n_snd_ = 0
        end if
        deallocate(n_snd_pio)

        ! distribute list of spots to be sent and list of target PEs
        allocate(ibuf(sum(n_snd_)), id_snd(sum(n_snd)), pe_snd(sum(n_snd)))
        if (dace%lpio) then
          j = 0
          do i = 1, dace%npe
            ibuf(j+1:j+n_snd_(i)) = id_snd_pio(i,1:n_snd_(i))
            j = j + n_snd_(i)
          end do
        end if
        call p_scatterv(ibuf, n_snd_, id_snd, dace%pio)
        if (dace%lpio) then
          j = 0
          do i = 1, dace%npe
            ibuf(j+1:j+n_snd_(i)) = pe_snd_pio(i,1:n_snd_(i))
            j = j + n_snd_(i)
          end do
        end if
        call p_scatterv(ibuf, n_snd_, pe_snd, dace%pio)
        deallocate(ibuf,n_snd_)
        if (dace%lpio) deallocate(pe_snd_pio, id_snd_pio)

        call destruct(r)

        ! 4. Fill t_radv and distribute content of t_radv as required for superobbing
        ! Fill t_radv with additional space for spots to be received
        call fill_rad(rad(:), rad_set(iset:iset), obs, y=y, lrttov=.true., nadd=(/sum(n_rcv)/), &
                      lpartial=.true.)
        r => rad(1)
        ! Distribute spots required on other PEs
        call alltoall_radv(r, id_snd, pe_snd, n_snd, n_rcv)
        if (associated(r%pe) .and. r%n_rec > 0) r%pe = r%pe + 1 ! pe from 0..npe-1, but here we count 1..npe
        ! 5. Superobbing, i.e. merge obs, in each so-box.
        allocate(mask(r%n_rec))
        call thin_superob_radv(r, ind, id, ibase_, mask, status=stat)
        if (stat /= 0) call finish('thin_superob_tovs', 'thin_superob_radv failed')

        ! 6. Store results (in t_radv) into general dace types, dismiss obs that are not a "base" obs.
        call radv_obs2(r, obs, y, mask)

        if (present(ldism)) ldism = ldism .or. .not.all(mask)

        if (dace%lpio) write(*,*)
      end if
    end do

  contains

    ! function which returns distance between two coordinates in km
    function distance(lat1,lon1,lat2,lon2) result(d)
      real(kind=wp), intent(in) :: lat1,lon1,lat2,lon2
      real(kind=wp) :: d
      real(kind=wp) :: x1(3), x2(3)
      x1 = cartesian(lat1,lon1)
      x2 = cartesian(lat2,lon2)
      d = acos(dot_product(x1,x2)) * 6371._wp
    end function distance

    !function which transfers the coordinates into a cartesian coordinate system
    function cartesian(lat,lon) result(x)
      real(kind=wp), intent(in) :: lat,lon
      real(kind=wp) :: x(3)
      x(1) = dcos(lat*d2r) * dsin(lon*d2r)
      x(2) = dcos(lat*d2r) * dcos(lon*d2r)
      x(3) = dsin(lat*d2r)
    end function cartesian

  end subroutine thin_superob_tovs


  !> Store t_radv into general dace obs type
  !! Possibly we might merge this with radv_obs in the future??
  subroutine radv_obs2(r, obs, y, mask)
    type(t_radv),    intent(in)    :: r
    type(t_obs_set), intent(inout) :: obs   ! observation data (3dvar type)
    type(t_vector),  intent(inout) :: y     ! background observation space
    logical,         intent(in)    :: mask(:)

    type(t_spot), pointer :: s => null()
    type(t_tovs) :: tovs
    integer :: is, ib, i, j, k, l

    do i = 1, r%n_rec
      if (r%pe(i)-1 /= dace%pe) cycle
      ib = r% i_box(i)
      is = r% i_reprt(i)
      s => obs% o(ib)% spot(is)

      if (s%hd%obstype /= OT_RAD) call finish('radv_obs2', 'bad obstype')
      if (r%obsnum(i)  /= s% hd% id) call finish('radv_obs2', 'mismatch of obsnum')
      if (mask(i)) then
        ! Keep observation, store stuff from r into general obs/spot types
        call load(obs% o(ib), s, tovs, tovs_io=TTOVS_BASE+TTOVS_EMIS+TTOVS_BTCS)
        if (associated(r%dlon   )) s% col% c% dlon = r% dlon   (i)
        if (associated(r%dlat   )) s% col% c% dlat = r% dlat   (i)
        ! s% hd% time = .. !TODO inverse of r% time(i) = imi(s% hd% time - ana_time)
        if (associated(r%nwc_flg)) tovs% nwc_flg   = r% nwc_flg(i)
        if (associated(r%fov    )) s% phase        = r% fov    (i)
        if (associated(r%scanl  )) tovs% scanl     = r% scanl  (i)
        if (associated(r%mdlsfc )) s% mdlsfc       = r% mdlsfc (i)
        if (associated(r%stzen  )) s% stzen        = r% stzen  (i)
        if (associated(r%stazi  )) tovs% boa       = r% stazi  (i)
        if (associated(r%sunzen )) tovs% soza      = r% sunzen (i)
        if (associated(r%sunazi )) tovs% soa       = r% sunazi (i)
        if (associated(r%shgt   )) s% gp_bg        = r% shgt (1,i) * (gacc*1.e3_wp)
        if (associated(r%stype  )) tovs%rt_stype(1)= r% stype(1,i)
        if (associated(r%r_state)) s%use%state     = r% r_state(i)

        k = 0
        do j = 1, r%i%n_chan
          if (r%valid(j,i)) then
            k = k + 1
            l = k + s%o%i
            if (k > s%o%i + s%o%n) call finish('rav_obs2','invalid indices')
            if (ldeb(s)) write(usd,*) dpref,'radv_obs2 chan',j,k
            if (associated(r%bcor_   )) obs% o(ib)% body(l)% bc         = -r% bcor_  (j,i)
            if (associated(r%bt_bcor )) obs% o(ib)% body(l)% o          =  r% bt_bcor(j,i)
            ! r%bt_obs is just r%bt_bcor + r%bcor_
            if (associated(r%state   )) obs% o(ib)% body(l)% use% state = r% state   (j,i)
            if (associated(r%flags   )) obs% o(ib)% body(l)% use% flags = r% flags   (j,i)
            if (associated(r%bt_fg   )) y% s(ib)% x(l)                  = r% bt_fg   (j,i)
            if (associated(r%emiss   )) tovs%emis(k)                    = r% emiss   (j,i)
            if (associated(r%bt_fg_cs)) tovs%bt_cs(k)                   = r% bt_fg_cs(j,i)
          end if
        end do

        ! Please add further variables if required ...

        call store(obs% o(ib), s, tovs, tovs_io=TTOVS_BASE+TTOVS_EMIS+TTOVS_BTCS)
        call destruct(tovs)
      else
        ! Do not keep observation
        call decr_use(s%use, STAT_DISMISS, CHK_MERGE)
        do k = s%o%i+1,s%o%i+s%o%n
          call decr_use(obs% o(ib)% body(k)% use, STAT_DISMISS, CHK_MERGE)
        end do
      end if

    end do

  end subroutine radv_obs2


  !> Routine to perform alltoall communication for t_radv content.
  !! The received spots are added at the end of the arrays in t_radv
  !! Currently not all possible content of t_radv is broadcasted. Please
  !! add further  variables if required.
  subroutine alltoall_radv(r, id_snd, pe_snd, n_snd, n_rcv)
    type(t_radv), intent(inout) :: r
    integer,      intent(in)    :: id_snd(:)
    integer,      intent(in)    :: pe_snd(:)
    integer,      intent(in)    :: n_snd(:)
    integer,      intent(in)    :: n_rcv(:)

    integer :: i_s(sum(n_snd)) ! indices of spots to be sent
    integer :: n(dace%npe)     ! number of spots to be sent to each PE
    integer :: displ(dace%npe) ! "displacement" of each PE in arrays to be sent
    integer :: pe, nr, ns, i, j, k, nch
    logical :: lfound(size(id_snd))
    real(wp), allocatable :: r_rcv(:)
    integer,  allocatable :: i_rcv(:)

    displ(1) = 0
    do i = 2, dace%npe
      displ(i) = displ(i-1) + n_snd(i-1)
    end do

    lfound = .false.
    n = 0
    do i = 1, r%n_rec
      k = 0
      do j = 1, size(id_snd)
        if (r%obsnum(i) == id_snd(j)) then
          k = j
          lfound(j) = .true.
          exit
        end if
      end do
      if (k > 0) then
        pe = pe_snd(k)
        n(pe) = n(pe) + 1
        i_s(displ(pe) + n(pe)) = i
      end if
    end do
    ! Consistency check
    do j = 1, size(id_snd)
      if (.not.lfound(j)) then
        write(0,*) 'j=',j,' id_snd(j)=',id_snd(j)
        call finish('alltoall_radv', 'Did not find spot to be distributed')
      end if
    end do
    nr = sum(n_rcv)
    ns = sum(n_snd)

    ! Make sure, that info on PE origin is set
    if (.not.associated(r%pe)) then
      allocate(r%pe(r%n_rec+nr))
      r%pe(1:r%n_rec) = dace%pe
    end if
    if (.not.associated(r%ind)) then
      allocate(r%ind(r%n_rec+nr))
      r%ind(1:r%n_rec) = (/ (j, j=1,r%n_rec) /)
    end if

    ! Scalar quantities
    if (associated(r%pe     )) call p_alltoall(r%pe     (i_s(:)), r%pe     (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%ind    )) call p_alltoall(r%ind    (i_s(:)), r%ind    (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%obsnum )) call p_alltoall(r%obsnum (i_s(:)), r%obsnum (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%i_box  )) call p_alltoall(r%i_box  (i_s(:)), r%i_box  (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%i_reprt)) call p_alltoall(r%i_reprt(i_s(:)), r%i_reprt(r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%dlat   )) call p_alltoall(r%dlat   (i_s(:)), r%dlat   (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%dlon   )) call p_alltoall(r%dlon   (i_s(:)), r%dlon   (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%time   )) call p_alltoall(r%time   (i_s(:)), r%time   (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%nwc_flg)) call p_alltoall(r%nwc_flg(i_s(:)), r%nwc_flg(r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%r_state)) call p_alltoall(r%r_state(i_s(:)), r%r_state(r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%scanl  )) call p_alltoall(r%scanl  (i_s(:)), r%scanl  (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%fov    )) call p_alltoall(r%fov    (i_s(:)), r%fov    (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%mdlsfc )) call p_alltoall(r%mdlsfc (i_s(:)), r%mdlsfc (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%op_na  )) call p_alltoall(r%op_na  (i_s(:)), r%op_na  (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%stzen  )) call p_alltoall(r%stzen  (i_s(:)), r%stzen  (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%stazi  )) call p_alltoall(r%stazi  (i_s(:)), r%stazi  (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%sunzen )) call p_alltoall(r%sunzen (i_s(:)), r%sunzen (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)
    if (associated(r%sunazi )) call p_alltoall(r%sunazi (i_s(:)), r%sunazi (r%n_rec+1:r%n_rec+nr), dace%comm, n_snd, n_rcv)

    ! Dimension n_chan
    nch = r%i%n_chan
    allocate(r_rcv(nr * nch), i_rcv(nr * nch))
    if (associated(r%bt_obs)) then
      call p_alltoall(reshape(r%bt_obs(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%bt_obs(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    if (associated(r%bt_bcor)) then
      call p_alltoall(reshape(r%bt_bcor(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%bt_bcor(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    if (associated(r%bcor_)) then
      call p_alltoall(reshape(r%bcor_(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%bcor_(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    if (associated(r%bt_fg)) then
      call p_alltoall(reshape(r%bt_fg(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%bt_fg(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    if (associated(r%bt_fg_cs)) then
      call p_alltoall(reshape(r%bt_fg_cs(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%bt_fg_cs(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    ! if (associated(r%not_rej)) then
    !   call p_alltoall(reshape(r%not_rej(:,i_s(:)), (/ns*nch/)), i_rcv, dace%comm, n_snd*nch, n_rcv*nch)
    !   r%not_rej(:,r%n_rec+1:r%n_rec+nr) = reshape(i_rcv, (/nch,nr/))
    ! end if
    if (associated(r%cloudy)) then
      call p_alltoall(reshape(r%cloudy(:,i_s(:)), (/ns*nch/)), i_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%cloudy(:,r%n_rec+1:r%n_rec+nr) = reshape(i_rcv, (/nch,nr/))
    end if
    if (associated(r%state)) then
      call p_alltoall(reshape(r%state(:,i_s(:)), (/ns*nch/)), i_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%state(:,r%n_rec+1:r%n_rec+nr) = reshape(i_rcv, (/nch,nr/))
    end if
    if (associated(r%flags)) then
      call p_alltoall(reshape(r%flags(:,i_s(:)), (/ns*nch/)), i_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%flags(:,r%n_rec+1:r%n_rec+nr) = reshape(i_rcv, (/nch,nr/))
    end if
    if (associated(r%emiss)) then
      call p_alltoall(reshape(r%emiss(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%emiss(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    ! if (associated(r%valid)) then
    !   call p_alltoall(reshape(r%valid(:,i_s(:)), (/ns*nch/)), i_rcv, dace%comm, n_snd*nch, n_rcv*nch)
    !   r%valid(:,r%n_rec+1:r%n_rec+nr) = reshape(i_rcv, (/nch,nr/))
    ! end if
    deallocate(r_rcv, i_rcv)


    ! Dimension n_sens
    nch = r%i%n_sens
    allocate(r_rcv(nr * nch), i_rcv(nr * nch))
    if (associated(r%shgt)) then
      call p_alltoall(reshape(r%shgt(:,i_s(:)), (/ns*nch/)), r_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%shgt(:,r%n_rec+1:r%n_rec+nr) = reshape(r_rcv, (/nch,nr/))
    end if
    if (associated(r%stype)) then
      call p_alltoall(reshape(r%stype(:,i_s(:)), (/ns*nch/)), i_rcv, dace%comm, n_snd*nch, n_rcv*nch)
      r%stype(:,r%n_rec+1:r%n_rec+nr) = reshape(i_rcv, (/nch,nr/))
    end if
    deallocate(r_rcv, i_rcv)

    r%n_rec = r%n_rec + nr

  end subroutine alltoall_radv

  elemental subroutine destruct_rtstat_spt(s)
    type(t_keep_rtstat_spt), intent(inout) :: s
    if (associated(s%i)) deallocate(s%i)
    s%id    = -1
    s%n     = 0
    s%nbits = 0
  end subroutine destruct_rtstat_spt

  elemental subroutine destruct_rtstat(r)
    type(t_keep_rtstat), intent(inout) :: r
    if (associated(r%spt)) then
      call destruct(r%spt)
      deallocate(r%spt)
    end if
    r%nspot  = 0
  end subroutine destruct_rtstat

  subroutine store_rtstat(s,larr1, larr2, larr3, larr4)
    type(t_keep_rtstat_spt), intent(inout)                 :: s
    logical,                 intent(in),   target, optional :: larr1(:,:)
    logical,                 intent(in),   target, optional :: larr2(:,:)
    logical,                 intent(in),   target, optional :: larr3(:,:)
    logical,                 intent(in),   target, optional :: larr4(:,:)

    logical, pointer     :: larr1p(:,:),larr2p(:,:),larr3p(:,:),larr4p(:,:)
    logical, allocatable :: larr_(:)
    logical              :: l_get
    integer :: nbits(0:4)
    integer :: n, i, j, k

    l_get = .false.
    nbits = s%nbits
    if (present(larr1)) then
      larr1p => larr1
      nbits(1) = size(larr1p)
    elseif (nbits(1) > 0) then
      allocate(larr1p(nbits(1)/2,2))
      l_get = .true.
    end if
    if (present(larr2)) then
      larr2p => larr2
      nbits(2) = size(larr2p)
    elseif (nbits(2) > 0) then
      allocate(larr2p(nbits(2),1))
      l_get = .true.
    end if
    if (present(larr3)) then
      larr3p => larr3
      nbits(3) = size(larr3p)
    elseif (nbits(3) > 0) then
      allocate(larr3p(nbits(3),1))
      l_get = .true.
    end if
    if (present(larr4)) then
      larr4p => larr4
      nbits(4) = size(larr4p)
    elseif (nbits(4) > 0) then
      allocate(larr4p(nbits(4),1))
      l_get = .true.
    end if
    nbits(0) = sum(nbits(1:))
    if (l_get) call load(s,larr1p,larr2p,larr3p,larr4p,&
         lget1=.not.present(larr1),lget2=.not.present(larr2),&
         lget3=.not.present(larr3),lget4=.not.present(larr4))

    n = ceiling(nbits(0)/(nbits_int*1.))
    if (n /= s%n) then
      if (associated(s%i)) deallocate(s%i)
      allocate(s%i(n))
    end if
    s%n = n
    s%nbits = nbits
    allocate(larr_(nbits(0)))
    if (nbits(1) > 0) larr_(1                :nbits(1)       ) = reshape(larr1p, (/nbits(1)/))
    if (nbits(2) > 0) larr_(nbits(1)+1       :sum(nbits(1:2))) = reshape(larr2p, (/nbits(2)/))
    if (nbits(3) > 0) larr_(sum(nbits(1:2))+1:sum(nbits(1:3))) = reshape(larr3p, (/nbits(3)/))
    if (nbits(4) > 0) larr_(sum(nbits(1:3))+1:sum(nbits(1:4))) = reshape(larr4p, (/nbits(4)/))
    j = 0
    i_loop: do i = 1, n
      s%i(i) = 0
      do k = 0, nbits_int-1
        j = j + 1
        if (j > nbits(0)) exit i_loop
        if (larr_(j)) s%i(i) = ibset(s%i(i), k)
      end do
    end do i_loop

    if (.not.present(larr1) .and. nbits(1) > 0) deallocate(larr1p)
    if (.not.present(larr2) .and. nbits(2) > 0) deallocate(larr2p)
    if (.not.present(larr3) .and. nbits(3) > 0) deallocate(larr3p)
    if (.not.present(larr4) .and. nbits(4) > 0) deallocate(larr4p)

  end subroutine store_rtstat

  subroutine load_rtstat(s,larr1,larr2,larr3,larr4,lget1,lget2,lget3,lget4)
    type(t_keep_rtstat_spt), intent(inout)          :: s
    logical,                 intent(out),  optional :: larr1(:,:)
    logical,                 intent(out),  optional :: larr2(:,:)
    logical,                 intent(out),  optional :: larr3(:,:)
    logical,                 intent(out),  optional :: larr4(:,:)
    logical,                 intent(in),   optional :: lget1
    logical,                 intent(in),   optional :: lget2
    logical,                 intent(in),   optional :: lget3
    logical,                 intent(in),   optional :: lget4

    logical, allocatable :: larr_(:)
    logical :: lget1_,lget2_,lget3_,lget4_
    integer :: nbits_out(0:4), is, ie
    integer :: i, j, k

    lget1_ = present(larr1)
    if (lget1_) then
      if (present(lget1)) lget1_ = lget1
    end if
    lget2_ = present(larr2)
    if (lget2_) then
      if (present(lget2)) lget2_ = lget2
    end if
    lget3_ = present(larr3)
    if (lget3_) then
      if (present(lget3)) lget3_ = lget3
    end if
    lget4_ = present(larr4)
    if (lget4_) then
      if (present(lget4)) lget4_ = lget4
    end if

    nbits_out = 0
    if (lget1_) nbits_out(1) = size(larr1)
    if (lget2_) nbits_out(2) = size(larr2)
    if (lget3_) nbits_out(3) = size(larr3)
    if (lget4_) nbits_out(4) = size(larr4)
    nbits_out(0) = sum(nbits_out(1:))

    if (nbits_out(0) == 0) return
    if (any(nbits_out(1:) > 0 .and. s%nbits(1:) > 0 .and. nbits_out(1:) /= s%nbits(1:))) &
         call finish('load_rtstat', 'nbits /= s%nbits')

    if (any(s%nbits > 0)) then
      allocate(larr_(s%nbits(0)))
      j = 0
      i_loop: do i = 1, s%n
        do k = 0, nbits_int-1
          j = j + 1
          if (j > s%nbits(0)) exit i_loop
          larr_(j) = btest(s%i(i), k)
        end do
      end do i_loop
    end if
    if (lget1_) then
      if (s%nbits(1) > 0) then
        is = 1
        ie = s%nbits(1)
        larr1 = reshape(larr_(is:ie), shape(larr1))
      else
        larr1 = .false.
      end if
    end if
    if (lget2_) then
      if (s%nbits(2) > 0) then
        is = s%nbits(1) + 1
        ie = sum(s%nbits(1:2))
        larr2 = reshape(larr_(is:ie), shape(larr2))
      else
        larr2 = .false.
      end if
    end if
    if (lget3_) then
      if (s%nbits(3) > 0) then
        is = sum(s%nbits(1:2)) + 1
        ie = sum(s%nbits(1:3))
        larr3 = reshape(larr_(is:ie), shape(larr3))
      else
        larr3 = .false.
      end if
    end if
    if (lget4_) then
      if (s%nbits(4) > 0) then
        is = sum(s%nbits(1:3)) + 1
        ie = sum(s%nbits(1:4))
        larr4 = reshape(larr_(is:ie), shape(larr4))
      else
        larr4 = .false.
      end if
    end if

  end subroutine load_rtstat

  subroutine print_t_jac(unit, tj)
     integer,     intent(in) :: unit
     type(t_jac), intent(in) :: tj

     write(unit,*) 'ps_llev        ',tj%ps_llev
     write(unit,*) 'il             ',tj%il
     write(unit,*) 'iu             ',tj%iu
     write(unit,*) 'wl             ',tj%wl
     write(unit,*) 'wu             ',tj%wu
     write(unit,*) 'd_lev          ',tj%d_lev
     write(unit,*) 'd_tv           ',tj%d_tv
     write(unit,*) 'd_rh           ',tj%d_rh
     if (associated(tj%dt_tv)) write(unit,*) 'dt_tv (:)      ',tj%dt_tv
     if (associated(tj%dt_rh)) write(unit,*) 'dt_rh (:)      ',tj%dt_rh
     if (associated(tj%dq_tv)) write(unit,*) 'dq_tv (:)      ',tj%dq_tv
     if (associated(tj%dq_rh)) write(unit,*) 'dq_rh (:)      ',tj%dq_rh
     if (associated(tj%de_dpc)) write(unit,*) 'de_dpc(:,:)    ',tj%de_dpc
     write(unit,*) 'dt_tv_2m       ',tj%dt_tv_2m
     write(unit,*) 'dt_rh_2m       ',tj%dt_rh_2m
     write(unit,*) 'dq_tv_2m       ',tj%dq_tv_2m
     write(unit,*) 'dq_rh_2m       ',tj%dq_rh_2m
     write(unit,*) 'sigma_var_tskin',tj%sigma_var_tskin
     write(unit,*) 'dclt_dum       ',tj%dclt_dum
     write(unit,*) 'dclf_dum       ',tj%dclf_dum
     write(unit,*) 'dsnf_dum       ',tj%dsnf_dum
  end subroutine print_t_jac

  !---------------------------------------------------------------
  subroutine prep_H_btcs (H_btcs, obs, op_flag)
  !---------------------------------------------------------
  ! gather model-equivalent clear-sky brightness temperatures in a t_vector
  !---------------------------------------------------------
  type (t_vector)   ,intent(inout) :: H_btcs
  type(t_obs_set)   ,intent(in)    :: obs      ! observation data
  integer           ,intent(out)   :: op_flag(:)
  !-------------------------------------------------------
  ! add verification data to feedback file
  ! (generic routine for deterministic and ensemble state)
  !-------------------------------------------------------
    integer                   :: i0, i1
    integer                   :: ib, is, isrc
    type(t_tovs)              :: ttovs
    type (t_rad_set), pointer :: rs => null()

    op_flag(:) = -1

    !----------------------------------------------------------------
    ! extract clear-sky brightness temperatures from t_tovs structure
    ! and save to t_vector (depending on namelist switch rad_out)
    !----------------------------------------------------------------
     do ib = 1, H_btcs% n_s
      if (obs% o(ib)% pe == dace% pe) then
        if (obs% o(ib)% n_obs > 0) then

          !write also clear-sky TBs in fdbk-file?
          do is = 1, obs% o(ib)% n_spot
            if (obs% o(ib)% spot(is)% hd% obstype == OT_RAD .and. &
                  obs% o(ib)% spot(is)% p% n > 0) then
              call load(obs% o(ib), obs% o(ib)% spot(is), ttovs, rs=rs, &
                           tovs_io=TTOVS_BTCS)
              if (any(btest(rs%iopts(1:rs%n_instr)%rad_out, OUT_CSB))) then
                i0=obs% o(ib)% spot(is)% o% i+1
                i1=i0 + obs% o(ib)% spot(is)% o %n -1
                H_btcs% s(ib)%x(i0:i1) = ttovs%bt_cs
                isrc = obs% o(ib)% spot(is)% hd% mon_file
                if (isrc > 0 .and. isrc <= n_source) op_flag(isrc) = OF_BT_CLEAR_SKY
              end if
            end if
          end do
        endif
      endif
    end do

    !----------------------
    ! gather data on I/O PE
    !----------------------
    op_flag = p_max(op_flag)
    if (any(op_flag == OF_BT_CLEAR_SKY)) call gather (H_btcs, dest=dace%pio)

  end subroutine prep_H_btcs
!-----------------------------------------------------------------------------------------


! specific MPI-bcast for derived type t_rtstat_opts
#undef  VECTOR
#undef  DERIVED
#define DERIVED type(t_rtstat_opts)
#define p_bcast_DERIVED bcast_rtstat_opts
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_bcast_DERIVED

! specific MPI-bcast for derived type t_trg_clim_trend
#undef  VECTOR
#undef  DERIVED
#define DERIVED type(t_trg_clim_trend)
#define p_bcast_DERIVED bcast_trg_clim_trend
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE

#if (_RTTOV_VERSION >= 12)
! specific p_alltoall routines for derived type t_wr_opd
#undef  VECTOR
#define DERIVED type(t_wr_opd)
#define p_alltoall_DERIVED p_alltoall_t_wr_opd
#include "p_alltoall_derived.incf"
#undef  DERIVED
#undef  p_alltoall_DERIVED
#undef  MPI_TYPE
#endif
end module mo_tovs
