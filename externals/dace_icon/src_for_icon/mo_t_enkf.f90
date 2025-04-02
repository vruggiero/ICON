!
!+ Ensemble Kalman Filter (EnKF) and Localised EnKF (LETKF) utility routines
!
MODULE mo_t_enkf
!
! Description:
!   Ensemble Kalman Filter (EnKF) and Localised EnKF (LETKF) utility routines.
!
! Current Code Owner: DWD, Hendrik Reich
!    phone: +49 69 8062 4943
!    fax:   +49 69 8062 3721
!    email: hendrik.reich@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Hendrik Reich
!  Changes for COSMO LETKF
! V1_5         2009/05/25 Hendrik Reich
!  Optimisations for SX-9; read observations from feedback file
! V1_6         2009/06/10 Andreas Rhodin
!  for GME perform LETKF calculations on common vertical grid
! V1_7         2009/08/24 Andreas Rhodin
!  Changes for COSMO; adapt subroutine nominal_height to AMSUB
! V1_8         2009/12/09 Hendrik Reich
!  Changes for COSMO
! V1_9         2010/04/20 Hendrik Reich
!  comp. ana. on C-grid;apply hum bounds to qr,qs,qg.
!  obs. diag. implemented (COSMO only).
!  remove factor 0.25 in random ensemble generation.
!  fix default values of ane_fname,fce_fname for GME.
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Hendrik Reich
!  grid point diagnostics implemented
! V1_13        2011/11/01 Hendrik Reich
!  split modules mo_t_enkf, mo_letkf; various changes
! V1_15        2011/12/06 Andreas Rhodin
!  additional diagnostic: print inventory of anaysis input file
! V1_19        2012-04-16 Hendrik Reich
!  revise vertical localisation
!  option to read multiple fof-files / member
!  vectorise ensemble mean and spread calculation for LETKF
! V1_20        2012-06-18 Andreas Rhodin
!  changes for forecast sensitivities
! V1_22        2013-02-13 Andreas Rhodin
!  changes for RADAR operator
! V1_23        2013-03-26 Andreas Rhodin
!  subroutine apply_H: changed to be called from verification program
! V1_26        2013/06/27 Andreas Rhodin
!  reorganise LETKF for multistep analysis: observation space analysis
! V1_27        2013-11-08 Andreas Rhodin
!  adapt LETKF to ICON
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF; changes for SEVIRI in LETKF
! V1_29        2014/04/02 Hendrik Reich
!  relaxation subroutines (pert./spread) added
!  remove invalid GPSRO rays from LETKF
! V1_31        2014-08-21 Andreas Rhodin
!  fix bugs ad consistency checks for fof-file input from COSMO.
!  letkf adaptive localisation: read 'rho_local.dat' from input/ directory.
! V1_35        2014-11-07 Andreas Rhodin
!  percent_fg_check<0: option to take fg check from deterministic run
! V1_37        2014-12-23 Andreas Rhodin
!  apply_H_member: make sure ph is present for ICON
! V1_42        2015-06-08 Andreas Rhodin
!  ICON coarse grid LETKF; MEC temporal interpolation; remove bad GPSRO rays
! V1_43        2015-08-19 Harald Anlauf
!  Add dealloc_fields feature to LETKF
! V1_44        2015-09-30 Andreas Rhodin
!  new option to handle missing members; fixes for MEC & GRIB2
! V1_45        2015-12-15 Andreas Rhodin
!  account for fillvalues; read_obs_ens: revised pressure level handling
! V1_46        2016-02-05 Andreas Rhodin
!  changes for COMETs (italian met. center) LETKF version
! V1_47        2016-06-06 Andreas Rhodin
!  read_det: new optional parameter gridfile,ierr,runtype,leveltypes
!  read_atm_ens: new optional parameter ierr
! V1_48        2016-10-06 Andreas Rhodin
!  apply_h: run observation operator on mean first
!  read_obs_ens: account for gross error check
!  modified mean/spread/talagrand-index calculation for DD
! V1_50        2017-01-09 Harald Anlauf
!  OpenMP optimization; fix for ICON ensemble MEC
! V1_51        2017-02-24 Andreas Rhodin
!  option to add ensemble perturbations from file (Valerio Cardinali, COMET)
!  changes in rppts versiond 1 and 2 (Valerio Cardinali, COMET)
!  LETKF: dont use spline interpolation for analysis in obs-space (overkill)
!  fix bug for 'percent_fg_check=-1.0' in namelist /fof_input/
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2005       First EnKF implementation
! Hendrik Reich   DWD  2007-2008  LETKF implementation
!==============================================================================
#include "tr15581.incf"
  !=============
  ! Modules used
  !=============
  use mo_kind,        only: wp,              &! working precision kind
                            sp                ! single  precision kind
  use mo_dace_string, only: char3,           &! conversion: integer -> char(3)
                            split,           &! char string -> array
                            eval_string       ! evaluate strings
  use mo_time,        only: t_time,          &! date/time derived type
                            operator (/=),   &! compare derived types t_time
                            cyyyymmddhhmm     ! derive string from time
  use mo_exception,   only: finish,          &! abort routine
                            message           ! issue a warning
  use mo_mpi_dace,    only: dace,            &! MPI group info
                            p_bcast,         &! generic MPI broadcast routine
                            p_send,          &! generic MPI send      routine
                            p_recv,          &! generic MPI receive   routine
                            p_min,           &! min over PEs
                            p_max,           &! max over PEs
                            p_sum,           &! sum over PEs
                            p_and,           &! and
                            p_or              ! or, ior
  use mo_obs_sndrcv,  only: p_bcast           ! bcast variable of type obs
  use mo_run_params,  only: input,           &! input  directory path
!                           output,          &! output directory path
                            obsinput,        &! observation input directory path
                            path_file,       &! concatenate pathname/basename
                            basename,        &! strip leading directories
                            nproc1,          &! partitioning parameter
                            nproc2,          &! partitioning parameter
!                           model,           &! 'COSMO' or 'GME'
                            method,          &! 'ENKF', 'LETKF', 'PSAS+LETKF'
                            fc_hours,        &! analysis cycle time
                            center,          &! generating center id
                            subcenter,       &! subcenter id
!                           proc_ana,        &! generating process
                            proc_ana_ens,    &! generating process
!                           ensemble_id,     &! ensemble id for GRIB output
                            flag_biasc_rad,  &! radiance bias correction flag
                            dealloc_fields    ! fields to deallocate
  use mo_p_output,    only: oline,           &! output line buffer
                            iol,             &! number of next line to write
                            nextline,        &! routine to increment line number
                            flush_buf         ! routine to write buffer
  use mo_t_datum,     only: t_datum,         &! observation          data type
                            rvind             ! invalid value
  use mo_t_obs,       only: t_obs,           &! observation 'box'    data type
                            t_spot,          &! observation 'report' data type
                            construct,       &! t_obs constructor routine
                            destruct,        &! t_obs destructor  routine
                            add_source,      &! extend observation file list
                            empty_source,    &! clear  observation file list
                            join_obs,        &! join observations from boxes
                            source,          &! list of source-files
                            ft_name,         &! get mnemonic of file type
                            release_mem,     &! release unused memory in t_obs
                            FT_FEEDBACK,     &! flag for feedback file
                            FT_MISSING,      &! flag for missing  file
                            FT_EMPTY,        &! flag for empty file
                            POC_FG
  use mo_obs_set,     only: t_obs_set,       &! observation data type
                            destruct          ! deallocate derived type compon.
  use mo_obs_tables,  only: decr_rpt_use,    &! change use-flags of report
                            write_pending,   &! write pending reports
                            derive_rept_stat,&! re-derive obs.statistics
                            gather_rept_stat,&! gather statistics on I/O PE
                            print_rept_stat   ! print observation statistics
  use mo_test_obs,    only: test_obs          ! test integrity of obs.data
  use mo_t_use,       only: CHK_CONSIST,     &! consistency     flag value
                            CHK_FG,          &! first guess     flag value
                            CHK_FG_LBC,      &! fg chk. vs. bc. flag value
                            CHK_NONE,        &! no check        flag value
                            CHK_BLACKLIST,   &! blacklist       flag value
                            CHK_HEIGHT,      &! height          flag value
                            CHK_GROSS,       &! gross error     flag value
                            CHK_THIN,        &! thinning        flag value
                            CHK_NOTUSED,     &! not used        flag value
                            CHK_RULE,        &! rule            flag value
                            CHK_BIASCOR,     &! biasc not poss. flag value
                            CHK_REDUNDANT,   &! redundant       flag value
                            CHK_CLOUD,       &! cloud           flag value
                            CHK_NO_OBS,      &! no obs          flag value
                            chk,             &! mnemonics, ... for 'checks'
                            n_chk,           &! max.value for CHK_...
                            STAT_DISMISS,    &! status          flag values
                            STAT_REJECTED,   &!
                            STAT_PAS_REJ,    &!
                            STAT_PASSIVE,    &!
                            STAT_ACTIVE       !
  use mo_boxes,       only: set_input_boxes   ! distribute files over PEs
  use mo_fdbk,        only: t_fdbk,          &! feedback file data type
                            open_fdbk_read,  &! open feedback file for reading
                            read_meta,       &! read verification meta data
                            print_fdbk,      &! print verification meta data
                            get_veri_index,  &! return indices of verif. runs
                            get_veri,        &! read verification entry
                            close_fdbk,      &! close a feedback file
                            cleanup_fdbk      ! deallocate t_fdbk components
  use mo_dec_matrix,  only: t_vector,        &! vector data type
                            construct,       &! construct vector or matrix
                            destruct,        &! destruct  vector or matrix
!                           global,          &! hold vector segments globally
                            local,           &! hold vector segments locally
                            operator  (+),   &! t_vector + t_vector
                            assignment(=),   &! t_vector = ...
                            size,            &!
                            minval,          &!
                            maxval,          &!
                            max,             &! maximum
                            count,           &!
                            operator  (-),   &! t_vector - t_vector
                            operator  (*),   &! t_vector * ...
                            operator  (**),  &! t_vector ** integer
                            operator  (/),   &! 'elemental' vector / vector
                            operator  (==)    ! t_vector == t_vector
  use mo_t_table,     only: name_value        ! find name     of table entry
  use mo_fdbk_tables, only: VT_FIRSTGUESS,   &! first guess flag value
                            VT_FORECAST,     &! forecast    flag value
                            VE_DETERM,       &! deterministic run flag value
                            OT_RAD,          &! RADiances observation type
                            OT_RADAR,        &! RADAR report id.
                            VN_DD,           &! variable Id for wind direction
                            VN_U,VN_V,       &!
                            VN_U10M,VN_V10M, &!
                            VN_RAWBT,        &! satellite brightness temperature
                            varno,           &! variable Id. table
                            obstype,         &! obstype  Id. table
                            OF_MISSING,      &! operator flag missing
                            OF_BT_CLEAR_SKY   ! clear-sky BT
  use mo_fdbk_in,     only: read_feedback,   &! read the feedback file
                            check_member,    &! check for correct member in fof
                            percent_fg_check,&! % of fg checks for rejection
                            read_nml_fof_input! read namelist /FOF_INPUT/
  use mo_fdbk_io,     only: get_var           ! read variable from fdbk-file
  use mo_tovs,        only: read_fdbk_rad     ! post-process fdbk-reading for TOVS
  use mo_grib_handling,only:t_inventory,     &! Inventory Data Type
                            get_inventory,   &! read inventory table
                            relabel_fields,  &! re-label fields from SST-ana
                            print_inventory, &! print inventory
                            ENS_MEAN,        &! ensemble mean
                            ENS_SPREAD,      &! ensemble spread
                            set_defaults,    &! set center, subcenter, process
                            set_dwd_def_ens   ! defaults for local part of PDB
  use mo_fortran_units,only:get_unit_number, &! request a free unit number
                            return_unit_number! return the unit number
  use mo_wmo_tables,  only: WMO6_ROTLL,      &! rotated lat/lon grid  type
                            WMO6_LATLON,     &! lat/lon         grid  type
!                           DWD6_ICON,       &! ICON            grid  type
                            DWD6_ICOSAHEDRON  ! GME             grid  type
  use mo_atm_state,   only: t_atm,           &! atmospheric state derived type
                            construct,       &! construct atmospheric state
                            destruct,        &! free atmospheric state
                            set_p,           &! set pressure levels
                            set_ph,          &! set pressure at half levels
                            set_ps,          &! set surface pressure (COSMO)
                            set_geo,         &! set geopotential
                            set_tv,          &! set virtual temperature
                            set_rh2m,        &! set 2m relative humidity
                            set_rh2m_land,   &! set 2m relative humidity (land)
                            C_to_A_grid,     &! transform from C- to A-grid
                            A_to_C_grid,     &! trafo from A-grid to C-grid
                            uvtrafo,         &! trafo of (u,v) geo<-->rot
                            assignment(=),   &! t_atm = ...
                            operator  (+),   &! t_atm + t_atm
                            operator  (-),   &! t_atm - t_atm
                            operator  (*),   &! t_atm * ...
                            operator  (/),   &! t_atm / ...
                            operator  (**),  &! t_atm ** ...
                            mean_stdev,      &! calculate mean and stdev
                            sqrt,            &! sqrt of t_atm
                            print,           &! print  atmospheric state
                            nm,              &! max. number of fields in t_atm
                            allocate,        &! allocate components of t_atm
                            deallocate,      &! deallocate components of t_atm
                            vert_intp,       &! vertical   interpolation
                            keep_diamond,    &! only keep one diamond (testing)
                            T,U,V,Q,W_SO      ! enumeration positions in t_atm
  use mo_atm_grid,    only: t_grid,          &! atmospheric grid derived type
                            construct,       &! construct grid
                            maxls,           &! max. number of soil layers
                            VCT_P_HYB,       &! GME/HRM Vertical coordinate
                            MO_ICON           ! ICON    model
  use mo_grib,        only: read,            &! read atmospheric state or grid
                            write_grib,      &! write atmospheric state
                            t_ctr_inv,       &! inventory container
                            IE_FILE           ! return code for file open error
  use mo_cpu_time,    only: stop_time         ! determine cpu and wall time
  use mo_grid_intpol, only: grid_indices,    &! indices + interp.coefficients
                            mp                ! max # of coefficients
  use mo_grads                                ! grads plots
  use mo_postmult,    only: setup_post,      &! set up model columns (t_cols)
                            finish_post,     &! finish postm.: t_cols to t_atm
                            finish_ana        ! finish analysis
  use mo_t_bg_err_op, only: covm,            &! covariance matrix meta info
                            random_gauss      ! Gaussian random numbers
  use mo_bg_err_2d,   only: apply_L_m         ! multiply by L (sqrt B)
  use mo_t_col,       only: t_cols,          &! model column data type
                            dealloc_cols,    &! deallocate model columns
                            get_cols,        &! redistribute model columns
                            set_cols_logp     ! store ln p instead of p
  use mo_obs,         only: process_obs       ! generic observation processing
  use mo_t_obs,       only: TSK_SETUP_COLS,  &! flags passed to process_obs
                            TSK_SETUP_FULL,  &! set up interpolation space
                            TSK_Y,           &!    evaluate operator
                            TSK_K             !    evaluate operator + adjoint
  use mo_fg_checks,   only: check_operator    ! check applicability of obsop.
  use mo_psasutil,    only: interpolate,     &! apply H operator
                            interpolate_strat
! use meteo_utilities,only: calps             ! compute surf.pres
  use utilities,      only: uvrot2uv
  use mo_physics,     only: pi,              &! 3.1415...
                            d2r, r2d,        &! convert degree <-> radians
                            t0c,             &! T [K] for T=0 degree Celsius
                            R_c => R,        &! gas constant of dry air
                            c_p,             &! specific heat [J/(kg*K)]
                            esw_t,           &! sat. vap. press.<-temperature
                            q_e,             &! spec.humid.<-water vap.press.
                            rh_q              ! rel.humidity <- spec.hum.
  use mo_radar,       only: thin_radar,      &! thinning of radar observations
                            shrink_par_radar  ! adapt table 'par'
  use mo_rad_enkf,    only: plevel_select,   &! how to determine plevel from ens
                            letkf_obserr_ir, &! obserr for different IR instruments
                            n_letkf_obserr_rad! # of entries in letkf_obserr_ir
  use mo_profile,     only: write_profiles    ! write atmospheric profiles
  use mo_letkf_util,  only: flat_copy         ! flat copy of atmospheric state
  use mo_soil,        only: n_sl,            &! number of soil layers
                            n_st,            &! number of soil types
                            soil_thickness,  &! soil thickness (m)
                            st,              &! soil thickness (mm)
                            Cadp,            &! air dryness point   (lower bound)
                            Cpwp,            &! plant wilting point (lower bound)
                            Cfcap,           &! field capacity      (upper bound)
                            Cporv,           &! pore volume         (upper bound)
                            ST_ICE, ST_ROCK, &!
                            ST_SEAWATER, ST_SEAICE
  use netcdf,         only: nf90_get_var,    &! read NetCDF file directly
                            nf90_inq_varid,  &! .. currently used for
                            nf90_noerr,      &!    weighting functions
                            nf90_nowrite,    &!
                            nf90_open,       &!
                            nf90_close        !
  use mo_wmoBufrIds,  only: SEVIRI            ! instrument ID for SEVIRI
  use mo_instrid,     only: hss_instr         ! to check whether an instrument is a hyperspectral sounder

  implicit none

  !================
  ! Public entities
  !================
  private
  public :: read_obs_ens     ! read ensemble observations (fof-files from COSMO)
  public :: read_det         ! read deterministic forecast or analysis
  public :: read_ens         ! read ensemble state
  public :: derive_obs_ana   ! derive estimate of analysis in observation space
  public :: rms_statistics   ! compute rms/spread
  public :: read_analysis    ! read (oper.) analysis file
  public :: compute_mean     ! compute mean and spread of ensemble
  public :: subtract_mean    ! subtract given mean from ensemble
  public :: add_mean         ! add given mean to ensemble
  public :: add_noise        ! add climatological noise to ensemble
  public :: relax_prior_pert ! relaxation to prior perturbation
  public :: relax_prior_spread ! relaxation to prior spread
  public :: apply_H          ! apply obs operator to ensemble
  public :: apply_H_m        ! apply obs operator to single state
  public :: scatter_coarse   ! scatter fields allocated on the 'coarse' grid
  public :: gather_coarse    ! gather  fields allocated on the 'coarse' grid
  public :: read_infl_factor ! read adaptive inflation factor
  public :: mean_spread_talag! compute mean,spread + talag vector (obs space)
  public :: sat_adj          ! saturation adjustment
  public :: apply_hum_bounds ! apply humidity bounds
  public :: apply_tke_bounds ! apply lower bound to TKE
  public :: set_d_ps         ! auxiliary routines for localization ...
  public :: apply_d_ps       ! ... on hybrid pressure grid (e.g., GME)
  public :: to_common_grid   ! interpolate to common vertical grid
  public :: from_common_grid ! interpolate to common vertical grid
  public :: adjust_soil      ! adjust ensemble soil moisture (w_so) to analysis
  public :: inflate_wso_spread ! inflate soil moisture spread

  !-----------
  ! Interfaces
  !-----------
  interface apply_H
    module procedure apply_H   ! apply observation-operator to ensemble
    module procedure apply_H2  ! apply operator with temporal interpolation
    module procedure apply_H_m ! apply operator to single state (with interp.)
  end interface

!==============================================================================
contains
!==============================================================================
  subroutine read_obs_ens (obs_out, H_x, k_enkf, prefix, l_det, H_det, &
                           cl_impact, members, num_vert_lev, weighting_func, &
                           ens_cld_fl)
  !============================================================
  ! read ensemble observations (from COSMO) from feedback files
  !============================================================
  type(t_obs_set)  ,intent(inout) :: obs_out    ! observation derived type
  type(t_vector)   ,pointer       :: H_x (:)    ! observation operator applied
  integer          ,intent(in)    :: k_enkf     ! ensemble size
  character(len=*) ,intent(in)    :: prefix(:)  ! file prefix
  logical          ,intent(in)    :: l_det      ! flag to read deterministic fc
  type(t_vector)   ,intent(out)   :: H_det      ! obs.op. applied to det. fc
  type(t_vector)   ,intent(out)   :: cl_impact  !cloud impact for error model
  type(t_vector)   ,intent(out)   :: ens_cld_fl ! everage of ensemble cloud flag
  integer, optional ,intent(in)   :: members(:) ! members to read
  integer, optional ,intent(in)   :: num_vert_lev   ! number of vertical levels
  logical, optional ,intent(in)   :: weighting_func ! use weighting functions for radiance localization

    character(len=128)          :: infile     ! input file name (date replaced)
    integer                     :: i, j, l
    integer                     :: ftype
    type(t_obs)                 :: o
    type(t_obs_set)             :: obs_in     ! observation container for input
    type(t_obs_set),allocatable :: obs_tmp(:) ! temporary observation container
    type(t_obs_set)             :: obs_src    ! temporary observation container
    type(t_vector) ,allocatable :: Hx_tmp(:,:)! temporary fc ensemble container
    type(t_vector) ,allocatable :: Hd_tmp  (:)! temporary determin.fc container
    type(t_vector) ,allocatable :: Hx_tmp_clear(:,:)! temporary fc ens container
    type(t_vector), allocatable :: cloud_impact(:)  ! temporary container
    type(t_vector), allocatable :: cloud_fl_mem(:) ! average of ensemble members' cloudiness
    type(t_vector), allocatable :: cloud_fl_tmp ! temporary cloud flag container
    real(wp)       ,allocatable :: pl_in   (:)! pressure levels (input)
    real(wp)       ,allocatable :: pl_tmp  (:)! pressure levels (re-ordered)
    real(wp)       ,allocatable :: pl_mean (:)! pressure levels (ensemble mean)
    logical        ,allocatable :: pl_valid(:)! valid pressure levels
    real(wp)       ,allocatable :: pl_high(:) ! highest plevel of members
    type(t_vector)              :: H_0        ! ensemble mean
    integer                     :: nfof       ! number of files to read
    type (t_fdbk)               :: fb         ! feedback file content
    integer                     :: ifg (1)    ! first guess index
    integer                     :: is         ! observation 'spot' index
    integer                     :: nidx       ! number of indices returned
    integer                     :: ifg_clear(1) !first guess index clear-sky TB
    integer                     :: nidx_clear ! number of indices
    integer                     :: i1, i2, i2x! header loop index
    integer                     :: j1, j2, j2x! body loop index
    integer                     :: k1, k2, l1, l2
    logical        ,allocatable :: mh (:)     ! mask for valid header entries
    logical        ,allocatable :: mb (:)     ! mask for valid body   entries
    logical        ,allocatable :: mbi(:)     ! mask for valid body   entries
    integer        ,allocatable :: mi (:)     ! body missing in member 'mi'
    integer                     :: nh, nb     ! number of valid entries
    integer        ,pointer     :: tmp(:)
    integer                     :: n_rep      ! number of reports
    integer                     :: n_obs      ! number of observations
    integer                     :: n          ! number of files to read
    integer                     :: ifi        ! file index to read / member
    integer                     :: if1        ! 1st file with content
    integer                     :: k          ! counter
    integer                     :: m          ! member index
    integer                     :: ens_member ! member# or determ.fc
    integer                     :: itrim      ! index of last char in file name
    logical                     :: fof_optional! true for optional files
    logical                     :: fof_ver    ! true for files from MEC   (ver...nc)
    logical                     :: fof_ekf    ! true for files from LETKF (ekf...nc)
    integer                     :: p_1        ! PE for 1st input file
    integer        ,allocatable :: s_state (:)! spot% state * 1000 + i_ens
    integer        ,allocatable :: b_state (:)! body% state * 1000 + i_ens
    integer        ,allocatable :: s_check (:)!
    integer        ,allocatable :: b_check (:)!
    integer        ,allocatable :: b_fg    (:)! # of members with fg check active
    integer                     :: flags
    integer                     :: n_fg       ! # of members required for fg check
    logical                     :: luvbug     ! flag to fix u,v rotated to model grid
    real(wp)                    :: ugeo, vgeo ! u,v rotated back to geographic coordinates
    logical                     :: query_member!check member when reading input?
    integer                     :: stride     ! PE stride factor
    integer        ,allocatable :: pio   (:)  ! I/O PE of members
    character(128) ,allocatable :: fname (:)  ! member file name
    logical        ,allocatable :: exists(:)  ! file exists?
    !----------------------------------
    ! variables for weighting functions
    !----------------------------------
    real(wp)       ,allocatable :: w_f_in    (:,:)
    real(wp)       ,allocatable :: w_f_mean  (:,:)
    real(wp)       ,allocatable :: p_wf_in   (:,:)
    real(wp)       ,allocatable :: p_wf_mean (:,:)
    integer                     :: status,ncid,varid_w_f,lll,varid_p_wf
    logical                     :: lrad
    logical                     :: lw_f
    integer                     :: vert_lev
    integer                     :: num_obs
    integer                     :: num_head
    integer                     :: runtype

    type t_pi
      integer,  pointer :: ih (:)               ! header index array
      integer,  pointer :: ib (:)               ! body   index array
      real(wp), pointer :: Hx (:)               ! first guess
      logical,  pointer :: mp (:)               ! mask of processed blocks
    end type t_pi
    type (t_pi)         :: ix (k_enkf+1)        ! index array
    type(t_spot)  ,pointer :: si       ! observation 'spot' pointer

    !----------------------------------------------------------------
    ! settings for input of vertical localisation weighting functions
    !----------------------------------------------------------------
    lrad = .false.
    lw_f = .false.
    vert_lev = 0
    if (present (weighting_func))  lw_f     = weighting_func
    if (present (num_vert_lev  ))  vert_lev = num_vert_lev

    !-----------------------------------------------------------
    ! determine array sizes, allocate temporaries, read namelist
    !-----------------------------------------------------------
    if1  = 0                                  ! first non-empty file
    n    = k_enkf; if (l_det) n = k_enkf + 1  ! no. ensembles + det.run
    nfof = count (prefix /= '')               ! no. files per member
!   if (nfof == 0) call finish ('read_obs_ens','nfof == 0')
!   if (nfof == 0) goto 999
    allocate (Hd_tmp     (nfof       ))
    allocate (Hx_tmp     (nfof,k_enkf))
    allocate (Hx_tmp_clear(nfof,k_enkf))
    allocate (obs_tmp    (nfof       ))
    allocate (cloud_impact(nfof      ))
    allocate (cloud_fl_mem(nfof))
    allocate (cloud_fl_tmp)

    do ifi = 1, nfof
      allocate       (obs_tmp(ifi)% o(1))
      allocate       (obs_tmp(ifi)% oi  )
      call construct (obs_tmp(ifi)% o   )
    end do
    call read_nml_fof_input  ! namelist /FOF_INPUT/

    !------------------------------------------------------------------
    ! Distribution of I/O (try stride > 1 for better I/O bandwidth use)
    !------------------------------------------------------------------
    stride = max (dace% npe / n, 1)
    allocate (pio(n))
    do i = 1, n
       pio(i) = mod ((i-1) * stride, dace% npe)
    end do
    allocate (fname (n))
    allocate (exists(n))

    !-----------------------------
    ! set up list of files to read
    !-----------------------------
    do ifi = 1, nfof
      call empty_source
      k            = 0
      query_member = check_member   ! preset with namelist value check_member
      fof_optional = index (prefix(ifi), '#OPTIONAL#') /= 0
      itrim        = index (prefix(ifi), ' ') - 1
      fof_ver      = index (prefix(ifi), 'ver')        == 1 &
                .or. index (prefix(ifi), 'mof')        == 1
      fof_ekf      = index (prefix(ifi), 'ekf')        == 1
      do i = 1, n
        if (fof_ver .or. fof_ekf) then
          !--------------------------------------
          ! one file for ver*.nc or ekf*.nc files
          !--------------------------------------
          infile = path_file ('', prefix(ifi)(1:itrim)//'.nc')
        else
          !-----------------------------------
          ! one file per member for fof* files
          !-----------------------------------
          if (i==1 .and. l_det) then
            infile = path_file ('', prefix(ifi)(1:itrim)//'__FCR_TIMEMMSS_.nc')
          else
            k = k + 1
            m = k; if (present (members)) m = members (k)
            infile = path_file ('', prefix(ifi)(1:itrim)//'__FCR_TIMEMMSS__ens',&
                                    char3(m)//'.nc'                             )
          endif
        endif
        fname(i) = infile
      end do
      !----------------------------------------------
      ! distributed check for presence of input files
      !----------------------------------------------
      do i = 1, n
        if (dace% pe == pio(i)) then
          inquire (file=path_file (obsinput, fname(i)), exist=exists(i))
        end if
      end do
      do i = 1, n
        call p_bcast (exists(i), pio(i))
        ftype = FT_MISSING; if (exists(i)) ftype = FT_FEEDBACK
        call add_source (obsinput, fname(i), ftype, complete=.true.)
      end do
      !--------------------------
      ! distribute files over PEs
      !--------------------------
      if(dace% lpio) write(6,'()')
      do i=1,n
        source(i)% pe = pio(i)
        if(dace% lpio) then
          write(6,'(i4,1x,a64,1x,a,i4)') i, source(i)% file, &
            ft_name(source(i)% filetype), source(i)% pe
        endif
      end do
      if(dace% lpio) write(6,'()')
      !-------------------
      ! skip missing files
      !-------------------
      if (any (source(1:n)% filetype == FT_MISSING))then
        if(dace% lpio) then
          if (fof_optional) then
             write(6,'(a/)') 'SKIPPING optional missing file !'
          else
             write(6,'(a/)') 'ABORTING: non-optional missing file !'
             call finish ("read_obs_ens",                           &
                          "file not found: "//trim(source(ifi)%file))
          endif
        endif
        cycle
      endif

      !------------------
      ! scan file content
      !------------------
!     allocate (obs_in% o (n))
!     call construct (obs_in% o)
!     call process_obs (TSK_INIT, obs_in, fc)
      call set_input_boxes

      !-----------------
      ! skip empty files
      !-----------------
      if (any (source(1:n)% filetype == FT_EMPTY))then
        if(dace% lpio) write(6,'(a/)') 'SKIPPING empty file !'
!       call destruct (obs_in% o)
!       deallocate    (obs_in% o)
        cycle
      endif

      allocate (obs_in% o (n))
      call construct (obs_in% o)

      !------------------
      ! read observations
      !------------------
      call read_feedback (pass=2, obs=obs_in% o,  &
                            id=.false.,           &
                           chk=.true.,            &
                           src=ifi,               &
!                          spec=flag_biasc_rad >= 0)
                          spec=.true.)
      call write_pending
!     call unique_spot (obs_in% o)
      call release_mem (obs_in% o, sort=.false.)
      !-------------------------------------------
      ! broadcast first ensemble member to all PEs
      !-------------------------------------------
      p_1 = source(1)% pe
      if (p_1 == dace% pe) o = obs_in% o(1)
      call p_bcast (o, p_1)

      !-------------------------------
      ! check consistency of ensembles
      !-------------------------------
      if (dace% lpio) then
        write (6,'()')
        write (6,'(a)') &
          ' check consistency of observations for ensemble members'
        write (6,'()')
      endif
      !-------------------------------------------------------------
      ! set up index arrays and masks for reordering of observations
      !-------------------------------------------------------------
      allocate (mh (o% n_spot)); mh = .true.
      allocate (mb (o% n_obs )); mb = .true.
      allocate (mi (o% n_spot)); mi = huge(mi)
      !------------------
      ! loop over members
      !------------------
      do i=2,n
        if (source(i)% pe == dace% pe) then
          allocate (ix(i)% ih (        o   % n_spot)); ix(i)% ih = 0
          allocate (ix(i)% ib (        o   % n_obs )); ix(i)% ib = 0
          allocate (ix(i)% Hx (obs_in% o(i)% n_obs ))
          allocate (ix(i)% mp (obs_in% o(i)% n_spot)); ix(i)% mp = .false.
          !-------------
          ! check header
          !-------------
          i1 = 1
          i2 = 1
          do
            if (i1 > o% n_spot) exit
            if (.not. mh(i1) .and. i2 <= obs_in% o(i)% n_spot) then
              if (same_head (i, i1, i2)) then
                ix(i)% mp (i2) = .true.
                i2 = i2 + 1
              endif
              i1 = i1 + 1
              cycle
            endif
!           if (i2 > obs_in% o(i)% n_spot) then
!             mh (i1:) = .false.
!             mi (i1:) = min (mi (i1:), i)
!             exit
!           endif
            if (i2 <= obs_in% o(i)% n_spot) then
              if (same_head (i, i1, i2)) then
                ix(i)% mp (i2) = .true.
                ix(i)% ih (i1) = i2
                i1 = i1 + 1
                i2 = i2 + 1
                cycle
              endif
            endif
!           i2x = i2 + 1
            i2x = 1
            do
              if (i2x > obs_in% o(i)% n_spot) then
                mh (i1) = .false.
                mi (i1) = min (mi (i1), i)
                i1 = i1 + 1
                exit
              endif
              if (ix(i)% mp (i2x)) then
                i2x = i2x + 1
                cycle
              endif
              if (same_head (i, i1, i2x)) then
                ix(i)% mp (i2x) = .true.
                ix(i)% ih (i1)  = i2x
                i1 = i1  + 1
                i2 = i2x + 1
                exit
              else
                i2x = i2x + 1
              endif
            end do
          end do
          !-----------
          ! check body
          !-----------
          do i1 = 1, o% n_spot
            k1 = o % spot(i1)% o% i + 1
            k2 = o % spot(i1)% o% i + o % spot(i1)% o% n
            j1 = k1

            if (.not. mh(i1)) then
              mb (k1:k2) = .false.
              cycle
            endif

            i2 = ix(i)% ih (i1)
            l1 = obs_in% o(i)% spot(i2)% o% i + 1
            l2 = obs_in% o(i)% spot(i2)% o% i + obs_in% o(i)% spot(i2)% o% n
            j2 = l1
            do
              if (j1 > k2) exit
              if (.not. mb(j1)) then
                j1 = j1 + 1
                cycle
              endif
              if (j2 > l2) then
                mb (j1:k2) = .false.
                exit
              endif
              if (same_body (i, i1, j1, j2)) then
                ix(i)% ib (j1) = j2
                j1 = j1 + 1
                j2 = j2 + 1
                cycle
              endif
              j2x = j2 + 1
              do
                if (j2x > l2) then
                  mb (j1) = .false.
                  j1 = j1 + 1
                  exit
                endif
                if (same_body (i, i1, j1, j2x)) then
                  ix(i)% ib (j1) = j2x
                  j1 = j1  + 1
                  j2 = j2x + 1
                  exit
                else
                  j2x = j2x + 1
                endif
              end do
            end do
            if (.not. any(mb (k1:k2))) mh (i1) = .false.
          end do
          deallocate (ix(i)% mp)
        endif
      end do

      mh = p_and (mh)
      mb = p_and (mb)
      mi = p_min (mi); if (l_det) mi = mi - 1
      !---------------------------------
      ! printout of non-matching entries
      !---------------------------------
      if (dace% lpio) then
        nh = count(mh)
        nb = count(mb)
        write(6,*)
        write(6,*) '  header entries in member_1 not matched:', &
                      (size(mh)-nh), ' / ', size(mh)
        write(6,*) '  body   entries in member_1 not matched:', &
                      (size(mb)-nb), ' / ', size(mb)
        write(6,*) '  header entries remaining              :', &
                      nh
        write(6,*) '  body   entries remaining              :', &
                      nb
        write(6,*)
        write(6,*) '  first 200 header entries not matched:'
        write(6,*)
        l = 0
        do i = 1, o% n_spot
          if (.not. mh (i)) then
            call decr_rpt_use (o% spot(i),  CHK_CONSIST,                     &
                                            STAT_DISMISS, 'ensemble mismatch')
            l = l + 1
            if (l > 200) cycle
            write(6,'(i10,1x,a12,a,2f8.2,i4)')                     &
              i, name_value (obstype, o% spot (i)% hd%   obstype), &
                                      o% spot (i)%       statid,   &
                                      o% spot (i)%col%c% dlat,     &
                                      o% spot (i)%col%c% dlon,     &
                                           mi (i)
          endif
        end do
        write(6,*)
        write(6,*) '  first 1000 body entries not matched:'
        write(6,*)
        l = 0
loop:   do i = 1, o% n_spot
          if (.not. mh (i)) then
            write(6,'(i10,1x,a12,a,2f10.2)')                     &
              i, name_value (obstype, o% spot (i)% hd% obstype), &
                                      o% spot (i)%       statid, &
                                      o% spot (i)%col%c% dlat,   &
                                      o% spot (i)%col%c% dlon
          endif
          k1 = o % spot(i)% o% i + 1
          k2 = o % spot(i)% o% i + o % spot(i)% o% n
          do l1 = k1, k2
            if (.not. mb (l1)) then
              l = l + 1
              if (l > 1000) exit loop
              write(6,'(i10,1x,a12,a,3f10.2,2i10,1x,a)')                      &
                 i,           name_value (obstype, o% spot (i)% hd% obstype), &
                                                   o% spot (i)%       statid, &
                                                   o% spot (i)%col%c% dlat,   &
                                                   o% spot (i)%col%c% dlon,   &
                                                   o% olev (l1),              &
                 l1-k1+1, l1, name_value (varno,   o% varno(l1))
            endif
          end do
        end do loop
        !--------------------------------------
        ! early thinning for specific operators
        !--------------------------------------
        call thin_radar  (o)
        mb = mb .and. o% body% use% state > STAT_DISMISS
        mh = mh .and. o% spot% use% state > STAT_DISMISS
        do i = 1, o% n_spot
          k1 = o % spot(i)% o% i + 1
          k2 = o % spot(i)% o% i + o % spot(i)% o% n
          if (.not. any(mb (k1:k2))) mh (i) = .false.
          if (.not. mh (i))      mb (k1:k2) = .false.
        end do
      endif
      if (p_1 /= dace% pio) then
        if (dace% lpio) then
          call p_send (o% spot% use% state            ,p_1  ,1)
          call p_send (o% body% use% state            ,p_1  ,2)
        endif
        if (dace% pe == p_1 ) then
          call p_recv (obs_in% o(1)% spot% use% state ,dace% pio ,1)
          call p_recv (obs_in% o(1)% body% use% state ,dace% pio ,2)
        endif
      endif

      !----------------------------------------------------------
      ! communicate number of reports/observations after thinning
      !----------------------------------------------------------
      mh    = p_and (mh)
      mb    = p_and (mb)
      n_rep = count (mh)
      n_obs = count (mb)
      if(dace% lpio) then
        write(6,*) '  header entries after first checks     :', &
                      n_rep, ' / ' ,nh
        write(6,*) '  body   entries after first checks     :', &
                      n_obs, ' / ' ,nb
        write(6,*)
      endif
      if (if1 == 0 .and. n_obs > 1) if1 = ifi
      call write_pending

      !------------------------------------
      ! reorganise observations
      ! enforce consistence for all members
      !------------------------------------
      call flush_buf
      if (p_1 == dace% pe) then
        allocate (ix(1)% ih (o% n_spot)); ix(1)% ih = (/(i,i=1,o% n_spot)/)
        allocate (ix(1)% ib (o% n_obs )); ix(1)% ib = (/(i,i=1,o% n_obs )/)
        allocate (ix(1)% Hx (obs_in% o(1)% n_obs))
      endif
      do i = 1, n
        if (source(i)% pe == dace% pe) then

          allocate (tmp (n_rep))
          tmp = pack (ix(i)% ih, mh)
          deallocate (ix(i)% ih)
          ix(i)% ih => tmp

          allocate (tmp (n_obs))
          tmp = pack (ix(i)% ib, mb)
          deallocate (ix(i)% ib)
          ix(i)% ib => tmp

          allocate (mbi (obs_in% o(i)% n_obs))
          mbi             = .false.
          mbi (ix(i)% ib) = .true.
          call shrink_par_radar (obs_in% o(i),                        mbi)
          call pack_header      (obs_in% o(i)% spot, ix(i)% ih, ix(i)% ib, i)
          obs_in% o(i)% n_spot = n_rep
          deallocate (mbi)

          obs_in% o(i)% n_obs = n_obs
          call pack_body (obs_in% o(i)% body,  ix(i)% ib)
          call pack_int  (obs_in% o(i)% varno, ix(i)% ib)
          call pack_real (obs_in% o(i)% bger,  ix(i)% ib)
          call pack_real (obs_in% o(i)% s_vqc, ix(i)% ib)
          call pack_int  (obs_in% o(i)% f_vqc, ix(i)% ib)
          call pack_real (obs_in% o(i)% olev,  ix(i)% ib)
          call pack_par  (obs_in% o(i)                  )

          if (obs_in% o(i)% n_lev > 0) call finish ('read_obs_ens','n_lev > 0')
          if (obs_in% o(i)% n_int > 0) call finish ('read_obs_ens','n_int > 0')

        endif
      end do
      call flush_buf

      !------------------------
      ! cross check consistence
      !------------------------
      if (p_1 == dace% pe) then
        o = obs_in% o(1)
      else
        call destruct (o)
      endif
      call p_bcast (o, p_1)
      do i=2,n
        if (source(i)% pe == dace% pe) then
          do i1 = 1, n_rep
            if (.not. same_head (i, i1, i1))                          &
              call finish ('read_obs_ens','inconsistent report header')
            do j1 = o% spot(i1)% o% i + 1 ,&
                    o% spot(i1)% o% i + o% spot(i1)% o% n
              if (.not. same_body (i, i1, j1, j1))                    &
                call finish ('read_obs_ens','inconsistent report body')
            end do
          end do
        endif
      end do

      !-----------------------------
      ! make status flags consistent
      !-----------------------------
      allocate (s_state (o% n_spot))
      allocate (s_check (o% n_spot))
      allocate (b_state (o% n_obs ))
      allocate (b_check (o% n_obs ))
      allocate (b_fg    (o% n_obs ))

      if (n > 999) call finish ('read_obs_ens','n > 999')

      !-------------------------
      ! handle first guess check
      !-------------------------
      s_state = o% spot% use% state * 1000 + 1   ! lowest state per report
      b_state = o% body% use% state * 1000 + 1   ! lowest state per observation
      b_fg    = 0                                ! # of fg-rejections / obs.
      n_fg    = floor (k_enkf * percent_fg_check /100._wp) ! # of fg-rej. rewquired

      if(dace% lpio) then
        write(6,'(a,i3,a,f4.0,a)') '  no.of hits required for fg-check      :', &
                                      n_fg, ' (',percent_fg_check,' % )'
        write(6,'()')
      endif

      !--------------------------------------------------------------
      ! mark observations fg-rejected in deterministic run for n_fg<0
      !--------------------------------------------------------------
      if ( (.not.l_det .or. n_fg<0)      &
          .and. source(1)% pe == dace% pe) then
        where (btest (o% body% use% flags, CHK_FG)    &
          .or. btest (o% body% use% flags, CHK_FG_LBC)) b_fg = b_fg + 1
      endif
      !-----------------------------------------------
      ! mark observations fg-rejected in ensemble runs
      !-----------------------------------------------
      do i=2,n
        if (source(i)% pe == dace% pe) then
          s_state = min (s_state, obs_in% o(i)% spot% use% state * 1000 + i)
          b_state = min (b_state, obs_in% o(i)% body% use% state * 1000 + i)
          if (n_fg >= 0) then
            where (btest (obs_in% o(i)% body% use% flags, CHK_FG)      &
              .or. btest (obs_in% o(i)% body% use% flags, CHK_FG_LBC)) &
              b_fg = b_fg + 1
          endif
          o% spot% use% flags = ior (o%            spot% use% flags, &
                                     obs_in% o(i)% spot% use% flags  )
          o% body% use% flags = ior (o%            body% use% flags, &
                                     obs_in% o(i)% body% use% flags  )
        endif
      end do
      s_state             = p_min (s_state)
      b_state             = p_min (b_state)
      b_fg                = p_sum (b_fg   )
      o% spot% use% flags = p_or  (o% spot% use% flags)
      o% body% use% flags = p_or  (o% body% use% flags)

      s_check = 0
      b_check = 0
      do i=1,n
        if (source(i)% pe == dace% pe) then
          where (mod (s_state, 1000) == i) &
            s_check = obs_in% o(i)% spot% use% check
          where (mod (b_state, 1000) == i) &
            b_check = obs_in% o(i)% body% use% check
        endif
      end do
      o% spot% use% check = p_sum (s_check)
      o% body% use% check = p_sum (b_check)

      o% spot% use% state = s_state / 1000
      o% body% use% state = b_state / 1000

      !--------------------------------------------
      ! clear first guess flag,
      ! set state consistently with remaining flags
      !--------------------------------------------
      do i=1,o% n_obs
        if    (btest (o% body(i)% use% flags, CHK_FG)   &              ! obs is fg rejected
          .or. btest (o% body(i)% use% flags, CHK_FG_LBC)) then
         if ((n_fg >= 0 .and. b_fg(i) < n_fg) .or. &                    ! by few members
             (n_fg <  0 .and. b_fg(i) == 0  )      ) then               ! not by det. run
          o% body(i)% use% flags = ibclr (o% body(i)% use% flags, CHK_FG)
          o% body(i)% use% flags = ibclr (o% body(i)% use% flags, CHK_FG_LBC)

          if (o% body(i)% use% check == CHK_FG .or. &
              o% body(i)% use% check == CHK_FG_LBC  ) then
            flags = o% body(i)% use% flags
            if (flags == 0) then
              o% body(i)% use% check = CHK_NONE
              select case (o% body(i)% use% state)
              case (STAT_PAS_REJ)
                o% body(i)% use% state = STAT_PASSIVE
              case (STAT_REJECTED)
                o% body(i)% use% state = STAT_ACTIVE
              end select
            else
              if (btest (flags, CHK_REDUNDANT)) then
                o% body(i)% use% check = CHK_REDUNDANT
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_REDUNDANT)
              endif
              if (btest (flags, CHK_BIASCOR)) then
                o% body(i)% use% check = CHK_BIASCOR
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_BIASCOR)
              endif
              if (btest (flags, CHK_RULE)) then
                o% body(i)% use% check = CHK_RULE
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_RULE)
              endif
              if (btest (flags, CHK_NOTUSED)) then
                o% body(i)% use% check = CHK_NOTUSED
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_NOTUSED)
              endif
              if (btest (flags, CHK_THIN)) then
                o% body(i)% use% check = CHK_THIN
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_THIN)
              endif
              if (btest (flags, CHK_GROSS)) then
                o% body(i)% use% check = CHK_GROSS
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_GROSS)
              endif
              if (btest (flags, CHK_CONSIST)) then
                o% body(i)% use% check = CHK_CONSIST
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_CONSIST)
              endif
              if (btest (flags, CHK_BLACKLIST)) then
                o% body(i)% use% check = CHK_BLACKLIST
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_BLACKLIST)
              endif
              if (btest (flags, CHK_HEIGHT)) then
                o% body(i)% use% check = CHK_HEIGHT
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_HEIGHT)
              endif
              if (btest (flags, CHK_NO_OBS)) then
                o% body(i)% use% check = CHK_NO_OBS
!               o% body(i)% use% state = STAT_REJECTED
                flags = ibclr (flags, CHK_NO_OBS)
              endif
              if (flags /= 0) then
                if(dace% lpio) then
                  write(0,'(a,b33)') 'flags =', flags
                  do i1=0,31
                    if (btest (flags, i1)) then
                     if (i1 >=1 .and. i1 <= n_chk) then
                      write(0,*) 'check flag set, bit=', i1, trim (chk(i1)% mnem)
                     else
                      write(0,*) 'check flag set, bit=', i1
                     end if
                    end if
                  end do
                end if
                call finish('read_obs_ens','unknown flags')
              endif
            endif
          endif
         endif
        endif
      end do

      call read_fdbk_rad('post',obs_in% o)

      !-------------------------------------------
      ! distribute observations to PEs
      ! (currently keep observations on PE 0 only)
      !-------------------------------------------
      obs_tmp(ifi)% o     = o
      obs_tmp(ifi)% o% pe = 0
      call construct (obs_tmp(ifi)% oi, o%n_obs, 1, nb=(/o%n_obs/), pe=(/0/))

      !-------------------------------------------------
      ! sum up cloud flagging result of ensemble members
      !-------------------------------------------------
      call construct (cloud_fl_mem(ifi), obs_tmp(ifi)% oi, global=.true.)
      call construct (cloud_fl_tmp,      obs_tmp(ifi)% oi, global=.true.)

      cloud_fl_mem(ifi)% s(1)% x = 0  !initialisation
      cloud_fl_tmp% s(1)% x   = 0
      do i=1,n
        if (source(i)% pe == dace% pe) then

          do is = 1, obs_in% o(i)% n_spot
            si => obs_in% o(i)% spot(is)
            if (si% hd% obstype /= OT_RAD) cycle   ! Skip non-radiance obstypes
            do i1 = si% o% i + 1, si% o% i + si% o% n
                if (obs_in% o(i)% varno(i1) /= VN_RAWBT) cycle
                if (hss_instr(INT(obs_in% o(i)% body(i1) %lev_sig), wmo = .true. )) then
                  ! only add one if cloud flag is set
                  if (btest (obs_in% o(i)% body(i1)% use% flags,CHK_CLOUD)) then
                    cloud_fl_tmp% s(1)% x(i1)  = cloud_fl_tmp% s(1)% x(i1) + 1
                  end if
                end if
            end do ! end do i1
          end do ! end do is

        end if ! end if source(i)% pe

      end do ! end do n

      cloud_fl_mem(ifi)% s(1)% x = p_sum (cloud_fl_tmp% s(1)% x)
      cloud_fl_mem(ifi)% s(1)% x = cloud_fl_mem(ifi)% s(1)% x / n

      !--------------------------------------------------
      ! only applicable to hyperspectral radiance
      ! set obs state to passive when cloud limit exceeds
      !--------------------------------------------------
      do is = 1, o% n_spot
        si => o% spot(is)
        if (si% hd% obstype /= OT_RAD) cycle   ! Skip non-radiance obstypes
        do i1 = si% o% i + 1, si% o% i + si% o% n
          if (o% varno(i1) /= VN_RAWBT) cycle
          if (hss_instr(INT(o% body(i1) %lev_sig), wmo = .true. )) then
            do i = 1, n_letkf_obserr_rad
              if ( letkf_obserr_ir(i)% rad_error_model > 0) then
                if (si% sttyp == letkf_obserr_ir(i)% instr) then
                  if (cloud_fl_mem(ifi)% s(1)% x(i1) >= letkf_obserr_ir(i)% ens_cld_limit ) then
                    o% body(i1)% use% state = STAT_PASSIVE
                  end if ! if cloud_fl_mem > ens_cld_limit
                end if ! if si% sttyp
              end if ! if rad_error_model > 0
            end do ! on i n_letkf_obserr_rad
          end if ! if hyperspectral
        end do ! end do i1
      end do ! end do is

      call destruct (cloud_fl_tmp)

      deallocate (s_state)
      deallocate (s_check)
      deallocate (b_state)
      deallocate (b_check)
      deallocate (b_fg   )

      !--------------------------------
      ! deallocate current input buffer
      !--------------------------------
      do i=2,n
        call destruct (obs_in% o(i))
      end do
      deallocate      (obs_in% o    )
      deallocate      (mh)
      deallocate      (mb)
      deallocate      (mi)


      !-----------------------------------
      ! read individual ensemble forecasts
      ! and mean pressure level
      !-----------------------------------
      call            construct (Hx_tmp(ifi,:),       obs_tmp(ifi)% oi, global=.true.)
      call            construct (Hx_tmp_clear(ifi,:), obs_tmp(ifi)% oi, global=.true.)
      call            construct (cloud_impact(ifi),   obs_tmp(ifi)% oi, global=.true.)
      if (l_det) call construct (Hd_tmp(ifi),         obs_tmp(ifi)% oi, global=.true.)
      allocate (pl_tmp   (n_obs))
      allocate (pl_mean  (n_obs)); pl_mean  = 0._wp
      allocate (pl_high(n_obs));   pl_high = 100000._wp
      allocate (pl_valid (n_obs)); pl_valid = .true.
      if (lw_f) then
        allocate (w_f_mean  (vert_lev, n_obs))
        allocate (p_wf_mean (vert_lev, n_rep))
        w_f_mean   = 0._wp
        p_wf_mean  = 0._wp
      endif

      k          = 1; if (l_det) k          = 0
      ens_member = 1; if (l_det) ens_member = VE_DETERM
      cloud_impact(ifi)%s(1)%x = 0  !initialisation
      do i=1,n
        if (source(i)% pe == dace% pe) then
          call open_fdbk_read (fb, path_file(source(i)% path, source(i)% file))
          !------------------------------
          ! search for runtype FIRSTGUESS
          !------------------------------
          call read_meta      (fb)
          runtype                                = VT_FIRSTGUESS
          if (fof_ver .or. fof_ekf) query_member = .true.
          do j = 1, 2
            if (query_member) then ! single file (for ver*.nc/mof*.nc/ekf*.nc)
              call get_veri_index (ifg, nidx, fb, run_type=runtype,     &
                                                  ens_member=ens_member )
              !clear-sky TB for seviri
              call get_veri_index (ifg_clear, nidx_clear, fb,   &
                                   run_type=runtype,            &
                                   ens_member=ens_member,       &
                                   operator_flag=OF_BT_CLEAR_SKY)
            else   ! one file per member for fof* files
              call get_veri_index (ifg, nidx, fb, run_type=runtype)

              !clear-sky TB for seviri
              call get_veri_index (ifg_clear, nidx_clear, fb,   &
                                   run_type=runtype,            &
                                   operator_flag=OF_BT_CLEAR_SKY)
            endif
            !-----------------------------------------
            ! search for runtype FORECAST if not found
            !-----------------------------------------
            if (j==1 .and. nidx /= 1 .and. fof_ver) then
              runtype = VT_FORECAST
            else
              exit
            endif
          end do
          !-------------------------------
          ! abort if runtype was not found
          !-------------------------------
          if (nidx /= 1) then
            call print_fdbk (fb)
            write(6,*)   'read_obs_ens: cannot find First Guess in '//&
                          trim(source(i)%file)
            write(6,*)   'required run_type  :',runtype
            if (check_member) &
              write(6,*) 'required ens_member:',ens_member
            write(6,*)   'get_veri_index returns       :',ifg,nidx
            call finish ('read_obs_ens',                                    &
                         'cannot find First Guess in '//trim(source(i)%file))
          endif
          !--------------------------
          ! actually read first guess
          !--------------------------
          ix(i)% Hx       = get_veri (fb, ifg(1), n=size(ix(i)% Hx))

          if (k==0) then
            Hd_tmp(ifi)  % s(1)% x = ix(i)% Hx (ix(i)% ib)
          else
            Hx_tmp(ifi,k)% s(1)% x = ix(i)% Hx (ix(i)% ib)
          endif

          if (nidx_clear == 1 .and. k /= 0) then  !get clear-sky first guess
              ix(i)% Hx       = get_veri (fb, ifg_clear(1), n=size(ix(i)% Hx))
              Hx_tmp_clear(ifi,k)% s(1)% x = ix(i)% Hx (ix(i)% ib)

              !!compute cloud impact needed for error model
              !! Sum up for all members to obtain average
              cloud_impact(ifi)%s(1)%x = cloud_impact(ifi)%s(1)%x + &
                   0.5_wp * (ABS(obs_tmp(ifi)% o(1)%body%o - Hx_tmp_clear(ifi,k)% s(1)% x) &
                          + (ABS (Hx_tmp(ifi,k)% s(1)% x   - Hx_tmp_clear(ifi,k)%s(1)% x)))
          end if



          if (k/=0) then
            allocate (pl_in (fb% n_body))
            call get_var (fb, 'plevel', pl_in, fill=-1._wp)
            pl_tmp   = pl_in (ix(i)% ib)
            pl_valid = pl_valid .and. pl_tmp > 0._wp
            pl_mean  = pl_mean + pl_tmp
            !find member with topmost plevel:
            pl_high = min(pl_tmp,pl_high)
            deallocate (pl_in)
          endif

          num_obs  = fb% n_body
          num_head = fb% n_hdr

          !----------------------------------------------
          ! fix bug (erroneously rotated wind components)
          !----------------------------------------------
          luvbug = ( fb% source      (1:5) =='COSMO' .and. &
                     fb% institution (1:5) /='CNMCA' .and. &
                     fb% institution (1:5) /='COMET' .and. &
                    (fb% version           ==' 1.01' .or.  &
                     fb% version           ==' 1.00').and. &
                    (fb% pole(1)           /=  0._sp .or.  &
                     fb% pole(2)           /=  0._sp)      )
          if (luvbug) then
            do i1 = 1, n_rep
              do j1 = o% spot(i1)% o% i + 1 ,&
                      o% spot(i1)% o% i + o% spot(i1)% o% n - 1
                select case (o% varno(j1))
                case (VN_U, VN_U10M)
                  select case (o% varno(j1+1))
                  case (VN_V, VN_V10M)
                    if (k==0) then
                      call uvrot2uv (Hd_tmp(ifi)% s(1)% x(j1), Hd_tmp(ifi)% s(1)% x(j1+1), &
                                     o%spot(i1)%col%c% dlat,   o% spot(i1)%col%c% dlon,    &
                                     real (fb% pole(1), wp),   real (fb% pole(2), wp),     &
                                     ugeo, vgeo                                            )
                      Hd_tmp(ifi)% s(1)% x(j1)   = ugeo
                      Hd_tmp(ifi)% s(1)% x(j1+1) = vgeo
                    else
                      call uvrot2uv (Hx_tmp(ifi,k)% s(1)% x(j1), Hx_tmp(ifi,k)% s(1)% x(j1+1), &
                                     o%spot(i1)%col%c% dlat,     o% spot(i1)%col%c% dlon,      &
                                     real (fb% pole(1), wp),     real (fb% pole(2), wp),       &
                                     ugeo, vgeo                                                )
                      Hx_tmp(ifi,k)% s(1)% x(j1)   = ugeo
                      Hx_tmp(ifi,k)% s(1)% x(j1+1) = vgeo
                    endif
                  end select
                end select
              end do
            end do
          endif

          call close_fdbk   (fb)
          call cleanup_fdbk (fb)

          !-------------------------
          ! read weighting functions
          !-------------------------
          if (lw_f) then
            lrad = .false.
            !--------------------------------------------------
            ! only ensemble members, not from deterministic run
            !--------------------------------------------------
            if (k/=0) then
              status = nf90_open (path_file(source(i)% path,                     &
                                            source(i)% file ), nf90_nowrite, ncid)
              if (status /= nf90_noerr)                                    &
                call finish('read_obs_ens',                                &
                            'cannot open file for weighting functions: '// &
                             trim(source(i)% file)                         )
              status = nf90_inq_varid (ncid,'w_f',varid_w_f)
              if (status == nf90_noerr) then
                lrad = .true.
                allocate (w_f_in  (vert_lev, num_obs ))
                allocate (p_wf_in (vert_lev, num_head))
                do lll = 1, num_obs
                  status = nf90_get_var (ncid,varid_w_f,w_f_in(:,lll),start=(/1,lll/))
                  if (status /= nf90_noerr)                       &
                    call finish("read_obs_ens","cannot read 'w_f'")
                    if (w_f_in(1,lll) == 9.96921e+36) w_f_in(:,lll) = 0._wp
                enddo
                status = nf90_inq_varid (ncid,'p_wf',varid_p_wf)
                if (status /= nf90_noerr)                           &
                  call finish("read_obs_ens","cannot inquire 'p_wf'")
                do lll = 1, num_head
                  status = nf90_get_var (ncid,varid_p_wf,p_wf_in(:,lll),start=(/1,lll/))
                  if (status /= nf90_noerr)                        &
                    call finish("read_obs_ens","cannot read 'p_wf'")
                    if (p_wf_in(1,lll) == 9.96921e+36) p_wf_in(:,lll) = 0._wp
                enddo
                do lll = 1, vert_lev
                  w_f_mean (lll,:) = w_f_mean (lll,:) + w_f_in  (lll,ix(i)% ib)
                  p_wf_mean(lll,:) = p_wf_mean(lll,:) + p_wf_in (lll,ix(i)% ih)
                end do
                deallocate (p_wf_in, w_f_in)
              endif
            endif
            status = nf90_close (ncid)
          endif

          deallocate (ix(i)% ih)
          deallocate (ix(i)% ib)
          deallocate (ix(i)% Hx)
        endif

        k          = k + 1
        ens_member = k!!
      end do    ! read forecasts

      !-----------------------------
      ! Distribute model equivalents
      !-----------------------------
      k = 1; if (l_det) k = 0
      do i = 1, n
        if (k==0) then
          call p_bcast (Hd_tmp(ifi)  % s(1)% x, source(i)% pe)
        else
          call p_bcast (Hx_tmp(ifi,k)% s(1)% x, source(i)% pe)
        endif
        k = k + 1
      end do
      cloud_impact(ifi)% s(1)% x = p_sum (cloud_impact(ifi)% s(1)% x)

      !----------------------------
      ! mean of weighting functions
      !----------------------------
      lrad = p_or (lrad)
      if (lw_f .and. lrad) then
        w_f_mean  = p_sum (w_f_mean)  / k_enkf
        p_wf_mean = p_sum (p_wf_mean) / k_enkf
!c      do lll = 1, obs_tmp(ifi)% o(1)% n_obs
!c        allocate (obs_tmp(ifi)% o(1)% body(lll)% w_f  (vert_lev))
!c        obs_tmp(ifi)% o(1)% body(lll)% w_f (:) = w_f_mean (:,lll)
!c      enddo
!c      do lll = 1, obs_tmp(ifi)% o(1)% n_spot
!c        allocate (obs_tmp(ifi)% o(1)% spot (lll)% p_wf (vert_lev))
!c        obs_tmp(ifi)% o(1)% spot(lll)% p_wf(:) = p_wf_mean(:,lll)
!c      enddo
      endif

      !--------------------
      ! mean pressure level
      !--------------------
      pl_mean  = p_sum (pl_mean) / k_enkf

      pl_high = p_min(pl_high)
      !ens mean cloud impact
      cloud_impact(ifi)%s(1)%x = cloud_impact (ifi)%s(1)%x / k_enkf



      pl_valid = p_and (pl_valid)
      if (l_det) then
        where (pl_valid) obs_tmp(ifi)% o(1)% body% plev = pl_mean
      else
        where (.not. pl_valid) pl_mean = -1._wp
        obs_tmp(ifi)% o(1)% body% plev = pl_mean
      endif


       !--------------------------------------------------
       ! assign highest plevel of all members for SEVIRI satellite radiances
       !----------------------------------------------------
       if (plevel_select == 1) then
         do is = 1, obs_tmp(ifi)%o(1)% n_spot
           si => obs_tmp(ifi)%o(1)% spot(is)
           do i = si% o% i + 1, si% o% i + si% o% n
             if (obs_tmp(ifi)%o(1)%varno(i) == VN_RAWBT .and. si% sttyp == SEVIRI) then
               obs_tmp(ifi)%o(1)%body(i)%plev = pl_high(i)
             end if
           end do
         end do
        end if







      deallocate (pl_tmp, pl_mean, pl_valid, pl_high)
      if (lw_f) deallocate (w_f_mean, p_wf_mean)
    end do      ! loop over files per member

!   if (if1 == 0) call finish('read_obs_ens','no observations present')
    if1 = max (if1,1)

    !-------------------------------------------
    ! collect observations in output variable
    ! (currently keep observations on PE 0 only)
    ! currently only first file
    !-------------------------------------------
999 continue
    allocate (obs_src% o (nfof))
    do ifi = 1, nfof
      obs_src% o(ifi) = obs_tmp(ifi)% o(1)
    end do
    allocate       (obs_out% o    (1))
    call construct (obs_out% o       )
    call join_obs  (obs_src% o (if1:), obs_out% o(1))
    n_rep = obs_out% o(1)% n_spot
    n_obs = obs_out% o(1)% n_obs
    if (n_rep == 0) allocate (obs_out% o(1)% spot (n_rep))
    if (n_obs == 0) allocate (obs_out% o(1)% body (n_obs))
    allocate       (obs_out% oi)
    call construct (obs_out% oi, n_obs, 1, nb=(/n_obs/), pe=(/0/))

    !--------------------------------------
    ! collect forecasts in output variables
    !--------------------------------------
    allocate (H_x (k_enkf))
    call            construct (H_x,        obs_out% oi,'Hi_x', global=.true.)
    call            construct (H_0,        obs_out% oi,'Hi_0', global=.true.)
    if (l_det) call construct (H_det,      obs_out% oi,        global=.true.)
    call            construct (cl_impact,  obs_out% oi,        global=.true.)
    call            construct (ens_cld_fl, obs_out% oi,        global=.true.)
    i   = 0
    H_0 = 0._wp
    do ifi = 1, nfof
      l = obs_src% o(ifi)% n_obs
      if (l == 0) cycle
      do k = 1, k_enkf
        H_x(k)%         s(1)% x(i+1:i+l) = Hx_tmp(ifi,k)% s(1)% x
        H_0   %         s(1)% x(i+1:i+l) = &
          H_0 %         s(1)% x(i+1:i+l) + Hx_tmp(ifi,k)% s(1)% x
      end do
      if (l_det) H_det% s(1)% x(i+1:i+l) = Hd_tmp(ifi)  % s(1)% x
      cl_impact% s(1)% x(i+1:i+l) = cloud_impact(ifi) % s(1)% x
      ens_cld_fl% s(1)% x(i+1:i+l) = cloud_fl_mem(ifi)% s(1)% x

      i = i + l
    end do
    obs_out% o(1)% body% bg = H_0% s(1)% x / k_enkf
    call local (H_x)
    call local (H_det)
    call local (cl_impact)
    call local (ens_cld_fl)

    !-----------------------
    ! deallocate temporaries
    !-----------------------
    call destruct (H_0)
    call destruct (Hd_tmp)
    call destruct (Hx_tmp)
    call destruct (Hx_tmp_clear)
    call destruct (cloud_impact)
    call destruct (cloud_fl_mem)
    call destruct (obs_tmp)
    deallocate    (obs_src% o)

    !------------------
    ! derive statistics
    !------------------
    call derive_rept_stat (obs_out% o)
    call gather_rept_stat
    call print_rept_stat ('after fof read')

  contains
  !--------------------------------------------------------------------------

    function same_head (i, i1, i2) result (ok)
    integer ,intent(in)  :: i, i1, i2
    logical              :: ok
    !--------------------------------------
    ! check consistence of observation head
    !--------------------------------------
      ok = .true.
      if (o% spot(i1)% hd% obstype  /= obs_in% o(i)% spot(i2)% hd% obstype)  ok = .false.
      if (o% spot(i1)% hd% codetype /= obs_in% o(i)% spot(i2)% hd% codetype) ok = .false.
      if (o% spot(i1)% hd% time     /= obs_in% o(i)% spot(i2)% hd% time)     ok = .false.
      if (o% spot(i1)% col% c% dlon /= obs_in% o(i)% spot(i2)% col% c% dlon) ok = .false.
      if (o% spot(i1)% col% c% dlat /= obs_in% o(i)% spot(i2)% col% c% dlat) ok = .false.
      if (o% spot(i1)% statid       /= obs_in% o(i)% spot(i2)% statid)       ok = .false.
      if (o% spot(i1)% actual_time  /= obs_in% o(i)% spot(i2)% actual_time)  ok = .false.
    end function same_head

  !--------------------------------------------------------------------------

    function same_body (i, i1, j1, j2) result (ok)
    !--------------------------------------
    ! check consistence of observation body
    !--------------------------------------
    integer ,intent(in)  :: i       ! box    index 1
    integer ,intent(in)  :: i1      ! report index 1
    integer ,intent(in)  :: j1, j2  ! body indices
    logical              :: ok      ! result variable

      ok = .true.
      select case (o% spot(i1)% hd% obstype)
      case default
        !-------------------------------
        ! default: check VARNO and level
        !-------------------------------
        if (o% varno (j1) /= obs_in% o(i)% varno (j2)) ok = .false.
        if (o% olev  (j1) /= obs_in% o(i)% olev  (j2)) ok = .false.
      case (OT_RADAR)
        !-------------------------------------------------
        ! RADAR: check VARNO and range, azimuth, elevation
        !-------------------------------------------------
        if (o% varno (j1)      /= obs_in% o(i)% varno (j2))      ok = .false.
!        if (o% body  (j1)% l2c /= obs_in% o(i)% body  (j2)% l2c) ok = .false.
        if (o% body  (j1)% spec_index /= obs_in% o(i)% body  (j2)% spec_index) ok = .false.
        if (o% body             (j1)% obs_par(1) /= &
            obs_in% o(i)% body  (j2)% obs_par(1)    ) ok = .false.
        if (o% body             (j1)% obs_par(2) /= &
            obs_in% o(i)% body  (j2)% obs_par(2)    ) ok = .false.
      end select
    end function same_body

  !--------------------------------------------------------------------------
    subroutine pack_header (x, ih, ib, m)
    type (t_spot) ,pointer    :: x  (:)
    integer       ,intent(in) :: ih (:)  ! header indices to keep
    integer       ,intent(in) :: ib (:)  ! body   indices to keep
    integer       ,intent(in) :: m       ! member index

      type (t_spot) ,pointer  :: y  (:)
      integer                 :: i, k, n, on, oi, nn
      integer                 :: ii, mem

      if (.not. associated (x)) return
      allocate (y (size(ih)))
      y = x (ih)
      i = 0
      do k = 1, size(y)                  ! loop over reports
        ii = - huge (0)                  ! previous body index
        oi = y(k)% o% i + 1              ! old first index in spot
        on = y(k)% o% i + y(k)% o% n     ! old last  index in spot
        y(k)% o% i = i                   ! new first index in spot
        do n = 1, y(k)% o% n             ! loop over old report size
          nn = i + n                     ! new last index in spot
          if (size(ib) < nn) then        ! index extends total body size
            y(k)% o% n = n-1             ! set new report size
            i      = i + n-1             ! set first index of next report
            exit                         ! exit report size loop
          endif
          if (on == ib(nn)) then         ! last body entry reached
            y(k)% o% n = n               ! set new report size
            i      = i + n               ! set first index of next report
            exit                         ! exit report size loop
          endif
          if (on <  ib(nn)) then         ! index beyond last body entry
            y(k)% o% n = n-1             ! set new report size
            i      = i + n-1             ! set first index of next report
            exit                         ! exit report size loop
          endif
          if (oi >  ib(nn)) then         ! index below first body entry
            mem = m; if (l_det) mem = m-1
            write (0,*) dace% pe,'read_obs_ens: pack_header WARNING member', mem
            write (0,*) dace% pe,'read_obs_ens: obstype,code,dbkz,t,statid =', &
                 y(k)% hd% obstype, y(k)% hd% codetype, y(k)% hd% dbkz,        &
                 cyyyymmddhhmm (y(k)% actual_time), "  ", trim (y(k)% statid)
            write (0,*) dace% pe,'read_obs_ens: truncating/shrinking report;', &
                 ' new,old =', n-1, y(k)% o% n
            write (0,*)
            call nextline
            write (oline(iol),'(i6,a,i6)'           ) dace% pe,   &
                 ' read_obs_ens: pack_header WARNING member =', mem
            call nextline
            write (oline(iol),'(i6,a,3i6,2x,a,2x,a)') dace% pe,         &
                 ' read_obs_ens: obstype,code,dbkz,t,statid =',         &
                 y(k)% hd% obstype, y(k)% hd% codetype, y(k)% hd% dbkz, &
                 cyyyymmddhhmm (y(k)% actual_time), trim (y(k)% statid)
            call nextline
            write (oline(iol),'(i6,a,2i6)'          ) dace% pe,           &
                 ' read_obs_ens: truncating/shrinking report; new,old =', &
                 n-1, y(k)% o% n
            call nextline
            y(k)% o% n = n-1             ! set new report size
            i      = i + n-1             ! set first index of next report
            exit                         ! exit report size loop
          end if
          if (ii >  ib(nn) .or. &        ! index out of order
              y(k)% o% n == n) then      ! report size reached, no end found
            write (*,*) dace% pe,'read_obs_ens: pack_header failed'
            write (*,*) dace% pe,'read_obs_ens: otyp,bft,sbt,ctype, statid =', &
                 y(k)% hd% obstype, y(k)% hd% buf_type, y(k)% hd% buf_subtype, &
                 y(k)% hd% codetype, y(k)% hd% dbkz,                           &
                 cyyyymmddhhmm (y(k)% actual_time), "  ", trim (y(k)% statid)
            write (0,*)
            write (0,*) dace% pe,'read_obs_ens: pack_header failed'
            write (0,*) dace% pe,'read_obs_ens: otyp,bft,sbt,ctype, statid =', &
                 y(k)% hd% obstype, y(k)% hd% buf_type, y(k)% hd% buf_subtype, &
                 y(k)% hd% codetype, y(k)% hd% dbkz,                           &
                 cyyyymmddhhmm (y(k)% actual_time), "  ", trim (y(k)% statid)
            write (0,*) dace% pe,'read_obs_ens: ny,nx =',size(y), size(x)
            write (0,*) dace% pe,'read_obs_ens: k, ih =',k, ih(k)
            write (0,*) dace% pe,'read_obs_ens: n, n  =',n, y(k)% o% n
            write (0,*) dace% pe,'read_obs_ens: ni,nn =',i+1,nn,n
            write (0,*) dace% pe,'read_obs_ens: oi,on =',oi,on,on-oi+1
            write (0,*) dace% pe,'read_obs_ens: ib    =',ib(i+1:nn)
            write (0,*)
            call finish('read_obs_ens','pack_header')
          endif
          ii = ib(nn)                    ! remember body index
        end do
      end do
      deallocate (x)
      x => y
    end subroutine pack_header
  !--------------------------------------------------------------------------
    subroutine pack_body (x, ix)
    type (t_datum) ,pointer    :: x  (:)
    integer        ,intent(in) :: ix (:)

      type (t_datum) ,pointer  :: y  (:)

      if (.not. associated (x)) return
      allocate (y (size(ix)))
      y = x (ix)
      deallocate (x)
      x => y
    end subroutine pack_body
  !--------------------------------------------------------------------------
    subroutine pack_int (x, ix)
    integer ,pointer    :: x  (:)
    integer ,intent(in) :: ix (:)

      integer ,pointer  :: y  (:)

      if (.not. associated (x)) return
      allocate (y (size(ix)))
      y = x (ix)
      deallocate (x)
      x => y
    end subroutine pack_int
  !--------------------------------------------------------------------------
    subroutine pack_real (x, ix)
    real(wp) ,pointer    :: x  (:)
    integer  ,intent(in) :: ix (:)

      real(wp) ,pointer  :: y  (:)

      if (.not. associated (x)) return
      allocate (y (size(ix)))
      y = x (ix)
      deallocate (x)
      x => y
    end subroutine pack_real
  !--------------------------------------------------------------------------
    subroutine pack_par (o)
    type(t_obs)  ,intent(inout) :: o

      type(t_spot) ,pointer :: s
      integer ,pointer :: y (:)
      integer          :: n
      integer          :: i, j

      n = sum (o% spot% p% n)
      o% n_par =  n
      if (n==0) then
        if (associated(o% par)) deallocate (o% par)
        return
      else
        allocate (y (n))
        j = 0
        do i = 1, o% n_spot
          s => o% spot(i)
          if (s%p% n > 0) then
            y (j + 1 : j + s%p% n) = o% par (s%p% i + 1 : s%p% i + s%p% n)
            s%p% i = j
            j      = j + s%p% n
          endif
        end do
        deallocate (o% par)
        o% par   => y
      endif
    end subroutine pack_par
  !--------------------------------------------------------------------------
  end subroutine read_obs_ens
!==============================================================================
  subroutine read_det (det, grid, time, file, par_read, runtype, optionals,   &
        relabel, reftime, leveltypes, gridfile, invar, geosp, lsm, hhl, ierr, &
        g_coarse, inv                                                         )
  !========================================
  ! read deterministic forecast or analysis
  !========================================
  type (t_atm)      ,intent(out) :: det       ! atmospheric state to read
  type (t_grid)     ,pointer     :: grid      ! atmospheric grid  to read
  character(len=*)  ,intent(in)  :: file      ! file name
  type (t_time)     ,intent(in)  :: time      ! valid time of atmospheric state
  character(len=*)  ,intent(in)  :: par_read  ! parameters to read
  character(len=*)                           &!
          ,optional ,intent(in)  :: runtype   ! forecast, analysis
  character(len=*)                           &!
          ,optional ,intent(in)  :: optionals ! parameters may be missing
  character(len=*)                           &!
          ,optional ,intent(in)  :: relabel   ! parameters to re-label
  type(t_time)                               &!
          ,optional ,intent(in)  :: reftime   ! start of forecast
  integer ,optional ,intent(in)  :: leveltypes(:) ! leveltypes to read
  character(len=*)                           &!
          ,optional ,intent(in)  :: gridfile  ! grid metadata (ICON)
  character(len=*)                           &!
          ,optional ,intent(in)  :: invar     ! invariant fields file
  logical ,optional ,intent(in)  :: geosp     ! passed to read_grid
  logical ,optional ,intent(in)  :: lsm       ! passed to read_grid
  logical ,optional ,intent(in)  :: hhl       ! passed to read_grid
  integer ,optional ,intent(out) :: ierr      ! error return code
  type (t_inventory)                         &!
          ,optional ,pointer     :: inv(:)    ! GRIB file inventory
  !-----------------------------
  ! LETKF coarse grid generation
  !-----------------------------
  type(t_grid)    ,pointer    ,optional :: g_coarse     ! coarse grid

    integer                    :: iunit   ! unit number for file presence check
    integer                    :: ios     ! iostat return variable
    type(t_inventory) ,pointer :: invt(:) ! GRIB file inventory

    !--------------------------------------------
    ! derive file name, check if file is readable
    !--------------------------------------------
    if (present (ierr)) ierr = 0
    if (dace% lpio) then
      iunit = get_unit_number()
      open (iunit, file=file, status='old', action='read', iostat=ios)
      if (ios/=0) then
        if (present (ierr)) then
          ierr = IE_FILE
        else
          call finish ('read_det','cannot open '//trim(file))
        endif
      endif
      close (iunit)
      call return_unit_number (iunit)
    end if
    if (present (ierr)) then
      ierr = p_min (ierr)
      if (ierr /= 0) return
    endif

    !---------------
    ! read inventory
    !---------------
    nullify (invt)
    if (present(inv)) then
      if (.not.associated(inv)) call finish ('read_det','inventory not associated')
      if ( size(inv) == 0 )     call finish ('read_det','inventory is empty')
      invt => inv
    else
      call get_inventory (invt, file)
    end if
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '    deterministic run file inventory'
      write(6,'()')
      call print_inventory (invt, first=.true.)
    endif
    !------------
    ! derive grid
    !------------
    if (.not. associated (grid)) then
      if (dace% lpio) then
        write(6,'()')
        write(6,'(a)') '    reading grid'
        write(6,'()')
      endif
      allocate  (grid)
      call read (grid, file, gridfile=gridfile, invt=invt, invar=invar, &
                 geosp=geosp, lsm=lsm, hhl=hhl, g_coarse=g_coarse,      &
                 nproc1=nproc1, nproc2=nproc2, comm=dace% comm          )
    else if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '    grid not read, already present'
      write(6,'()')
    endif
    !-----------------------
    ! read atmospheric state
    !-----------------------
    if (dace% lpio) then
      write(6,'(a,i4)') '    reading file: '//trim(file)
    endif
    call construct (det, grid, alloc=trim(par_read))
    if (present (relabel)) call relabel_fields (invt, par=relabel)
    call read (det, file, invt,         &
               time       = time,       &
               reftime    = reftime,    &
               fields     = par_read,   &
               leveltypes = leveltypes, &
               runtype    = runtype,    &
               optionals  = optionals,  &
               ierr       = ierr        )
    !---------------
    ! error handling
    !---------------
    if (present(inv)) then
      nullify (invt)
    else
      if (associated(invt)) deallocate (invt)
    end if
    if (present (ierr)) then
      if (ierr /= 0) return
    endif

    !---------
    ! printout
    !---------
    if (dace% lpio) write(6,'()')
    call print (det, comment='deterministic forecast ')

  end subroutine read_det
!==============================================================================
  subroutine read_ens (ens, grid, n_ens, file, time, c2a, par_read, runtype,  &
                       optionals, unsp_type, members, reftime, leveltypes,    &
                       gridfile, suffix, invar, geosp, lsm, hhl, verbose,     &
                       read_seq, diamond, ierr,                               &
                       rf_coarse, ni_coarse, nzr_coarse, g_coarse, ctr        )
  !====================
  ! read ensemble state
  !====================
  type(t_atm)     ,pointer              :: ens (:)      ! ensemble state
  type(t_grid)    ,pointer              :: grid         ! ensemble grid
  integer         ,intent(in)           :: n_ens        ! ensemble size
  character(len=*),intent(in)           :: file         ! file name
  type (t_time)   ,intent(in)           :: time         ! valid time of atm. state
  logical         ,intent(in)           :: c2a          ! transform to A grid
  character(len=*),intent(in)           :: par_read     ! parameters to read
  character(len=*),intent(in) ,optional :: runtype      ! forecast, analysis
  character(len=*),intent(in) ,optional :: optionals    ! parameters may missed
  character(len=*),intent(in) ,optional :: unsp_type    ! unspecified runtype
  integer         ,intent(in) ,optional :: members(:)   ! members to read
  type (t_time)   ,intent(in) ,optional :: reftime      ! start of forecast
  integer         ,intent(in) ,optional :: leveltypes(:)! leveltypes to read
  character(len=*),intent(in) ,optional :: gridfile     ! grid metadata (ICON)
  character(len=*),intent(in) ,optional :: suffix       ! suffix
  character(len=*),intent(in) ,optional :: invar        ! invariant fields file
  logical         ,intent(in) ,optional :: geosp        ! Add orography
  logical         ,intent(in) ,optional :: lsm          ! Add land-sea-mask
  logical         ,intent(in) ,optional :: hhl          ! Add HHL
  integer         ,intent(in) ,optional :: verbose      !1.timings 2:inventory 3:fields
  logical         ,intent(in) ,optional :: read_seq     ! read ensemble sequentially
  integer         ,intent(in) ,optional :: diamond      ! keep only this diamond
  integer         ,intent(out),optional :: ierr         ! 0: OK, >0: missing member id
  type(t_ctr_inv) ,intent(in) ,optional :: ctr(:)       ! inventory container
  !-----------------------------
  ! LETKF coarse grid generation
  !-----------------------------
  integer         ,intent(in) ,optional :: rf_coarse    ! coarse grid factor
  integer         ,intent(in) ,optional :: ni_coarse    ! coarse grid GME resolution
  integer         ,intent(in) ,optional :: nzr_coarse   ! coarse grid no.levels
  type(t_grid)    ,pointer    ,optional :: g_coarse     ! coarse grid

    integer                    :: k       ! ensemble member index
    integer                    :: ke      ! ensemble member number
    character(len=128)         :: fname   ! full file name
    character(len=128)         :: iname   ! invariant fields file name
    integer                    :: iunit   ! unit number for file presence check
    integer                    :: ios     ! iostat return variable
    type(t_inventory) ,pointer :: invt(:) ! GRIB file inventory
    integer                    :: verb    ! verbosity
    logical                    :: lseq    ! read sequentially
    integer                    :: ldia    ! keep only this diamond

    lseq = .false.; if (present (read_seq)) lseq = read_seq
    ldia = 0;       if (present (diamond))  ldia = diamond
    verb = 2;       if (present (verbose))  verb = verbose
                    if (present (ierr))     ierr = 0
    if (verb>0) call stop_time ('read forecast ensemble')

    !------------------------------------
    ! loop over forecast ensemble members
    !------------------------------------
    if (n_ens == 0) call finish('read_ens','n_ens == 0')
    do k = 1, n_ens
      ke = k; if(present(members)) ke = members(k)
      !--------------------------------------------
      ! derive file name, check if file is readable
      !--------------------------------------------
      fname = path_file (input, file, suffix=suffix, iens=ke)
      if (present(ctr)) then
        if (fname /= ctr(k)% name) then
          if (present (ierr)) then
            ierr = k
          else
            call finish ('read_ens','name mismatch: '//trim(fname) &
                                      //' /= '//trim(ctr(k)% name) )
          endif
        end if
      else
        if (dace% lpio) then
          iunit = get_unit_number()
          open (iunit, file=fname, status='old', action='read', iostat=ios)
          if (ios/=0) then
            if (present (ierr)) then
              ierr = k
            else
              call finish ('read_ens','cannot open '//trim(fname))
            endif
          endif
          close (iunit)
          call return_unit_number (iunit)
        end if
      end if
      if (present (ierr)) then
        ierr = p_max (ierr)
        if (ierr /= 0) return
      endif
      if (k==1) then
        !---------------------------
        ! allocate pointer variables
        !---------------------------
        nullify  (invt)
        if (.not. associated (ens )) allocate (ens (n_ens))
        if (.not. associated (grid)) then
          allocate (grid)
          !---------------------------------------
          ! derive grid from first ensemble member
          !---------------------------------------
          if (present(ctr)) then
            invt => ctr(1)% invt
          else
            call get_inventory (invt, fname)
          end if
          if (dace% lpio) then
            if (verb>1) then
              write(6,'()')
              call print_inventory (invt, first=.true.)
              call print_inventory (invt, first=.true., liname=.true.)
            endif
            write(6,'()')
            write(6,'(a)') '    reading grid'
            write(6,'()')
          endif
          iname = ''; if (present(invar)) iname = path_file (input, invar)
          call read (grid, fname, gridfile=gridfile, invt=invt,       &
                     invar=iname, lsm=lsm, geosp=geosp, hhl=hhl,      &
                     optionals=optionals, rni=ni_coarse,              &
                     rf=rf_coarse, nzr=nzr_coarse, g_coarse=g_coarse, &
                     nproc1=nproc1, nproc2=nproc2, comm=dace% comm    )
        else if (dace% lpio) then
          write(6,'()')
          write(6,'(a)') '    grid not read, already present'
          write(6,'()')
        endif
      endif
    end do

    !---------------------------------------------
    ! read atmospheric state ensemble sequentially
    !---------------------------------------------
    if (lseq) then
      do k = 1, n_ens
        if (verb>0) call stop_time('read forecast ensemble member '//char3(k))
        fname = path_file (input, file, suffix=suffix, iens=k)
        if (dace% lpio) then
          write(6,'(a,i4)') '    reading file: '//trim(fname)
        endif
        if (present(ctr)) then
          invt => ctr(k)% invt
        else
          call get_inventory   (invt, fname)
        end if
        call construct (ens(k), grid, alloc=trim(par_read))
        call read (ens(k), fname, invt,     &
                   time    = time,          &
                   fields  = par_read,      &
                   runtype = runtype,       &
                 optionals = optionals,     &
                 unsp_type = unsp_type,     &
                   reftime = reftime,       &
               leveltypes  = leveltypes,    &
                      ierr = ierr           )
        call write_profiles (ens(k),'enkf_fc_'//char3(k),       &
                             'LETKF forecast member '//char3(k) )
        call keep_diamond (ldia, ens(k))
      end do
      call keep_diamond (ldia, grid=grid)
    else
    !--------------------------------------------
    ! read atmospheric state ensemble in parallel
    !--------------------------------------------
      if (verb>0) call stop_time('read forecast ensemble in parallel')
      fname = path_file (input, file, suffix=suffix)
      call construct (ens, grid, alloc=trim(par_read))
      call read (ens, fname,                  &
                 members     = members,       &
                 time        = time,          &
                 fields      = par_read,      &
                 runtype     = runtype,       &
                 optionals   = optionals,     &
                 unsp_type   = unsp_type,     &
                 reftime     = reftime,       &
                 leveltypes  = leveltypes,    &
                 ctr         = ctr,           &
                 ierr        = ierr           )
    endif
    if (present (ierr)) then
      if (ierr /= 0) return
    endif
    !-------------------
    ! print the ensemble
    !-------------------
    if (verb>2) then
      do k = 1, n_ens
        if (dace% lpio) write (6,'()')
        call print (ens(k),                                       &
                    comment='forecast ensemble member '//char3(k),&
                    verbose=verb>3                                )
      end do
    endif
    !------------------------------------------------------------
    ! rotate (u,v)-winds and transform to Arakawa-A grid
    ! transform wind from rotated coordinate system
    ! not required for 'pure' LETKF (no obs.operator application)
    !------------------------------------------------------------
    if (method /= 'LETKF') then
      if (verb>0) call stop_time('read_ens: c to a grid')
      do k = 1, n_ens
        select case (grid% gridtype)
        case (WMO6_LATLON, WMO6_ROTLL)
          if (c2a) then
            if (grid% arakawa == 'C') &
                 call C_to_A_grid (ens(k), save=.true.)
            if (grid% rot)            &
                 call uvtrafo(ens(k),1)            ! rotate wind
          end if
!         if (grid% vct /= VCT_P_HYB) then
!           call set_p  (ens(k))
!           if (.not. associated(ens(k)% ps)) call set_ps (ens(k))
!         endif
        case default
        end select
      end do
    end if
    if (present(ctr)) then
      nullify (invt)
    else
      if (associated(invt)) deallocate (invt)
    end if
  end subroutine read_ens
!=============================================================================

  subroutine derive_obs_ana (obs, H_x0_ana, zlp, grid,                       &
                             lbi, ubi, lbj, ubj, nend, k_enkf, n_ana, h_red, &
                             v_red, fc_mean,                                 &
                             pr, w_q, rmask, write_gain, Wr, rh,             &
                             H_x_ana, ldet, H_det_ana, wq_det, adap_loc      )
  !----------------------------------------------------------
  ! derive (linear) estimate of analysis in observation space
  !----------------------------------------------------------
  type(t_obs_set),intent(inout):: obs     ! observation derived type
  type(t_vector) ,intent(inout):: H_x0_ana! mean of H appl. to ana-ens
  type(t_vector) ,intent(in)   :: zlp     ! log of press levels of obs
  type(t_grid)   ,intent(in)   :: grid    ! coarse grid for weights
  integer        ,intent(in)   :: lbi     ! lower bound x direction
  integer        ,intent(in)   :: ubi     ! upper bound x direction
  integer        ,intent(in)   :: lbj     ! lower bound y direction
  integer        ,intent(in)   :: ubj     ! upper bound y direction
  integer        ,intent(in)   :: nend    ! number of vertical levels
  integer        ,intent(in)   :: k_enkf  ! number of ens members
  integer        ,intent(in)   :: n_ana   ! analysis 'step'
  logical        ,intent(in)   :: h_red   ! horizontal coarse grid used
  logical        ,intent(in)   :: v_red   ! vertical coarse grid used
  type(t_atm)    ,intent(in)   :: fc_mean ! forecast mean (for pressure)
  real(wp)       ,intent(in)   :: pr (:)  ! press. values on reduced grid (Pa)
  real(wp)       ,intent(in)   :: w_q(lbi:ubi,lbj:ubj,nend,k_enkf)!
                                          ! weight vector to be interpolated
  type(t_vector) ,intent(in)   :: rmask   ! = 1 if obs valid
  logical        ,intent(in)   :: write_gain         ! flag to write gain (WH)
  real(wp)       ,intent(in)   :: Wr(lbi:,lbj:,:,:,:)! weight matrix
  real(wp)       ,intent(in)   :: rh(lbi:,lbj:,:    )! hor.loc.length scale
  type(t_vector) ,intent(inout):: H_x_ana (:)        ! WH
  logical        ,intent(in)   :: ldet               ! update det.run
  type(t_vector) ,intent(inout):: H_det_ana          ! deterministic run
  real(wp)       ,intent(in)   :: wq_det(lbi:ubi,lbj:ubj,nend,k_enkf)
                                          ! weight vector for det. run
  logical        ,intent(in)   :: adap_loc! adaptive localisation used
    !----------------
    ! local variables
    !----------------
    integer                  :: nb          ! number of 'boxes'
    type(t_obs) ,pointer     :: o           ! ptr. to observation derived type
    type(t_spot),pointer     :: spt         ! ptr. to spot
    real(wp)    ,allocatable :: w  (:,:)    ! weights
    integer     ,allocatable :: ix (:,:,:)  ! indices
    integer     ,allocatable :: ixl(:,:,:)  ! indices (locally on PE)
    integer                  :: i,j,k,l,m   ! report index, loop indices
    integer                  :: np          ! number of neighbours
    integer                  :: ns          ! number of spots on this PE
    integer                  :: nobs        ! number of obs in spots
    integer                  :: nlev        ! number of levs in spots
    integer                  :: nlev_max    ! max. number of levs in all spots
    real(wp)    ,allocatable :: w_q_n(:,:)  ! weight vector at obs loc (at
                                            ! model levels)
    real(wp)    ,allocatable :: wq_nd(:,:)  ! weight vector at obs loc (at
                                            ! model levels)
    real(wp)    ,allocatable :: lh_n   (:)  ! hor.loc.length at obs loc
    real(wp)    ,allocatable :: w_q_o(:,:)  ! weight vector at obs loc
    real(wp)    ,allocatable :: wq_od(:,:)  ! weight vector at obs loc
    real(wp)    ,allocatable :: lh_o   (:)  ! hor.loc.length at obs loc
    real(wp)    ,allocatable :: x  (:,:)    ! pressure values on reduced grid
    real(wp)    ,allocatable :: xf (:,:)    ! pressure values in obs space
    integer     ,allocatable :: ixo(:,:)    ! index array, relates obs to levels
    integer     ,allocatable :: kf (:,:)    ! index array, relates obs to levels
    integer     ,allocatable :: jf   (:)    ! index array, relates obs to levels
    real(wp)                 :: z           ! level in log(p)
    integer                  :: i1,i2,i3,i4
    integer                  :: k1,k2
    real(wp)                 :: H_x_ana_mean
    real(wp)                 :: Hxt(k_enkf) ! temporary
    type(t_vector)           :: tmp         ! temporary
    type(t_vector)           :: tmph        ! temporary
    type(t_vector)           :: tmpd        ! temporary
    real(wp)                 :: wv          ! vertical interpolation weight
    logical                  :: ok

    !==========================
    ! derive horizontal weights
    !==========================
    call stop_time ('derive horizontal weights')

    nb = size (obs% o)
    if (nb /= 1) call finish ('derive_obs_ana','nb /= 1')
    o => obs% o(1)
    allocate (w  (   mp, o% n_spot))
    allocate (ix (mp, 4, o% n_spot))
!   w  = 0._wp
!   ix = 0
!   do i = 1+dace% pe, o% n_spot, dace% npe
    do i = 1,      o% n_spot
      call Grid_Indices                   &
           (o% spot(i)% col%c% dlon,      & ! longitude
            o% spot(i)% col%c% dlat,      & ! latitude
                               grid,      & ! grid data type
                               ix(:,:,i), & ! Grid point indices
                               w (:,  i), & ! Weight
                               np         ) ! # of coefficients
      if (.not. h_red) then
        !-----------------------------------------------
        ! for native grid fall back to nearest neighbour
        !-----------------------------------------------
        j = sum (maxloc (w(1:np,i)))
        w (1   ,i) = w (j  ,i)
        w (2:  ,i) = 0._wp
        ix(1 ,:,i) = ix(j,:,i)
        ix(2:,:,i) = 0
      endif
    end do
!   w  = p_sum (w)
!   ix = p_sum (ix)
    np = mp; if (.not. h_red) np = 1

    !---------------------------------
    ! relate observations to processor
    ! +++ will fail on GME grid +++
    !---------------------------------
    call stop_time ('relate observations to processors')

    ok        = .true.
    ix(:,4,:) = -1
    do i = 1, o% n_spot
      if (all (ix(:np,1,i) >= lbi .and. &
               ix(:np,1,i) <= ubi .and. &
               ix(:np,2,i) >= lbj .and. &
               ix(:np,2,i) <= ubj .and. &
               ix(:np,3,i) == ix(1,3,i))) ix(1,4,i) = dace% pe
    end do
    !-------------------------------------------------
    ! check if all neighbour gridpoints are on same PE
    ! for valid observations
    !-------------------------------------------------
    ix(1,4,:) = p_max (ix(1,4,:))
    do i = 1, o% n_spot
      if (ix(1,4,i) < 0 .and. o% spot(i)% use% state >= STAT_PAS_REJ) then
        if (ok) write(0,*) "derive_obs_ana: cannot relate obs to pe:"
        write(0,*) "derive_obs_ana: obstype codetype ident statid=", &
             o% spot(i)% hd% obstype, o% spot(i)% hd% codetype,      &
             o% spot(i)% ident, trim (o% spot(i)% statid)
        ok = .false.
      endif
    end do
    if (.not. ok) call finish ('derive_obs_ana','cannot relate obs to pe')

    !---------------------------------------------------------------
    ! pack index array ix:
    !   1st index: (1:3 or 1:4): neighbour grid points
    !   2nd index: 1           : 1st horizontal grid index
    !              2           : 2nd horizontal grid index
    !              3           : diamond        grid index
    !              4           : PE index (same for all neighbours)
    !   3rd index:             : reports (field of views)
    !     entries 1..ns hold reports handled on this PE
    !   misuse entry (2,4,:)   : original report index
    ! pack weights w:
    !   1st index: (1:3 or 1:4): neighbour grid points
    !   2nd index:             : reports (field of views)
    !--------------------------------------------------------------
    ns = 0
    do i = 1, o% n_spot
      if (ix(1,4,i) == dace% pe) then
        ns = ns + 1
        ix(:,1:3,ns) = ix(:,1:3,i)
        ix(:, 4 ,ns) = i
        w (:    ,ns) = w (:    ,i)
      endif
    end do
    allocate (ixl (mp, 4, ns))
    ixl = ix (:,:,1:ns)
    deallocate (ix)

    !--------------------------------------
    ! allocate for horizontal interpolation
    !--------------------------------------
    allocate   (jf    (ns))
    allocate   (w_q_n (k_enkf, nend))
    allocate   (lh_n          (nend))
    if (ldet) allocate (wq_nd (k_enkf, nend))

    !-------------------------------------------------
    ! allocate for vertical interpolation with splines
    !-------------------------------------------------
    if (nend > 1) then
      allocate (x   (k_enkf, nend))
      do k=1,k_enkf
        x(k,:) = log(pr(:)) !coarse grid pressure values
      end do
    endif

    !--------------------------------------------------------
    ! determine nlev_max (only weights for each level needed)
    ! nlev_max corresponds to nobs_max here
    !--------------------------------------------------------
    nlev_max = 0
    do i=1,ns
      spt => o % spot(ixl(2,4,i))
      nlev = spt% o% n
      if (nlev > nlev_max) nlev_max = nlev
    end do
    !-----------------------
    ! nlev_max = nobs_max...
    !-----------------------
    allocate (xf    (k_enkf, nlev_max))
    allocate (w_q_o (k_enkf, nlev_max))
    allocate (lh_o          (nlev_max))
    allocate (ixo   (ns,     nlev_max))
    allocate (kf    (ns,     nlev_max))
    xf    = 0._wp
    ixo   = 0
    kf    = 0
    if (ldet) allocate (wq_od (k_enkf, nlev_max))

    call   construct (tmp,  obs% oi, 'tmp',  global=.true.)
    call   construct (tmph, obs% oi, 'tmph', global=.true.)
    tmp  = 0._wp
    tmph = 0._wp
    if (ldet) then
      call construct (tmpd, obs% oi, 'tmpd', global=.true.)
      tmpd = 0._wp
    endif

    call stop_time ('estimate observation-analysis')
    !-----------------------------------------
    ! loop over reports (horizontal locations)
    !-----------------------------------------
    do i = 1, ns
      w_q_o = 0._wp
      lh_o  = 0._wp
      if (ldet) then
        wq_od = 0._wp
      endif
      !=========================
      ! horizontal interpolation (interpolate columns to obs locations)
      !=========================
      jf (i:i) = maxloc (w(1:np,i))
!$omp parallel do private(j,k)
      do k = 1, nend
        w_q_n (:,k) = 0._wp
        lh_n    (k) = 0._wp
        if (ldet) then
          wq_nd(:,k) = 0._wp
        endif
        do j = 1,np
          w_q_n(:,k) = w_q_n(:,k) +                           &
                        w(j,i) * w_q(ixl(j,1,i),ixl(j,2,i),k,:)
          lh_n   (k) = lh_n   (k) +                         &
                        w(j,i) * rh (ixl(j,1,i),ixl(j,2,i),k)
          if (ldet)                                               &
            wq_nd(:,k) = wq_nd(:,k) +                             &
                         w(j,i) * wq_det(ixl(j,1,i),ixl(j,2,i),k,:)
        end do
      end do
!$omp end parallel do

      !=======================
      ! vertical interpolation
      !=======================
      if (nend > 1) then
        if (.not. v_red) then
          !------------------------------------------
          ! for native grid use exact pressure values
          !------------------------------------------
          do k=1,k_enkf
            i1 = ixl  (1,1,i)
            i2 = ixl  (1,2,i)
            i4 = ixl  (1,3,i)
            x(k,:) = log (fc_mean% pf (i1,i2,:,i4))
          end do
        endif
        spt => o % spot(ixl(2,4,i))
        z = huge(z)
        nobs = spt% o% n
        nlev = 0
        do j = 1, nobs!loop over obs in spot
          k = spt% o% i + j
          if (zlp% s(1)% x(k) /= z) then
            z = zlp% s(1)% x(k)      ! z = log (plev)
            nlev=nlev+1
            xf( : ,nlev) = z
            kf(i:i,nlev) = minloc (abs(z-x(1,:)))
          end if
          !---------------------------------------------
          ! relate obs to level and store in index array
          !---------------------------------------------
          ixo(i,j) = nlev
        end do
        do k = 1, nlev
          j = kf(i,k)
          if      (xf(1,k) >= x(1,nend))    then
            j  = j - 1
            wv = 1._wp
          else if (xf(1,k) <= x(1,1))       then
            wv = 0._wp
          else
            if (xf(1,k) < x(1,j)) j = j - 1
            wv = (xf(1,k) - x(1,j)) / (x(1,j+1) - x(1,j))
          endif
          w_q_o  (:,k) = wv * w_q_n(:,j+1) + (1._wp-wv) * w_q_n(:,j)
          lh_o     (k) = wv * lh_n   (j+1) + (1._wp-wv) * lh_n   (j)
          if (ldet) &
            wq_od(:,k) = wv * wq_nd(:,j+1) + (1._wp-wv) * wq_nd(:,j)
        end do
      !==========================
      ! no vertical interpolation
      !==========================
      else
        w_q_o (:,1) = w_q_n(:,1)
        lh_o    (1) = lh_n   (1)
        wq_od (:,1) = wq_nd(:,1)
        ixo   (i,:) = 1
        kf    (i,:) = 1
      endif

      !================
      ! derive analysis
      !================

      !-------------------------------------------------------------
      ! now relate w_q_o (wq at obs spots on levels) to observations
      ! (stored in H_x)
      !-------------------------------------------------------------
      spt => o % spot(ixl(2,4,i))
      nobs = spt% o% n
      do j = 1, nobs
        m     = spt% o% i + j     ! observation index
        l     = ixo(i,j)          ! level       index
!       if (rmask% s(1)% x(m) /= 1._wp) cycle  ! check if rmask == 1
        do k  = 1, k_enkf         ! member      index
          tmp% s(1)% x(m)  = tmp%        s(1)% x(m)            &
                           + H_x_ana(k)% s(1)% x(m) * w_q_o(k,l)
          if (ldet)                                              &
            tmpd% s(1)% x(m) = tmpd%       s(1)% x(m)            &
                             + H_x_ana(k)% s(1)% x(m) * wq_od(k,l)
        end do
        tmph% s(1)% x(m) = tmph% s(1)% x(m) + lh_o (l)
      end do
    end do

    tmp%            s(1)% x = p_sum (tmp  % s(1)% x)
    tmph%           s(1)% x = p_sum (tmph % s(1)% x)
    if (ldet) tmpd% s(1)% x = p_sum (tmpd % s(1)% x)

    where (H_x0_ana% s(1)% x /= rvind) &
      H_x0_ana% s(1)% x = H_x0_ana% s(1)% x + tmp% s(1)% x
    if (adap_loc) then
      where (obs% o(1)% body% set% ekf_pass == n_ana      &
       .and. obs% o(1)% body% set% h_loc    <= 0._sp      &
       .and. obs% o(1)% body% use% state    >= STAT_ACTIVE) &
      obs% o(1)% body% set% h_loc = - tmph% s(1)% x
    endif
    if (ldet) then
      where (H_x0_ana%  s(1)% x /= rvind .and. &
             H_det_ana% s(1)% x /= rvind       )
        H_det_ana% s(1)% x = H_det_ana% s(1)% x + tmpd% s(1)% x
      elsewhere
        H_det_ana% s(1)% x = rvind
      endwhere
    endif

    !==============================
    ! derive Kalman gain (i.e. W H)
    !==============================
    if (write_gain) then
      tmp = 0._wp
      call stop_time ('derive Kalman gain')
      do i = 1, ns
        spt => o % spot(ixl(2,4,i))
        i1 = ixl  (jf(i),1,i)
        i2 = ixl  (jf(i),2,i)
        do j = 1, spt% o% n
          m  = spt% o% i + j
          l  = ixo (i,j)
          i3 = kf  (i,l)
          if (H_x0_ana% s(1)% x(m) /= rvind) then
            tmp% s(1)% x(m) = 1._wp
            Hxt             = 0._wp
            H_x_ana_mean    = 0._wp
!           if (rmask% s(1)% x(m) /= 1._wp) cycle  ! check if rmask == 1
            do k1 = 1, k_enkf
              do k2 = 1, k_enkf
                Hxt (k1) = Hxt    (k1)             &
                         + H_x_ana(k2)% s(1)% x(m) &
                         * Wr (i1,i2,i3,k2,k1)
              end do
              H_x_ana_mean = H_x_ana_mean + Hxt (k1)
            end do
            H_x_ana_mean = H_x_ana_mean / k_enkf
            Hxt = Hxt    - H_x_ana_mean
            do k1 = 1, k_enkf
              H_x_ana(k1)% s(1)% x(m) = Hxt (k1)
            end do
          endif
        end do
      end do
      do k1 = 1, k_enkf
        H_x_ana(k1)% s(1)% x = tmp% s(1)% x * H_x_ana(k1)% s(1)% x
        H_x_ana(k1)% s(1)% x = p_sum         (H_x_ana(k1)% s(1)% x)
      end do
    endif
    call destruct           (tmp)
    call destruct           (tmph)
    if (ldet) call destruct (tmpd)

  end subroutine derive_obs_ana

!=============================================================================

 subroutine rms_statistics     (enkf_mean, enkf_spread, enkf_ens, analysis,   &
                                enkf_grid, grid_weight, nend, fct, diag_mode, &
                                xdiag, n_diag, pr, v_red, par_stat, rf, nzr,  &
                                norm_flag, k_enkf, model, input, output,      &
                                ana_time)
  !--------------------------------------------------------------------------
  ! Calculate rms difference between LETKF and 3dVar analysis.
  ! Statistics are calculated on a coarse regular grid
  ! of size nend x n_lon x n_lat.
  ! Output as GRADS readable (ieee) files + GRADS descriptor (.ctl) files.
  !--------------------------------------------------------------------------

  type (t_atm), intent(inout)     ::  enkf_mean       ! LETKF mean   (fg/ana)
  type (t_atm), intent(in)        ::  enkf_spread     ! LETKF spread (fg/ana)
  type (t_atm), intent(in)        ::  enkf_ens (:)    ! LETKF ensemble (fg/ana)
  type (t_atm), intent(in)        ::  analysis        ! 'truth'
  type (t_grid),intent(in)        ::  enkf_grid       ! reduced grid
  type (t_grid),intent(in)        ::  grid_weight     ! points to enkf_grid / grid_coarse
  integer,      intent(in)        ::  nend            ! number of levels
  integer,      intent(in)        ::  fct             ! 1: fg, 2: ana
  integer,      intent(in)        ::  diag_mode       ! 0: det run, 1: mean
  real(wp),     pointer, optional ::  xdiag(:,:,:,:,:)! diagnostic array (IN)
  integer,      intent(in),optional  ::  n_diag       ! number of elements
                                                      ! of xdiag
  real(wp),     pointer           :: pr(:)            ! press. values on
                                                      ! diagnostic grid (Pa)
  logical,           intent(in)   :: v_red            ! coarse grid used
  character(len=*)  ,intent(in)   :: par_stat ! parameters to read
  integer,           intent(in)   :: rf       !horizontal reducing factor
  integer,           intent(in)   :: nzr      !number of vertical levels
  integer,           intent(in)   :: norm_flag
  integer,           intent(in)   :: k_enkf
  character(len=*)  ,intent(in)   :: model    ! model: GME / COSMO
  character(len=*)  ,intent(in)   :: input    ! input path
  character(len=*)  ,intent(in)   :: output   ! output path
  type (t_time)     ,intent(in)   :: ana_time ! analysis time
    !----------------
    ! local variables
    !----------------
    integer                      ::  n_lon         ! number of longitudes
    integer                      ::  n_lat         ! number of latitudes
    integer                      ::  nzd           ! number of diagn. press. levels
    type (t_ctl)                 ::  ctl                 ! GRADS descriptor
    real(wp)                     ::  lo0, la0, dlo, dla  ! GRADS coordinates
    real(wp)                     ::  lonran,latran,lon_v,lat_v
    real(wp) , pointer           ::  rms_ges   (:,:,:,:) ! rms       +
    real(wp) , pointer           ::  spread_ges(:,:,:,:) ! spread    |
    real(wp) , pointer           ::  bias_ges  (:,:,:,:) ! bias      | on lon,lat,z grid
    real(wp) , pointer           ::  diag_ges  (:,:,:,:) ! grid-diag +
    real(wp) , pointer           ::  add_diag  (:)
    real(wp),  pointer           ::  d2        (:,:,:,:) ! layer thickness
!   real(wp),  pointer           ::  pht(:,:,:,:)        !half level pressure at top, interpolated
!   real(wp),  pointer           ::  qrs(:,:)
    integer ,  pointer           ::  n_ges(:,:,:)       ! number of grid points in lon,lat,z box (diagnostic)
    integer ,  pointer           ::  n_ges_ps(:,:)      ! number of grid points in lon,lat box (diagnostic, ps)
    integer ,  pointer           ::  indx (:,:,:,:)     ! GME/COSMO grid point --> index lon,lat box (diagnostic)
    integer ,  pointer           ::  indz (:,:,:,:)     ! GME/COSMO grid point --> z box (diagnostic)
    integer                      ::  m_max              ! max number of variables
    real(wp)                     ::  add, add_spread, add_bias
    real(wp)                     ::  delta_p     ! delta p in log(p)
    real(wp)                     ::  pmin        ! min. pressure
    real(wp)                     ::  pmax        ! max. pressure
    real(wp)                     ::  pmin_diag,pmax_diag
    integer                      ::  m2
    integer                      ::  i,j,n,d,m,k,l,index_lon,index_lat,index_z,m_ps,m_par
    character (len=5)            ::  names(10)
    character (len=100)          ::  filename
    real(wp) ,parameter          ::  eps = 1.e-6_wp ! fix rounding errors
    character(len=16)            :: pars(nm) ! names of variables to derive
    logical                      :: diag

    diag = .false. ; if(present(xdiag)) diag = .true.

    !--------------------------------
    !set par_diagn in enkf_fc for GME
    !--------------------------------
    !if (model=='GME') then
    !   call derive_params1   (enkf_mean,      par_stat, dealloc=.true.)
    !end if
    !and for COSMO??
    !not necessary???
    !set pressure
    if(enkf_grid% gridtype == DWD6_ICOSAHEDRON) then
       call set_p (enkf_mean)
    else if(enkf_grid% gridtype==WMO6_ROTLL .or. enkf_grid% gridtype==WMO6_LATLON) then
       call set_p (enkf_mean)
    end if

    !-------------------------------------
    ! calculate index-array for diagnostic
    !-------------------------------------
    select case(enkf_grid% gridtype)
    case(DWD6_ICOSAHEDRON)
       n_lat  = 36            ! number of latitude-bins (GME)
       n_lon  = 72            ! number of longitude-bins (GME)
       nzd    = 25            ! number of diagn. press. lev.
       pmax_diag = 105000._wp ! max press. on diagn. levels
       pmin_diag = 1000._wp   ! min press. on diagn. levels
       !determine diagnostic pressure levels
       if (dace% lpio) write(6,'(a)') ' derive diagnostic pressure levels'
       allocate (pr(nzd))
       pmax = maxval(enkf_mean% pf(:,:,enkf_grid% nz,:))
       pmax = p_max(pmax)
       pmax = max(pmax,pmax_diag) !check if true pmax is larger
       pmin = minval(enkf_mean% pf(:,:,1,:))
       pmin = p_min(pmin)
       pmin = min(pmin,pmin_diag) !check if true pmin is smaller
       delta_p = (log(pmax)-log(pmin))/(nzd-1)

       pr(1) = log(pmin)
       do i=1,nzd-1
          pr(i+1) = log(pmin) + i*delta_p
       end do
       pr = exp(pr)

!      if(dace% lpio) then
!         print *, 'diagnostic pressure levels pr = ', pr
!      end if

    case(WMO6_ROTLL,WMO6_LATLON)
       !n_lat  = 100            ! number of latitude-bins (COSMO)
       !n_lon  = 100            ! number of longitude-bins (COSMO)
       n_lon=enkf_grid% nx/rf   ! number of longitude-bins (COSMO)
       n_lat=enkf_grid% ny/rf   ! number of latitude-bins (COSMO)
       lonran = enkf_grid% dlon(enkf_grid% ubg(1))-enkf_grid% dlon(enkf_grid% lbg(1))
       latran = enkf_grid% dlat(enkf_grid% ubg(2))-enkf_grid% dlat(enkf_grid% lbg(2))
       !reduced_grid was not used or no vertical localization
       if (.not. v_red) then
          nzd    = 15            ! number of diagn. press. lev.
          pmax_diag = 105000._wp ! max press. on diagn. levels
          pmin_diag = 4000._wp   ! min press. on diagn. levels
          !determine diagnostic pressure levels
          if (dace% lpio) write(6,'(a)') ' derive diagnostic pressure levels'
          allocate (pr(nzd))
          pmax = maxval(enkf_mean% pf(:,:,enkf_grid% nz,:))
          pmax = p_max(pmax)
          pmax = max(pmax,pmax_diag) !check if true pmax is larger
          pmin = minval(enkf_mean% pf(:,:,1,:))
          pmin = p_min(pmin)
          pmin = min(pmin,pmin_diag) !check if true pmin is smaller
          delta_p = (log(pmax)-log(pmin))/(nzd-1)

          pr(1) = log(pmin)
          do i=1,nzd-1
             pr(i+1) = log(pmin) + i*delta_p
          end do
          pr = exp(pr)
       else
          !otherwise, pr was already computed in SR reduced_grid
          !use nzd = nzr (namelist variable)
          nzd  = nzr
          pmin = pr(1)
          pmax = pr(nzd)
          if(dace% lpio) then
             print *, 'pmin, pmax, nzd, nzr = ', pmin, pmax, nzd, nzr
          end if
       end if
!      if(dace% lpio) then
!         print *, 'diagnostic pressure levels pr = ', pr
!      end if
    case default
      write(6,'(a)') ' rms_statistic: grid not supported'
    end select


    allocate(n_ges_ps(n_lon,n_lat))
    n_ges_ps = 0

    allocate(n_ges(n_lon,n_lat,nzd))
    n_ges = 0

    allocate(indx(enkf_grid% lb(1):enkf_grid% ub(1), &
                  enkf_grid% lb(2):enkf_grid% ub(2), &
                  enkf_grid% lb(4):enkf_grid% ub(4),2))
    indx = 0

    allocate(indz(enkf_grid% lb(1):enkf_grid% ub(1), &
                  enkf_grid% lb(2):enkf_grid% ub(2), &
                  enkf_grid% lb(4):enkf_grid% ub(4), &
                  enkf_grid% lb(3):enkf_grid% ub(3)))
    indz = 0

    !loop over model grid to determine indx (horizontal)
!NEC$ nomove
    do d=enkf_grid% lb(4), enkf_grid% ub(4)
!NEC$ nomove
       do j=enkf_grid% lb(2), enkf_grid% ub(2)
!NEC$ nomove
          do i=enkf_grid% lb(1), enkf_grid% ub(1)
             select case(enkf_grid% gridtype)
             case(DWD6_ICOSAHEDRON)
                lat_v = enkf_grid% rlat(i,j,1,d)
                lon_v = enkf_grid% rlon(i,j,1,d)
                index_lat = ceiling(real(n_lat/pi,wp)*(lat_v+(pi/2._wp)))
                index_lon = ceiling(real(n_lon/(2.0*pi),wp)*(lon_v+pi))
             case(WMO6_ROTLL,WMO6_LATLON)
                lat_v = enkf_grid% dlat(j)
                lon_v = enkf_grid% dlon(i)
                index_lat = nint(real((n_lat-1)/latran,wp)*(lat_v+abs(enkf_grid% dlat(enkf_grid% lbg(2))))) + 1
                index_lon = nint(real((n_lon-1)/lonran,wp)*(lon_v+abs(enkf_grid% dlon(enkf_grid% lbg(1))))) + 1
             end select
             if (index_lat>n_lat) then
                index_lat = index_lat - 1
!               print *, 'i,j lat > ', i,j
             else if (index_lat==0) then
                index_lat = index_lat + 1
!               print *, 'i,j lat == ', i,j
             endif
             if (index_lon>n_lon) then
                index_lon = index_lon - 1
!               print *, 'i,j lon > ', i,j
             else if (index_lon==0) then
                index_lon = index_lon + 1
!               print *, 'i,j lon == ', i,j
             endif
             indx(i,j,d,1) = index_lon
             indx(i,j,d,2) = index_lat
             n_ges_ps(index_lon,index_lat) = n_ges_ps(index_lon,index_lat) + 1
             !now vertical loop
             do k=enkf_grid% lb(3), enkf_grid% ub(3)
                !index_z = nint(real((nzd-1)/(log(pmax)-log(pmin)),wp)*(log(enkf_mean% pf(i,j,k,d))-log(pmin))) + 1
                index_z = sum(minloc(abs(enkf_mean% pf(i,j,k,d)-pr)))
                if (index_z>nzd) then
                   index_z = index_z - 1
                else if (index_z==0) then
                   index_z = index_z + 1
                endif
                indz(i,j,d,k) = index_z
                n_ges(index_lon,index_lat,index_z) = n_ges(index_lon,index_lat,index_z) + 1
             end do
          end do
       end do
    end do
!   print *, 'finish index'
!   print *, 'lb,ub(indx,1,2) in letkf, proc =', lbound(indx,1),ubound(indx,1), &
!        lbound(indx,2),ubound(indx,2),dace% pe,indx(lbound(indx,1),lbound(indx,2),1,1)
    n_ges    = p_sum(n_ges)
    n_ges_ps = p_sum(n_ges_ps)

    !---------------------------------
    ! allocate arrays for regular grid
    !---------------------------------
    m_max = size(analysis% m) ! = nm = 62?
    call split (pars, par_stat, n)
    m_par = n !n = number of variables in par_stat
!   print *, 'n = ', n
    allocate (rms_ges    (n_lon, n_lat, m_par, nzd)) !nend
    allocate (spread_ges (n_lon, n_lat, m_par, nzd))
    allocate (bias_ges   (n_lon, n_lat, m_par, nzd))
    if (diag) &
      allocate (diag_ges   (n_lon, n_lat, n_diag, nzd))
    if (model == 'GME' .and. diag) then !more general: model grid = letkf grid
       allocate (add_diag   (n_diag))
    end if
    allocate (d2(enkf_grid% lb(1):enkf_grid% ub(1),&
                 enkf_grid% lb(2):enkf_grid% ub(2),&
                 enkf_grid% lb(3):enkf_grid% ub(3),&
                 enkf_grid% lb(4):enkf_grid% ub(4)))
    d2 = 0 ! Need to fix code below!

!   print *, 'min indx,n_ges = ', minval(indx), minval(n_ges)
!   print *, 'max indx,n_ges = ', maxval(indx), maxval(n_ges)
!   print *, 'min indz,max indz = ', minval(indz), maxval(indz)
!   print *, 'n_lon,n_lat = ', n_lon, n_lat
!   print *, 'lb,ub(indx,1,2) in rms, proc =', lbound(indx,1),ubound(indx,1), &
!        lbound(indx,2),ubound(indx,2),dace% pe,indx(lbound(indx,1),lbound(indx,2),1,1)
!   print *, 'm_max = ', m_max
!   print *, 'n = number of vars = ', n
!   if (diag) &
!     print *, 'n_diag = ', n_diag
    !print *, 'rms stat: lbound(xdiag,1),ubound(xdiag,1),lbound(xdiag,4),ubound(xdiag,4), dace% pe = ', &
    !                    lbound(xdiag,1),ubound(xdiag,1),lbound(xdiag,4),ubound(xdiag,4), dace% pe
    !------------------------
    ! compute layer thickness
    !------------------------

    call print (enkf_mean, comment='rms stat: enkf_mean')
    !call print (enkf_ens(1), comment='rms stat: enkf_ens(1)')
    call print (analysis, comment='rms stat: analysis')

!!$    if (enkf_grid% ivctype /= 0) then
!!$       ie = size(enkf_mean% pp,1)
!!$       je = size(enkf_mean% pp,2)
!!$       ke = size(enkf_mean% pp,3)
!!$       allocate (qrs (ie,je))
!!$       allocate (pht(enkf_grid% lb(1):enkf_grid% ub(1),&
!!$                 enkf_grid% lb(2):enkf_grid% ub(2),&
!!$                 1:1,&
!!$                 enkf_grid% lb(4):enkf_grid% ub(4)))
!!$       qrs = 0._wp
!!$       if (associated(enkf_mean% qr)) qrs = &
!!$            qrs + enkf_mean% qr(:,:,ke,1) ! rain
!!$       if (associated(enkf_mean% qs)) qrs = &
!!$            qrs + enkf_mean% qs(:,:,ke,1) ! snow
!!$       ! ice,graupel?
!!$       call calps(pht(:,:,1,1),enkf_mean% pp(:,:,ke,1),&
!!$            enkf_mean% t(:,:,ke,1),  &
!!$            enkf_mean% q(:,:,ke,1),  &
!!$            enkf_mean% qcl(:,:,ke,1),&
!!$            qrs,                        &
!!$            enkf_grid% rho0(:,:,ke,1),  &
!!$            enkf_grid% p0(:,:,ke,1),    &
!!$            enkf_grid% dp0(:,:,ke,1),   &
!!$            ie,je,                      &
!!$            rddrm1,                     &
!!$            R_c,                        &
!!$            1,ie,1,je)
!!$    end if

!!$    do d=enkf_grid% lb(4), enkf_grid% ub(4)
!!$       do j=enkf_grid% lb(2), enkf_grid% ub(2)
!!$          do i=enkf_grid% lb(1), enkf_grid% ub(1)
!!$             do n=nend,1,-1
!!$                select case(enkf_grid% ivctype)
!!$                case(0)
!!$                   d2(i,j,n,d) = ((enkf_mean% ph (i,j,n+1,d)&
!!$                        -enkf_mean% ph (i,j,n,d))&
!!$                        / enkf_mean% ps (i,j,1,d))
!!$                case default
!!$                   if(n > 1 .and. n < nend) then
!!$                      d2(i,j,n,d) = ((((enkf_grid% dp0(i,j,n+1,d)-enkf_grid% dp0(i,j,n-1,d))&
!!$                           *enkf_mean% pf(i,j,n,d)*enkf_grid% dp0(i,j,n,d))/&
!!$                           ((enkf_grid% dp0(i,j,n+1,d)+enkf_grid% dp0(i,j,n,d))&
!!$                            *(enkf_grid% dp0(i,j,n-1,d)+enkf_grid% dp0(i,j,n,d))))+&
!!$                            (enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n+1,d)/&
!!$                            (enkf_grid% dp0(i,j,n+1,d)+enkf_grid% dp0(i,j,n,d)))-&
!!$                            (enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n-1,d)/&
!!$                            (enkf_grid% dp0(i,j,n-1,d)+enkf_grid% dp0(i,j,n,d))))/enkf_mean% ps(i,j,1,d)
!!$                   else if(n==nend) then
!!$                      d2(i,j,n,d) = (enkf_mean% ps(i,j,1,d) - ((enkf_grid% dp0(i,j,n-1,d)*enkf_mean% pf(i,j,n,d)+&
!!$                                     enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n-1,d))/&
!!$                                     (enkf_grid% dp0(i,j,n-1,d)+enkf_grid% dp0(i,j,n,d)))) &
!!$                                     /enkf_mean% ps(i,j,1,d)
!!$                   else !(n==1)
!!$                      !d2(i,j,n,d) = d2(i,j,n+1,d)
!!$                       d2(i,j,n,d) = (((enkf_grid% dp0(i,j,n+1,d)*enkf_mean% pf(i,j,n,d)+&
!!$                                     enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n+1,d))/&
!!$                                     (enkf_grid% dp0(i,j,n+1,d)+enkf_grid% dp0(i,j,n,d)))-pht(i,j,1,1))/&
!!$                                     enkf_mean% ps(i,j,1,d)
!!$                       !(ph(k+1/2)-pht)/ps
!!$                   end if
!!$                end select
!!$             end do
!!$          end do
!!$       end do
!!$    end do

    !------------------
    ! sum up statistics
    !------------------
    m2         = 0
    rms_ges    = 0._wp
    spread_ges = 0._wp
    bias_ges   = 0._wp
    if (diag) &
      diag_ges   = 0._wp

!   print *, 'local ub(j)', enkf_grid% ub(2)
    do d=enkf_grid% lb(4), enkf_grid% ub(4)
       do j=enkf_grid% lb(2), enkf_grid% ub(2)
          do i=enkf_grid% lb(1), enkf_grid% ub(1)
             do n=1, nend
                !-------------
                !fill diag_ges
                !-------------
                if (model=='GME' .and. diag) then !more general: model grid = letkf grid
                   add_diag = xdiag(i,j,n,d,:)
                   index_lon = indx(i,j,d,1)
                   index_lat = indx(i,j,d,2)
                   index_z   = indz(i,j,d,n)
                   diag_ges(index_lon,index_lat,:,index_z) = &
                            diag_ges(index_lon,index_lat,:,index_z) + add_diag
                end if
                m2 = 0
                do m=1, m_max
                   if (enkf_mean% m(m)% i% alloc .and.     &
                        analysis% m(m)% i% alloc .and.         &
                        .not. enkf_mean% m(m)% i% ref .and. &
                        .not. analysis% m(m)% i% ref           ) then
                      m2 = m2 + 1
                      names(m2) = enkf_mean% m(m)% i% name
                      select case (enkf_mean% m(m)% i% name)
                      !-----------------
                      ! surface pressure
                      !-----------------
                      case ('ps')
                         if (n==nend) then
                            add_spread = 0._wp
                            select case (norm_flag)
                            case (1)
                               do k=1, k_enkf
                                  add_spread = add_spread &
                                       +(log(enkf_ens(k)%  m(m)% ptr(i,j,1,d))&
                                       -log(enkf_mean%  m(m)% ptr(i,j,1,d)))**2
                               end do
                               add_spread = add_spread * R_c * t0c
                            case (0)
                               !do k=1, k_enkf
                               !   add_spread = add_spread &
                               !        +(enkf_ens(k)%  m(m)% ptr(i,j,1,d)&
                               !        -enkf_mean%  m(m)% ptr(i,j,1,d))**2
                               !end do
                               add_spread = (enkf_spread%  m(m)% ptr(i,j,1,d))**2
                            end select
                            !add_spread = add_spread*(1._wp/(k_enkf-1))
                            select case (norm_flag)
                            case (1)
                               add = (log(enkf_mean%  m(m)% ptr(i,j,1,d))&
                                    -log(analysis%  m(m)% ptr(i,j,1,d)))**2
                               add_bias = (log(enkf_mean%  m(m)% ptr(i,j,1,d))&
                                    -log(analysis%  m(m)% ptr(i,j,1,d)))
                               add = add * R_c * t0c
                               add_bias = add_bias * R_c * t0c
                            case (0)
                               add = (enkf_mean%  m(m)% ptr(i,j,1,d)&
                                    -analysis%  m(m)% ptr(i,j,1,d))**2
                               add_bias = (enkf_mean%  m(m)% ptr(i,j,1,d)&
                                    -analysis%  m(m)% ptr(i,j,1,d))
                            end select
                            index_lon = indx(i,j,d,1)
                            index_lat = indx(i,j,d,2)
                            index_z   = nzd !indz(i,j,d,n) store ps in pmax-level
                            !if (d==5 .and. j==5 .and. i==5) then
                            !   print *, enkf_mean% m(m)% i% name, enkf_ens(12)% m(m)% i% name, analysis% m(m)% i% name, m2
                            !end if
                            rms_ges(index_lon,index_lat,m2,nzd) = &
                                 rms_ges(index_lon,index_lat,m2,nzd) + add
                            spread_ges(index_lon,index_lat,m2,nzd) = &
                                 spread_ges(index_lon,index_lat,m2,nzd) + add_spread
                            bias_ges(index_lon,index_lat,m2,nzd) = &
                                 bias_ges(index_lon,index_lat,m2,nzd) + add_bias
                         end if
                      case default
                         add_spread = 0._wp
                         select case (norm_flag)
                         case (1)
                            do k=1, k_enkf
                               add_spread = add_spread &
                                    +d2(i,j,n,d)*(enkf_ens(k)%  m(m)% ptr(i,j,n,d)&
                                    -enkf_mean%  m(m)% ptr(i,j,n,d))**2
                            end do
                            select case (enkf_mean% m(m)% i% name)
                            case ('t')
                               add_spread = add_spread * (c_p/t0c)
                            end select
                         case (0)
                            !do k=1, k_enkf
                            !   add_spread = add_spread &
                            !        +(enkf_ens(k)%  m(m)% ptr(i,j,n,d)&
                            !        -enkf_mean%  m(m)% ptr(i,j,n,d))**2
                            !end do
                            add_spread = (enkf_spread%  m(m)% ptr(i,j,n,d))**2
                         end select
                         !add_spread = add_spread*(1._wp/(k_enkf-1))
                         select case (norm_flag)
                         case (1)
                            add = d2(i,j,n,d)*(enkf_mean%  m(m)% ptr(i,j,n,d)&
                                 - analysis%  m(m)% ptr(i,j,n,d))**2
                            add_bias = d2(i,j,n,d)*(enkf_mean%  m(m)% ptr(i,j,n,d)&
                                 - analysis%  m(m)% ptr(i,j,n,d))
                            select case (enkf_mean% m(m)% i% name)
                            case ('t')
                               add = add * (c_p/t0c)
                               add_bias = add_bias * (c_p/t0c)
                            end select
                         case (0)
                            add = (enkf_mean%  m(m)% ptr(i,j,n,d)&
                                 - analysis%  m(m)% ptr(i,j,n,d))**2
                            add_bias = (enkf_mean%  m(m)% ptr(i,j,n,d)&
                                 - analysis%  m(m)% ptr(i,j,n,d))
                         end select
                         index_lon = indx(i,j,d,1)
                         index_lat = indx(i,j,d,2)
                         index_z   = indz(i,j,d,n)
                         !if (d==1 .and. j==5 .and. i==5 .and. n==1) then
                         !   print *, 'names of vars = ', &
                         !        enkf_mean% m(m)% i% name, enkf_ens(1)% m(m)% i% name, analysis% m(m)% i% name, m2
                         !end if
                         rms_ges(index_lon,index_lat,m2,index_z) = &
                              rms_ges(index_lon,index_lat,m2,index_z) + add
                         spread_ges(index_lon,index_lat,m2,index_z) = &
                              spread_ges(index_lon,index_lat,m2,index_z) + add_spread
                         bias_ges(index_lon,index_lat,m2,index_z) = &
                              bias_ges(index_lon,index_lat,m2,index_z) + add_bias
                      end select
                   endif
                end do !m-loop
             end do !n-loop
          end do !i-loop
       end do !j-loop
    end do !d-loop


    !----------------------------
    ! sum over processor elements
    !----------------------------
    rms_ges    = p_sum (rms_ges)
    spread_ges = p_sum (spread_ges)
    bias_ges   = p_sum (bias_ges)
    if (diag) &
      diag_ges   = p_sum (diag_ges)

    !----------
    ! normalize
    !----------
    !determine variable index of surface pressure ps (GME only)
    if (model == 'GME') then
       do i = 1,m_max
          if (names(i) == 'ps') exit
       end do
       m_ps = i
!      print *, 'm_ps = ', m_ps
    end if

!   print *, 'names = ', names

    do k=1,n_lon
       do l=1,n_lat
          if (model == 'GME') then
             if (n_ges_ps  (k,l) .ne. 0) then
                rms_ges    (k,l,m_ps,nzd) = rms_ges    (k,l,m_ps,nzd) / n_ges_ps(k,l)
                spread_ges (k,l,m_ps,nzd) = spread_ges (k,l,m_ps,nzd) / n_ges_ps(k,l)
                bias_ges   (k,l,m_ps,nzd) = bias_ges   (k,l,m_ps,nzd) / n_ges_ps(k,l)
             else
                !set to undef value
                rms_ges    (k,l,m_ps,nzd) = -huge(1._sp)
                spread_ges (k,l,m_ps,nzd) = -huge(1._sp)
                bias_ges   (k,l,m_ps,nzd) = -huge(1._sp)
             end if
          end if
          do n=1,nzd
             if (n_ges     (k,l,n) .ne. 0 .and. diag ) then
                diag_ges   (k,l,:,n)  = diag_ges   (k,l,:,n) / n_ges(k,l,n)
             end if
             do m=1,m_par
                if (names(m) == 'ps') cycle
                if (n_ges     (k,l,n) .ne. 0) then
                   rms_ges    (k,l,m,n) = rms_ges    (k,l,m,n) / n_ges(k,l,n)
                   spread_ges (k,l,m,n) = spread_ges (k,l,m,n) / n_ges(k,l,n)
                   bias_ges   (k,l,m,n) = bias_ges   (k,l,m,n) / n_ges(k,l,n)
                else
                   !set to undef value
                   rms_ges    (k,l,m,n) = -huge(1._sp)
                   spread_ges (k,l,m,n) = -huge(1._sp)
                   bias_ges   (k,l,m,n) = -huge(1._sp)
                end if
             end do
          end do
       end do
    end do

    !--------------------------------------------
    ! derive coordinates of coarse grid for GRADS
    !--------------------------------------------
    if (dace% lpio) then
      select case(enkf_grid% gridtype)
      case(DWD6_ICOSAHEDRON)
        lo0 = -180._wp + 180._wp / n_lon
        la0 =  -90._wp +  90._wp / n_lat
        dlo =            360._wp / n_lon
        dla =            180._wp / n_lat
      case(WMO6_ROTLL,WMO6_LATLON)
        lo0 = grid_weight% lo1 !-4.9475_wp
        la0 = grid_weight% la1 !-4.9425_wp
        dlo = grid_weight% di  !0.105_wp
        dla = grid_weight% dj  !0.115_wp
      end select
      !------------------------
      ! write rms to GRADS-file
      !------------------------
      filename = path_file (output,'^rms_ges__ANA_TIME_')
      call init_ctl (ctl, file=filename, nx=n_lon, ny=n_lat, ke=nzd,&
                     zlev=pr/100._wp, lo1=lo0, di=dlo,la1=la0, dj=dla)
!     print *, 'ke = ', nzd
!     print *, 'm_par, names(1), names(m_par)', m_par, names(1), names(m_par)
!     print *, 'filename = ', filename
!     print *, 'output = ', output
      !print *, 'ana_time = ', ana_time% days, ana_time% secs
      do i=1,m_par
         if (names(i) == 'ps') then
            call add_var (ctl,names(i),1,names(i))
         else
            call add_var (ctl,names(i),nzd,names(i))
         end if
         if (names(i) == 'ps') then
            call write_var (ctl,rms_ges(:,:,i,nzd),names(i))
         else
            call write_var (ctl,rms_ges(:,:,i,:),names(i))
         end if
      end do
      call write_ctl (ctl)
      call destruct(ctl)
      !---------------------------
      ! write spread to GRADS file
      !---------------------------
      filename = path_file (output,'^spread_ges__ANA_TIME_')
      call init_ctl (ctl, file=filename, nx=n_lon, ny=n_lat, ke=nzd,&
                     zlev=pr/100._wp, lo1=lo0, di=dlo,la1=la0, dj=dla)
      do i=1,m_par
         if (names(i) == 'ps') then
            call add_var (ctl,names(i),1,names(i))
         else
            call add_var (ctl,names(i),nzd,names(i))
         end if
         if (names(i) == 'ps') then
            call write_var (ctl,spread_ges(:,:,i,nzd),names(i))
         else
            call write_var (ctl,spread_ges(:,:,i,:),names(i))
         end if
      end do
      call write_ctl (ctl)
      call destruct(ctl)
      !-------------------------
      ! write bias to grads file
      !-------------------------
      filename = path_file (output,'^bias_ges__ANA_TIME_')
      call init_ctl (ctl, file=filename, nx=n_lon, ny=n_lat, ke=nzd,&
                     zlev=pr/100._wp, lo1=lo0, di=dlo,la1=la0, dj=dla)
      do i=1,m_par
         if (names(i) == 'ps') then
            call add_var (ctl,names(i),1,names(i))
         else
            call add_var (ctl,names(i),nzd,names(i))
         end if
         if (names(i) == 'ps') then
            call write_var (ctl,bias_ges(:,:,i,nzd),names(i))
         else
            call write_var (ctl,bias_ges(:,:,i,:),names(i))
         end if
      end do
      call write_ctl (ctl)
      call destruct(ctl)

      !-------------------------
      ! write diag to grads file
      !-------------------------
      if (model == 'GME' .and. diag) then
         filename = path_file (output,'^diag_ges__ANA_TIME_')
         call init_ctl (ctl, file=filename, nx=n_lon, ny=n_lat, ke=nzd,&
                        zlev=pr/100._wp, lo1=lo0, di=dlo,la1=la0, dj=dla)
         call add_var (ctl,'nobs',nzd,'nobs')
         call add_var (ctl,'obswq',nzd,'obswq')
         call add_var (ctl,'x2',nzd,'x2')
         call add_var (ctl,'edimz',nzd,'edimz')
         call add_var (ctl,'edimn',nzd,'edimn')
         call add_var (ctl,'maxwens',nzd,'maxwens')
         call add_var (ctl,'rhoadap',nzd,'rhoadap')
         call add_var (ctl,'rholoc',nzd,'rholoc')
         call add_var (ctl,'dob_sum',nzd,'dob_sum')
         call add_var (ctl,'HPH_sum',nzd,'HPH_sum')
         call add_var (ctl,'bc_weight',nzd,'bc_weight')
         call add_var (ctl,'x2_adap',nzd,'x2_adap')
         call write_var (ctl,diag_ges(:,:,1,:),'nobs')
         call write_var (ctl,diag_ges(:,:,2,:),'obswq')
         call write_var (ctl,diag_ges(:,:,3,:),'x2')
         call write_var (ctl,diag_ges(:,:,4,:),'edimz')
         call write_var (ctl,diag_ges(:,:,5,:),'edimn')
         call write_var (ctl,diag_ges(:,:,6,:),'maxwens')
         call write_var (ctl,diag_ges(:,:,7,:),'rhoadap')
         call write_var (ctl,diag_ges(:,:,8,:),'rholoc')
         call write_var (ctl,diag_ges(:,:,9,:),'dob_sum')
         call write_var (ctl,diag_ges(:,:,10,:),'HPH_sum')
         call write_var (ctl,diag_ges(:,:,11,:),'bc_weight')
         call write_var (ctl,diag_ges(:,:,12,:),'x2_adap')
         call write_ctl (ctl)
         call destruct (ctl)
      end if
   endif
   deallocate(n_ges)
   deallocate(n_ges_ps)
   deallocate(indx)
   deallocate(indz)

  end subroutine rms_statistics

!=============================================================================
!!$  subroutine ens_statistics (enkf_mean, enkf_ens, analysis, enkf_grid, &
!!$                             nend, par_stat, rf, nzr, k_enkf, model)
!!$  !--------------------------------------------------------------------------
!!$  ! Calculates explained variance, ensemble dimension etc. on a set
!!$  ! of model levels; output as GRADS readable (ieee) files + GRADS
!!$  ! descriptor (.ctl) files on a horizontal coarsed grid.
!!$  !--------------------------------------------------------------------------
!!$
!!$  type (t_atm), intent(inout)     ::  enkf_mean       ! LETKF mean   (fg/ana)
!!$  type (t_atm), intent(in)        ::  enkf_ens (:)    ! LETKF ensemble (fg/ana)
!!$  type (t_atm), intent(in)        ::  analysis        ! 'truth'
!!$  type (t_grid),intent(in)        ::  enkf_grid       ! reduced grid
!!$  integer,      intent(in)        ::  nend            ! number of levels
!!$  character(len=*)  ,intent(in)   :: par_stat         ! parameters to read
!!$  integer,           intent(in)   :: k_enkf           ! number of ens members
!!$  character(len=*)  ,intent(in)   :: model            ! model: GME / COSMO
!!$  integer,           intent(in)   :: rf               ! horiz. reduc. factor (COSMO only)
!!$
!!$  !---------------
!!$  !Local variables
!!$  !---------------
!!$    real(wp) , pointer     ::  expl_var(:,:,:,:)   ! explained variance
!!$    real(wp) , pointer     ::  e_dim(:,:,:,:)      ! ensemble-dimension
!!$    real(wp) , pointer     ::  ml(:,:,:,:)         ! max. likelihood
!!$    real(wp) , pointer     ::  rms_e(:,:,:,:)      ! rms in energy-norm
!!$    real(wp) , pointer     ::  equ(:,:,:,:)        ! E^(-1/2)*Q^t*(x_b-x_t)
!!$    real(wp) , pointer     ::  equ_2(:,:,:,:)      ! square of elements of equ
!!$    real(wp) , pointer     ::  ew(:,:,:,:)         ! eigenvalues
!!$    real(wp)               ::  pi, add, add2, add_spread, add_bias, add_bias2, ev_add, edim_add, ml_add, rms_e_add
!!$    real(wp),  pointer     ::  equ_add(:)          ! holds column of equ at certain grid point
!!$    real(wp),  pointer     ::  ew_add(:)           ! holds eigenvalues at certain grid point
!!$    real(wp) , pointer     ::  lat(:)
!!$    real(wp) , pointer     ::  lon(:)
!!$    real(wp)               ::  lonran,latran
!!$    real(wp)               ::  lat_v
!!$    real(wp)               ::  lon_v
!!$    real(wp),  pointer     ::  v_xb(:)        ! contains (x_b-x_t)
!!$    real(wp),  pointer     ::  m_xb(:,:)      ! contains X_b
!!$    real(wp),  pointer     ::  v_xb_nps(:)    ! contains (x_b-x_t) (levels without ps)
!!$    real(wp),  pointer     ::  m_xb_nps(:,:)  ! contains X_b       (levels without ps)
!!$    real(wp),  pointer     ::  p_t(:,:)       ! contains (X_b)^t(X_b)
!!$    real(wp),  pointer     ::  u(:,:)         ! contains eigenvec of p_t
!!$    real(wp),  pointer     ::  e(:)           ! contains eigenval of p_t
!!$    real(wp),  pointer     ::  e_i(:)         ! contains inverse eigenvals of p_t
!!$    real(wp),  pointer     ::  u_qt(:)        ! contains Q^t*U
!!$    integer                ::  n_n, n_v       ! number of neighbour grid points (total) (hor./ver.)
!!$    integer                ::  n_i, n_v_i     ! number of neighbour grid points (one direction) (hor./ver.)
!!$    integer                ::  n_grid         ! number of grid points in loc vol (total)
!!$    integer                ::  n_grid_m       ! number of grid points in loc vol (per variable)
!!$    integer                ::  n_grid_nps     !
!!$    integer                ::  n_lat
!!$    integer                ::  n_lon
!!$    integer                ::  n_c            ! correction for ps
!!$    integer                ::  ps_flag        ! flag for surface pressure
!!$    integer , pointer      ::  n_ges(:,:)
!!$    integer , pointer      ::  n_ges_e(:,:)
!!$    integer , pointer      ::  index(:,:,:,:)
!!$    integer                ::  m_max
!!$    integer                ::  nh
!!$    integer                ::  m2,m_v,n_z    ! counts variables in rms/exp.var calculations or levels
!!$    integer                ::  i,j,n,d,k,l,m,ie,je,ke !(loop) indices
!!$    integer                ::  index_lon,index_lat,ni,nj,nd
!!$    logical, pointer       ::  mask(:,:,:)
!!$    character (len=100)    ::  filename
!!$    real(wp), allocatable  ::  d2(:,:,:,:)
!!$    real(wp),  pointer     ::  pht(:,:,:,:)        !half level pressure at top, interpolated
!!$    real(wp),  pointer     ::  qrs(:,:)
!!$    type (t_grid) ,pointer :: g
!!$    real(wp), parameter    :: R = 287.0
!!$    real(wp), parameter    :: T = 273.0
!!$    real(wp), parameter    :: C_p = 1004.0
!!$    type (t_ctl)           :: ctl
!!$
!!$    !-------------
!!$    ! choose model
!!$    !-------------
!!$    select case(enkf_grid% gridtype)
!!$    case(DWD6_ICOSAHEDRON)
!!$       n_lat  = 36            ! number of latitude-bins (GME)
!!$       n_lon  = 72            ! number of longitude-bins (GME)
!!$    case(WMO6_ROTLL,WMO6_LATLON)
!!$       n_lon=enkf_grid% nx/rf   ! number of longitude-bins (COSMO)
!!$       n_lat=enkf_grid% ny/rf   ! number of latitude-bins (COSMO)
!!$       lonran = enkf_grid% dlon(enkf_grid% ubg(1))-enkf_grid% dlon(enkf_grid% lbg(1))
!!$       latran = enkf_grid% dlat(enkf_grid% ubg(2))-enkf_grid% dlat(enkf_grid% lbg(2))
!!$    end select
!!$
!!$
!!$    !-------------
!!$    ! set values
!!$    !-------------
!!$    pi = 3.14159265359
!!$    n_n = 5                                   ! number of hor. neighbours
!!$    n_i = (n_n-1)/2                           ! neighb. in  hor. direction
!!$    n_v = 11                                  ! number of ver. neighbours
!!$    n_v_i = (n_v-1)/2                         ! neighb. in  vert. direction
!!$    m_max = max(size(enkf_mean1% m), size(analysis% m))
!!$    n_grid = ((n_n*n_n)-6)*(3*n_v+1)
!!$    n_grid_nps = ((n_n*n_n)-6)*(3*n_v)
!!$    n_grid_m = ((n_n*n_n)-6)*n_v
!!$
!!$    !------------------------------
!!$    ! allocate diagnostic variables
!!$    !------------------------------
!!$
!!$    allocate (expl_var(n_lon,n_lat,1,4))
!!$    allocate (e_dim(n_lon,n_lat,1,4))
!!$    allocate (equ(n_lon,n_lat,6,k_enkf-1))
!!$    allocate (equ_2(n_lon,n_lat,3,k_enkf-1))
!!$    allocate (ew(n_lon,n_lat,3,k_enkf-1))
!!$    allocate (equ_add(k_enkf-1))
!!$    allocate (ew_add(k_enkf-1))
!!$    allocate (ml(n_lon,n_lat,1,4))
!!$    allocate (rms_e(n_lon,n_lat,1,4))
!!$    allocate (v_xb(n_grid))
!!$    allocate (v_xb_nps(n_grid_nps))
!!$    allocate (m_xb(n_grid,k_enkf))
!!$    allocate (m_xb_nps(n_grid_nps,k_enkf))
!!$    allocate (p_t(k_enkf,k_enkf))
!!$    allocate (u(k_enkf,k_enkf))
!!$    allocate (e(k_enkf))
!!$    allocate (e_i(k_enkf-1))
!!$    allocate (u_qt(k_enkf-1))
!!$    allocate (mask(n_n,n_n,n_v))
!!$
!!$    !------------------------
!!$    !Set pressure half levels
!!$    !------------------------
!!$    call set_ph (enkf_mean)
!!$
!!$    !-------------------------------------
!!$    ! calculate index-array for diagnostic
!!$    !-------------------------------------
!!$    allocate(n_ges(n_lon,n_lat))
!!$    allocate(n_ges_e(n_lon,n_lat))
!!$
!!$    allocate(index(enkf_grid% lb(1):enkf_grid% ub(1), &
!!$                  enkf_grid% lb(2):enkf_grid% ub(2), &
!!$                  enkf_grid% lb(4):enkf_grid% ub(4),2))
!!$    index   = 0
!!$    n_ges   = 0
!!$    n_ges_e = 0
!!$
!!$    !loop over model grid to determine indx (horizontal)
!!$!CDIR NOMOVE
!!$    do d=enkf_grid% lb(4), enkf_grid% ub(4)
!!$!CDIR NOMOVE
!!$       do j=enkf_grid% lb(2), enkf_grid% ub(2)
!!$!CDIR NOMOVE
!!$          do i=enkf_grid% lb(1), enkf_grid% ub(1)
!!$             select case(enkf_grid% gridtype)
!!$             case(DWD6_ICOSAHEDRON)
!!$                lat_v = enkf_grid% rlat(i,j,1,d)
!!$                lon_v = enkf_grid% rlon(i,j,1,d)
!!$                index_lat = ceiling(real(n_lat/pi,wp)*(lat_v+(pi/2._wp)))
!!$                index_lon = ceiling(real(n_lon/(2.0*pi),wp)*(lon_v+pi))
!!$             case(WMO6_ROTLL,WMO6_LATLON)
!!$                lat_v = enkf_grid% dlat(j)
!!$                lon_v = enkf_grid% dlon(i)
!!$                index_lat = nint(real((n_lat-1)/latran,wp)*(lat_v+abs(enkf_grid% dlat(enkf_grid% lbg(2))))) + 1
!!$                index_lon = nint(real((n_lon-1)/lonran,wp)*(lon_v+abs(enkf_grid% dlon(enkf_grid% lbg(1))))) + 1
!!$             end select
!!$             if (index_lat>n_lat) then
!!$                index_lat = index_lat - 1
!!$!               print *, 'i,j lat > ', i,j
!!$             else if (index_lat==0) then
!!$                index_lat = index_lat + 1
!!$!               print *, 'i,j lat == ', i,j
!!$             endif
!!$             if (index_lon>n_lon) then
!!$                index_lon = index_lon - 1
!!$!               print *, 'i,j lon > ', i,j
!!$             else if (index_lon==0) then
!!$                index_lon = index_lon + 1
!!$!               print *, 'i,j lon == ', i,j
!!$             endif
!!$             index(i,j,d,1) = index_lon
!!$             index(i,j,d,2) = index_lat
!!$             n_ges(index_lon,index_lat) = n_ges(index_lon,index_lat) + 1
!!$             if (j>(enkf_grid% lb(2)+1) .and. i>(enkf_grid% lb(1)+1) &
!!$                .and. j<(enkf_grid% ub(2)-1) .and. i<(enkf_grid% ub(1)-1)) then
!!$                n_ges_e(index_lon,index_lat) = n_ges_e(index_lon,index_lat)+1
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    n_ges = p_sum(n_ges)
!!$    n_ges_e = p_sum(n_ges_e)
!!$
!!$    !*************
!!$     !------------------------
!!$    ! compute layer thickness
!!$    !------------------------
!!$    allocate (d2(enkf_grid% lb(1):enkf_grid% ub(1),&
!!$                 enkf_grid% lb(2):enkf_grid% ub(2),&
!!$                 enkf_grid% lb(3):enkf_grid% ub(3),&
!!$                 enkf_grid% lb(4):enkf_grid% ub(4)))
!!$
!!$    if (enkf_grid% ivctype /= 0) then
!!$       ie = size(enkf_mean% pp,1)
!!$       je = size(enkf_mean% pp,2)
!!$       ke = size(enkf_mean% pp,3)
!!$       allocate (qrs (ie,je))
!!$       allocate (pht(enkf_grid% lb(1):enkf_grid% ub(1),&
!!$                 enkf_grid% lb(2):enkf_grid% ub(2),&
!!$                 1:1,&
!!$                 enkf_grid% lb(4):enkf_grid% ub(4)))
!!$       qrs = 0._wp
!!$       if (associated(enkf_mean% qr)) qrs = &
!!$            qrs + enkf_mean% qr(:,:,ke,1) ! rain
!!$       if (associated(enkf_mean% qs)) qrs = &
!!$            qrs + enkf_mean% qs(:,:,ke,1) ! snow
!!$       ! ice,graupel?
!!$       call calps(pht(:,:,1,1),enkf_mean% pp(:,:,ke,1),&
!!$            enkf_mean% t(:,:,ke,1),  &
!!$            enkf_mean% q(:,:,ke,1),  &
!!$            enkf_mean% qcl(:,:,ke,1),&
!!$            qrs,                        &
!!$            enkf_grid% rho0(:,:,ke,1),  &
!!$            enkf_grid% p0(:,:,ke,1),    &
!!$            enkf_grid% dp0(:,:,ke,1),   &
!!$            ie,je,                      &
!!$            rddrm1,                     &
!!$            R_c,                        &
!!$            1,ie,1,je)
!!$    end if

!!$    do d=enkf_grid% lb(4), enkf_grid% ub(4)
!!$       do j=enkf_grid% lb(2), enkf_grid% ub(2)
!!$          do i=enkf_grid% lb(1), enkf_grid% ub(1)
!!$             do n=nend,1,-1
!!$                select case(enkf_grid% ivctype)
!!$                case(0)
!!$                   d2(i,j,n,d) = ((enkf_mean% ph (i,j,n+1,d)&
!!$                        -enkf_mean% ph (i,j,n,d))&
!!$                        / enkf_mean% ps (i,j,1,d))
!!$                case default
!!$                   if(n > 1 .and. n < nend) then
!!$                      d2(i,j,n,d) = ((((enkf_grid% dp0(i,j,n+1,d)-enkf_grid% dp0(i,j,n-1,d))&
!!$                           *enkf_mean% pf(i,j,n,d)*enkf_grid% dp0(i,j,n,d))/&
!!$                           ((enkf_grid% dp0(i,j,n+1,d)+enkf_grid% dp0(i,j,n,d))&
!!$                            *(enkf_grid% dp0(i,j,n-1,d)+enkf_grid% dp0(i,j,n,d))))+&
!!$                            (enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n+1,d)/&
!!$                            (enkf_grid% dp0(i,j,n+1,d)+enkf_grid% dp0(i,j,n,d)))-&
!!$                            (enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n-1,d)/&
!!$                            (enkf_grid% dp0(i,j,n-1,d)+enkf_grid% dp0(i,j,n,d))))/enkf_mean% ps(i,j,1,d)
!!$                   else if(n==nend) then
!!$                      d2(i,j,n,d) = (enkf_mean% ps(i,j,1,d) - ((enkf_grid% dp0(i,j,n-1,d)*enkf_mean% pf(i,j,n,d)+&
!!$                                     enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n-1,d))/&
!!$                                     (enkf_grid% dp0(i,j,n-1,d)+enkf_grid% dp0(i,j,n,d)))) &
!!$                                     /enkf_mean% ps(i,j,1,d)
!!$                   else !(n==1)
!!$                      !d2(i,j,n,d) = d2(i,j,n+1,d)
!!$                       d2(i,j,n,d) = (((enkf_grid% dp0(i,j,n+1,d)*enkf_mean% pf(i,j,n,d)+&
!!$                                     enkf_grid% dp0(i,j,n,d)*enkf_mean% pf(i,j,n+1,d))/&
!!$                                     (enkf_grid% dp0(i,j,n+1,d)+enkf_grid% dp0(i,j,n,d)))-pht(i,j,1,1))/&
!!$                                     enkf_mean% ps(i,j,1,d)
!!$                       !(ph(k+1/2)-pht)/ps
!!$                   end if
!!$                end select
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$    !*************
!!$
!!$    !----------------------------------------------------------
!!$    ! Calculate rms difference between LETKF and 3dVar analysis
!!$    !----------------------------------------------------------
!!$
!!$    rms_ges    = 0._wp
!!$    rms_ges2   = 0._wp
!!$    spread_ges = 0._wp
!!$    bias_ges   = 0._wp
!!$    bias_ges2  = 0._wp
!!$    expl_var   = 0._wp
!!$    e_dim      = 0._wp
!!$    ml         = 0._wp
!!$    rms_e      = 0._wp
!!$    equ        = 0._wp
!!$    equ_2      = 0._wp
!!$    ew         = 0._wp
!!$    v_xb       = 0._wp
!!$    m_xb       = 0._wp
!!$    v_xb_nps   = 0._wp
!!$    m_xb_nps   = 0._wp
!!$    p_t        = 0._wp
!!$    u          = 0._wp
!!$    e          = 0._wp
!!$    u_qt       = 0._wp
!!$
!!$
!!$    select case (model)
!!$    case ('COSMO')
!!$       mask(:,:,:) = .true.
!!$    case('GME')
!!$       do i=1,n_n
!!$          do j=1,n_n
!!$             if  ((i==1     .and. j==1)   .or. &
!!$                  (i==n_n   .and. j==n_n) .or. &
!!$                  (i==1     .and. j==2)   .or. &
!!$                  (i==2     .and. j==1)   .or. &
!!$                  (i==n_n-1 .and. j==n_n) .or. &
!!$                  (i==n_n   .and. j==n_n-1)) then
!!$                mask(i,j,:) = .false.
!!$             else
!!$                mask(i,j,:) = .true.
!!$             end if
!!$          end do
!!$       end do
!!$    end select
!!$
!!$    do d=enkf_grid% lb(4), enkf_grid% ub(4)
!!$       do j=enkf_grid% lb(2), enkf_grid% ub(2)
!!$          do i=enkf_grid% lb(1), enkf_grid% ub(1)
!!$             n_z = 0
!!$             do n=1, nend
!!$                v_xb       = 0._wp
!!$                m_xb       = 0._wp
!!$                v_xb_nps   = 0._wp
!!$                m_xb_nps   = 0._wp
!!$                p_t        = 0._wp
!!$                u          = 0._wp
!!$                e          = 0._wp
!!$                u_qt       = 0._wp
!!$                m2 = 0
!!$                m_v = 0
!!$                ps_flag = 0
!!$                n_c = 0
!!$                if (n==6 .or. n==12 .or. n==19 .or. n==35) then
!!$                   n_z = n_z + 1
!!$                end if
!!$                do m=1, m_max
!!$                   if (enkf_mean1% m(m)% i% alloc .and.     &
!!$                        analysis% m(m)% i% alloc .and.         &
!!$                        .not. enkf_mean1% m(m)% i% ref .and. &
!!$                        .not. analysis% m(m)% i% ref           ) then
!!$                      if (spread_flag == 1) then
!!$                      select case(fc_flag)
!!$                      case(1)
!!$                      case(0)
!!$                      select case (enkf_mean1% m(m)% i% name)
!!$                      case('u','v','t','ps')
!!$                         if ((n==6 .or. n==12 .or. n==19 .or. n==35) .and. &
!!$                              (j>(enkf_grid% lb(2)+1)      .and. &
!!$                              i>(enkf_grid% lb(1)+1)      .and. &
!!$                              j<(enkf_grid% ub(2)-1)      .and. &
!!$                              i<(enkf_grid% ub(1)-1))) then
!!$                            if (n==35) then
!!$                               select case (enkf_mean1% m(m)% i% name)
!!$                               case('ps')
!!$                                  m_v = m_v + 1
!!$                                  v_xb(((m_v-1)*n_grid_m)+1:m_v*((n_n*n_n)-6))                                         &
!!$                                    = pack((log(enkf_mean2%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,1,d))                 &
!!$                                    - log(analysis%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,1,d))),mask(:,:,1)) * sqrt(R * T)
!!$                                  do k=1, k_enkf
!!$                                     m_xb(((m_v-1)*n_grid_m)+1:m_v*((n_n*n_n)-6),k)                                         &
!!$                                       = pack((log(enkf_ens2(k)%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,1,d))                   &
!!$                                       - log(enkf_mean2%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,1,d))),mask(:,:,1)) * sqrt(R * T)
!!$                                  end do
!!$                                  ps_flag = 1
!!$                               case('u','v','t')
!!$                                  m_v = m_v + 1
!!$                                  if (ps_flag == 1) then
!!$                                     n_c = (n_v-1)*((n_n*n_n)-6)
!!$                                  end if
!!$                                  v_xb(((m_v-1)*n_grid_m)+1-n_c:(m_v*n_grid_m)-n_c)                             &
!!$                                       = pack((enkf_mean2%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d) &
!!$                                       - analysis%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d))          &
!!$                                       * sqrt(d2(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d)),mask)
!!$                                  do k=1, k_enkf
!!$                                     m_xb(((m_v-1)*n_grid_m)+1-n_c:(m_v*n_grid_m)-n_c,k)                         &
!!$                                       = pack((enkf_ens2(k)%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d) &
!!$                                       - enkf_mean2%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d))    &
!!$                                          * sqrt(d2(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d)),mask)
!!$                                  end do
!!$                                  if (enkf_mean1% m(m)% i% name == 't') then
!!$                                     v_xb(((m_v-1)*n_grid_m)+1-n_c:(m_v*n_grid_m)-n_c) = &
!!$                                          v_xb(((m_v-1)*n_grid_m)+1-n_c:(m_v*n_grid_m)-n_c) * sqrt(C_p/T)
!!$                                     m_xb(((m_v-1)*n_grid_m)+1-n_c:(m_v*n_grid_m)-n_c,:) = &
!!$                                          m_xb(((m_v-1)*n_grid_m)+1-n_c:(m_v*n_grid_m)-n_c,:) * sqrt(C_p/T)
!!$                                  end if
!!$                               end select
!!$                            else if (n==6 .or. n==12 .or. n==19) then
!!$                               select case (enkf_mean1% m(m)% i% name)
!!$                               case('u','v','t')
!!$                                  m_v = m_v + 1
!!$                                  v_xb_nps(((m_v-1)*n_grid_m)+1:(m_v*n_grid_m))                                 &
!!$                                    = pack((enkf_mean2%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d) &
!!$                                    - analysis%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d))          &
!!$                                       * sqrt(d2(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d)),mask)
!!$                                  do k=1, k_enkf
!!$                                     m_xb_nps(((m_v-1)*n_grid_m)+1:(m_v*n_grid_m),k)                             &
!!$                                       = pack((enkf_ens2(k)%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d) &
!!$                                       - enkf_mean2%  m(m)% ptr(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d))    &
!!$                                          * sqrt(d2(i-n_i:i+n_i,j-n_i:j+n_i,n-n_v_i:n+n_v_i,d)),mask)
!!$                                  end do
!!$                                  if (enkf_mean1% m(m)% i% name == 't') then
!!$                                     v_xb_nps(((m_v-1)*n_grid_m)+1:(m_v*n_grid_m)) = &
!!$                                          v_xb_nps(((m_v-1)*n_grid_m)+1:(m_v*n_grid_m)) * sqrt(C_p/T)
!!$                                     m_xb_nps(((m_v-1)*n_grid_m)+1:(m_v*n_grid_m),:) = &
!!$                                          m_xb_nps(((m_v-1)*n_grid_m)+1:(m_v*n_grid_m),:) * sqrt(C_p/T)
!!$                                  end if
!!$                               case default
!!$                               end select
!!$                            end if !inner n-level selection
!!$                         end if !outer n_level selection
!!$                      case default
!!$                      end select !variable selection
!!$                   end select !fc_flag selection
!!$                end if !spread-flag selection
!!$             endif !allocation check
!!$          enddo !m-loop
!!$          select case(fc_flag)
!!$          case(1)
!!$          case(0)
!!$                if (spread_flag == 1) then
!!$                if ((n==6 .or. n==12 .or. n==19 .or. n==35) .and. &
!!$                    (j>(enkf_grid% lb(2)+1)      .and. &
!!$                     i>(enkf_grid% lb(1)+1)      .and. &
!!$                     j<(enkf_grid% ub(2)-1)      .and. &
!!$                     i<(enkf_grid% ub(1)-1))) then
!!$                   !calculate eigen-values/vectors
!!$                   if (n==35) then
!!$                      call eigen((1._wp/(k_enkf-1))*matmul(transpose(m_xb),m_xb),e,u)
!!$                      m_xb=matmul(m_xb,u)!m_xb = Q
!!$                      do k=1, k_enkf
!!$                         m_xb(:,k) = m_xb(:,k)/sqrt(dot_product(m_xb(:,k),m_xb(:,k)))
!!$                      end do
!!$                      u_qt = matmul(v_xb,m_xb(:,2:k_enkf)) !(x_b-x_t)^t*Q hat dim k_enkf-1
!!$                      rms_e_add = dot_product(v_xb,v_xb)
!!$                      ev_add = dot_product(u_qt,u_qt)/dot_product(v_xb,v_xb)
!!$                      do k=1,k_enkf-1
!!$                         e_i(k) = 1._wp/e(k+1)
!!$                      enddo
!!$                      ml_add = log(product(e(2:k_enkf)))!+dot_product(u_qt,(e_i*u_qt))+(k_enkf-1)*log(2._wp*pi)
!!$                      equ_add = sqrt(e_i)*u_qt
!!$                      ew_add =  e(2:k_enkf)
!!$                   else if (n==6 .or. n==12 .or. n==19) then
!!$                      call eigen((1._wp/(k_enkf-1))*matmul(transpose(m_xb_nps),m_xb_nps),e,u)
!!$                      m_xb_nps=matmul(m_xb_nps,u)!m_xb = Q
!!$                      do k=1, k_enkf
!!$                         m_xb_nps(:,k) = m_xb_nps(:,k)/sqrt(dot_product(m_xb_nps(:,k),m_xb_nps(:,k)))
!!$                      end do
!!$                      u_qt = matmul(v_xb_nps,m_xb_nps(:,2:k_enkf)) !(x_b-x_t)^t*Q hat dim k_enkf-1
!!$                      rms_e_add = dot_product(v_xb_nps,v_xb_nps)
!!$                      ev_add = dot_product(u_qt,u_qt)/dot_product(v_xb_nps,v_xb_nps)
!!$                      do k=1,k_enkf-1
!!$                         e_i(k) = 1._wp/e(k+1)
!!$                      enddo
!!$                      ml_add = log(product(e(2:k_enkf)))!+dot_product(u_qt,(e_i*u_qt))+(k_enkf-1)*log(2._wp*pi)
!!$                      equ_add = sqrt(e_i)*u_qt
!!$                      ew_add =  e(2:k_enkf)
!!$                   end if
!!$                   do k=1,k_enkf
!!$                      if (e(k) < 0) then
!!$                         e(k) = 0
!!$                      end if
!!$                   end do
!!$                   edim_add = dot_product(u_qt,(e_i*u_qt)) !sum(sqrt(e(2:k_enkf)))*sum(sqrt(e(2:k_enkf)))/sum(e(2:k_enkf))
!!$                   index_lon = index(i+1,j,d,1)
!!$                   index_lat = index(i+1,j,d,2)
!!$                   expl_var(index_lon,index_lat,1,n_z) = &
!!$                        expl_var(index_lon,index_lat,1,n_z) + ev_add
!!$                   e_dim(index_lon,index_lat,1,n_z) = &
!!$                        e_dim(index_lon,index_lat,1,n_z) + edim_add
!!$                   ml(index_lon,index_lat,1,n_z) = &
!!$                        ml(index_lon,index_lat,1,n_z) + ml_add
!!$                   rms_e(index_lon,index_lat,1,n_z) = &
!!$                        rms_e(index_lon,index_lat,1,n_z) + rms_e_add
!!$                   if (n==35) then
!!$                      equ(index_lon,index_lat,1,:) = &
!!$                           equ(index_lon,index_lat,1,:) + equ_add
!!$                      ew(index_lon,index_lat,1,:) = &
!!$                           ew(index_lon,index_lat,1,:) + ew_add
!!$                      equ_2(index_lon,index_lat,1,:) = &
!!$                           equ_2(index_lon,index_lat,1,:) + equ_add*equ_add
!!$                   else if (n==19) then
!!$                      equ(index_lon,index_lat,3,:) = &
!!$                           equ(index_lon,index_lat,3,:) + equ_add
!!$                      ew(index_lon,index_lat,2,:) = &
!!$                           ew(index_lon,index_lat,2,:) + ew_add
!!$                      equ_2(index_lon,index_lat,2,:) = &
!!$                           equ_2(index_lon,index_lat,2,:) + equ_add*equ_add
!!$                   else if (n==12) then
!!$                      equ(index_lon,index_lat,5,:) = &
!!$                           equ(index_lon,index_lat,5,:) + equ_add
!!$                      ew(index_lon,index_lat,3,:) = &
!!$                           ew(index_lon,index_lat,3,:) + ew_add
!!$                      equ_2(index_lon,index_lat,3,:) = &
!!$                           equ_2(index_lon,index_lat,3,:) + equ_add*equ_add
!!$                   end if
!!$                end if !n-level selection
!!$             end if !spread-flag selection
!!$           end select !fc-flag selection
!!$        enddo !n-loop
!!$      enddo !i-loop
!!$   enddo !j-loop
!!$ enddo ! d-loop
!!$    do k=1,n_lon
!!$       do l=1,n_lat
!!$          select case(fc_flag)
!!$          case(1)
!!$          case(0)
!!$          if (spread_flag == 1) then
!!$          if (n_ges_e(k,l) .ne. 0) then
!!$             expl_var(k,l,:,:) = expl_var(k,l,:,:)/n_ges_e(k,l)
!!$             e_dim(k,l,:,:) = e_dim(k,l,:,:)/n_ges_e(k,l)
!!$             ml(k,l,:,:) = ml(k,l,:,:)/n_ges_e(k,l)
!!$             rms_e(k,l,:,:) = rms_e(k,l,:,:)/n_ges_e(k,l)
!!$             equ(k,l,1,:) = equ(k,l,1,:)/n_ges_e(k,l)
!!$             equ(k,l,3,:) = equ(k,l,3,:)/n_ges_e(k,l)
!!$             equ(k,l,5,:) = equ(k,l,5,:)/n_ges_e(k,l)
!!$             ew(k,l,:,:) = ew(k,l,:,:)/n_ges_e(k,l)
!!$             equ(k,l,2,:) = (equ_2(k,l,1,:)/n_ges_e(k,l))-(equ(k,l,1,:)*equ(k,l,1,:))
!!$             equ(k,l,4,:) = (equ_2(k,l,2,:)/n_ges_e(k,l))-(equ(k,l,3,:)*equ(k,l,3,:))
!!$             equ(k,l,6,:) = (equ_2(k,l,3,:)/n_ges_e(k,l))-(equ(k,l,5,:)*equ(k,l,5,:))
!!$             if (sum(equ(k,l,2,:))>1.0e010) then
!!$                print *, 'k,l = ', k,l, equ(k,l,2,:), equ_2(k,l,1,:)/n_ges_e(k,l),equ(k,l,1,:)*equ(k,l,1,:)
!!$             end if
!!$          end if
!!$          end if
!!$         end select
!!$       enddo
!!$    enddo
!!$    if (spread_flag == 1) then
!!$       ml = p_sum(ml)
!!$       rms_e = p_sum(rms_e)
!!$       equ = p_sum(equ)
!!$       ew  = p_sum(ew)
!!$       expl_var = p_sum(expl_var)
!!$       e_dim = p_sum(e_dim)
!!$    end if
!!$  !-------------------------------
!!$  !write diagnostic fields (grads)
!!$  !-------------------------------
!!$  if (dace% lpio) then
!!$  ! write exp var, ml etc.
!!$     select case(fc_flag)
!!$     case(1)
!!$     case(0)
!!$     if (spread_flag == 1) then
!!$  ! write expl. var. to grads file
!!$     filename = path_file (output,'^exp_var__ANA_TIME_')
!!$     call init_ctl (ctl, file=filename, nx=72, ny=36, ke=4,&
!!$                    lo1=-177.5_wp, di=5._wp,la1=-87.5_wp, dj=5._wp)
!!$     call add_var (ctl,'ev',4,'ev')
!!$     call write_var (ctl,expl_var(:,:,1,:),'ev')
!!$     call write_ctl (ctl)
!!$     call destruct(ctl)
!!$
!!$  ! write ens. dim. to grads file
!!$     filename = path_file (output,'^ens_dim__ANA_TIME_')
!!$     call init_ctl (ctl, file=filename, nx=72, ny=36, ke=4,&
!!$                    lo1=-177.5_wp, di=5._wp,la1=-87.5_wp, dj=5._wp)
!!$     call add_var (ctl,'ed',4,'ed')
!!$     call write_var (ctl,e_dim(:,:,1,:),'ed')
!!$     call write_ctl (ctl)
!!$     call destruct(ctl)
!!$
!!$  ! write max. likelihood to grads file
!!$     filename = path_file (output,'^ml__ANA_TIME_')
!!$     call init_ctl (ctl, file=filename, nx=72, ny=36, ke=4,&
!!$                    lo1=-177.5_wp, di=5._wp,la1=-87.5_wp, dj=5._wp)
!!$     call add_var (ctl,'ml',4,'ml')
!!$     call write_var (ctl,ml(:,:,1,:),'ml')
!!$     call write_ctl (ctl)
!!$     call destruct(ctl)
!!$
!!$  ! write equ to grads file
!!$     filename = path_file (output,'^equ__ANA_TIME_')
!!$     call init_ctl (ctl, file=filename, nx=72, ny=36, ke=k_enkf-1,&
!!$                    lo1=-177.5_wp, di=5._wp,la1=-87.5_wp, dj=5._wp)
!!$     call add_var (ctl,'equ',k_enkf-1,'equ')
!!$     call add_var (ctl,'equ_s',k_enkf-1,'equ_s')
!!$     call add_var (ctl,'ew',k_enkf-1,'ew')
!!$     call add_var (ctl,'equ500',k_enkf-1,'equ500')
!!$     call add_var (ctl,'equ_s500',k_enkf-1,'equ_s500')
!!$     call add_var (ctl,'ew500',k_enkf-1,'ew500')
!!$     call add_var (ctl,'equ250',k_enkf-1,'equ250')
!!$     call add_var (ctl,'equ_s250',k_enkf-1,'equ_s250')
!!$     call add_var (ctl,'ew250',k_enkf-1,'ew250')
!!$     call write_var (ctl,equ(:,:,1,:),'equ')
!!$     call write_var (ctl,equ(:,:,2,:),'equ_s')
!!$     call write_var (ctl,ew(:,:,1,:),'ew')
!!$     call write_var (ctl,equ(:,:,3,:),'equ500')
!!$     call write_var (ctl,equ(:,:,4,:),'equ_s500')
!!$     call write_var (ctl,ew(:,:,2,:),'ew500')
!!$     call write_var (ctl,equ(:,:,5,:),'equ250')
!!$     call write_var (ctl,equ(:,:,6,:),'equ_s250')
!!$     call write_var (ctl,ew(:,:,3,:),'ew250')
!!$     call write_ctl (ctl)
!!$     call destruct(ctl)
!!$  end if
!!$  end select
!!$
!!$endif
!!$
!!$ end subroutine ens_statistics

!=============================================================================
!SR read_analysis

  subroutine read_analysis  (analysis,grid,ana_time,par_stat,fname)

  type (t_atm)      ,intent(inout) :: analysis ! deterministic analysis
  type (t_grid)     ,pointer       :: grid     ! atmospheric grid  to read
  type (t_time)     ,intent(in)    :: ana_time ! valid time of atmos. state
  character(len=*)  ,intent(in)    :: par_stat ! parameters to read
  character(len=*)                 :: fname    ! file name
  type (t_inventory),pointer       :: invt(:) => NULL()
  integer                          :: iunit    ! unit number (file pres. check)
  integer                          :: ios      ! iostat return variable

  !--------------------
  ! read analysis
  !--------------------

   !--------------------------------------------
   ! derive file name, check if file is readable
   !--------------------------------------------
    if (dace% lpio) then
      iunit = get_unit_number()
      open (iunit, file=fname, status='old', action='read', iostat=ios)
      if (ios/=0) call finish ('read_analysis','cannot open '//trim(fname))
      close (iunit)
      call return_unit_number (iunit)
      if (dace% lpio) then
        write(6,'(a,i4)') '    reading file: '//trim(fname)
      endif
    end if
    call get_inventory   (invt, fname)
    !----------------
    ! read 3dVar file
    !----------------
    if (dace% lpio) then
       write (6,'(a)')
       write (6,'(a)') 'reading analysis'
       write (6,'(a)')
       call print_inventory (invt, first=.true.)
       write (6,'(a)')
       write (6,'(a,a)') 'par_stat = ', par_stat
       write (6,'(a)')
    endif

!   print *, 'par_stat = ', par_stat

    call read (analysis, fname, invt, &
         time    = ana_time,          &
         fields  = trim(par_stat))

  end subroutine read_analysis
!==============================================================================
  subroutine compute_mean (enkf_mean, enkf_spread, enkf_ens, enkf_grid,       &
                           ana_time, k_enkf, iopath, ens_fname, write_files,  &
                           src_fc, par_mean, det_fname, verbose)
  !------------------------------------------------
  ! determine mean and spread of ensemble members
  ! optional: adjust ensemble around (deterministic)
  ! field
  !------------------------------------------------
    type(t_atm)      ,intent(inout)       :: enkf_mean  ! mean to compute
    type(t_atm)      ,intent(inout)       :: enkf_spread! spread to compute
    type(t_atm)      ,intent(inout)       :: enkf_ens(:)! ens to process
    type(t_grid)     ,pointer             :: enkf_grid  ! atmosph. grid to read
    type(t_time)     ,intent(in)          :: ana_time   ! veri time of fields
    integer          ,intent(in)          :: k_enkf     ! number ens members
    character(len=*) ,intent(in)          :: iopath     ! path (input/output)
    character(len=*) ,intent(in)          :: ens_fname  ! ens file name
    logical          ,intent(in)          :: write_files! write mean / spread
    integer          ,intent(in)          :: src_fc     ! 1: adjust ensemble
                                                        ! 0: do nothing
    character(len=*) ,intent(in),optional :: par_mean   ! param. in mean/spread
    character(len=*) ,intent(in),optional :: det_fname  ! (det) file name to
                                                        ! adjust ens around
    logical          ,intent(in),optional :: verbose    ! print mean / spread?

    !---------------
    !local variables
    !---------------
    integer                     :: k       ! ensemble member index
    integer                     :: m       ! parameter index
    character(len=128)          :: fname   ! full file name
    character(len=256)          :: par_m   ! copy of par_mean
    logical                     :: verb    ! copy of verbose
    type (t_atm)                :: fc      ! deterministic field
    type (t_inventory), pointer :: invt(:) => NULL()

    verb = .true.; if (present (verbose)) verb = verbose
    !--------------------------
    ! determine mean and spread
    !--------------------------
    call stop_time ('compute_mean: determine mean and spread')
    if (present (par_mean)) then
      !--------------------------------------------------
      ! remove certain fields from list
      ! aerosol constituents: no GRIB2 template defined ?
      !--------------------------------------------------
      par_m = par_mean
      call eval_string (par_m, '- aer_so4 aer_dust aer_org aer_bc aer_ss')
      call construct (enkf_mean,   template=enkf_ens(1), alloc=par_m)
      call construct (enkf_spread, template=enkf_ens(1), alloc=par_m)
    else
      call construct (enkf_mean,   template=enkf_ens(1))
      call construct (enkf_spread, template=enkf_ens(1))
    endif
    do m = 1, size (enkf_ens(1)% m)
      if (.not. enkf_ens(1)% m(m)% i% alloc) then
        call deallocate (enkf_mean,   enkf_ens(1)% m(m)% i% name)
        call deallocate (enkf_spread, enkf_ens(1)% m(m)% i% name)
      endif
    end do
    !--------------
    ! set meta data
    !--------------
    enkf_mean  % runtype     = enkf_ens(1)% runtype
    enkf_mean  % runclass    = enkf_ens(1)% runclass
    enkf_mean  % expid       = enkf_ens(1)% expid
!   enkf_mean  % ensemble_id = proc_ana_ens       ! or enkf_ens(1)% ensemble_id
    enkf_mean  % member      = ENS_MEAN
    enkf_mean  % members     = k_enkf
    enkf_spread% runtype     = enkf_ens(1)% runtype
    enkf_spread% runclass    = enkf_ens(1)% runclass
    enkf_spread% expid       = enkf_ens(1)% expid
!   enkf_spread% ensemble_id = enkf_mean  % ensemble_id
    enkf_spread% member      = ENS_SPREAD
    enkf_spread% members     = k_enkf

!   enkf_mean = enkf_ens(1)
!   do k = 2, k_enkf
!     enkf_mean = enkf_ens(k) + enkf_mean
!   end do
!   enkf_mean = enkf_mean * (1._wp/k_enkf)
!
!   enkf_spread = 0._wp
!   do k = 1, k_enkf
!     enkf_spread =  (enkf_ens(k) - enkf_mean)**2 + enkf_spread
!   end do
!   enkf_spread = enkf_spread*(1._wp/(k_enkf-1))
!   enkf_spread = sqrt(enkf_spread)

    call mean_stdev (enkf_mean, enkf_spread, enkf_ens)

    !---------
    ! printout
    !---------
    if (verb) then
      call stop_time ('compute_mean: printout')
      call print (enkf_mean,   verbose=.true.,              &
                  comment='ensemble mean: '//trim(ens_fname))
      call print (enkf_spread, verbose=.true., grid=.false.,  &
                  comment='ensemble spread: '//trim(ens_fname))
    end if

    !---------------------------
    ! adjust ensemble (optional)
    !---------------------------
    select case (src_fc)
    case (0)
    case (1)
      call stop_time ('compute_mean: adjust ensemble')
      !---------------
      ! read inventory
      !---------------
      if(dace% lpio) then
        write (6,'(a)') repeat('-',79)
        write (6,'(a)')
        write (6,'(a)') 'inventory of grib file: '//trim(det_fname)
        write (6,'(a)')
      endif
      call get_inventory (invt, det_fname)
      call print_inventory (invt, first=.true.)
      !----------------------------
      ! read deterministic forecast
      !----------------------------
      if (dace% lpio) then
        write (6,'(a)')
        write (6,'(a)') 'reading background from grib file: '//trim(det_fname)
        write (6,'(a)')
      endif
      call construct (fc ,enkf_grid, alloc=trim(par_mean))
      call read (fc, det_fname, invt,          &
                 !runtype = runtype,          &
                 time    = ana_time,         &
                 fields  = trim(par_mean)    )
      call print (fc,verbose=.true.,comment='deterministic field')
      !----------------------------------------------
      ! restrict specific humidity to positive values
      ! (appeared in ECMWF analyses)
      !----------------------------------------------
      fc% q = max(fc% q, 0._wp)
      !-------------------------
      ! determine offset, adjust
      !-------------------------
      fc = fc - enkf_mean
      call print(fc,verbose=.true.,comment='det field minus ensemble mean')
      enkf_mean = enkf_mean + fc
      call print(enkf_mean,verbose=.false.,&
                   comment='adjusted mean')
      do k = 1, k_enkf
        enkf_ens(k) = enkf_ens(k) + fc
        if (associated(enkf_ens(k)% psr)) enkf_ens(k)% psr = enkf_ens(k)% ps
        call print(enkf_ens(k),verbose=.false.,&
                   comment='adjusted ensemble member')
      end do
      !--------
      ! cleanup
      !--------
      call destruct (fc)
      deallocate  (invt)
    case default
      call finish ('compute_mean','invalid value of src_fc')
    end select

    !+++++++++++++++++++++++++++++++++++++++++++++++
    ! check: removal changes deterministic GRIB-file
    !+++++++++++++++++++++++++++++++++++++++++++++++
    call set_defaults    (center, subcenter, proc_ana_ens)
    call set_dwd_def_ens (k_enkf, -1, enkf_ens(1)% ensemble_id)

    !-----------------
    ! write GRIB files
    !-----------------
    if (write_files) then
      !-----
      ! mean
      !-----
      call stop_time ('compute_mean: write mean (GRIB)')
      fname = path_file (iopath, basename (ens_fname, '.m_ENS_'), '.mean')
!     call deallocate (enkf_mean,   'psr')
!     call deallocate (enkf_mean,   'pf')
      call deallocate (enkf_mean, dealloc_fields)
      call write_grib (enkf_mean, file=fname ,mode='w', &
                       ref=.false.                      )
      !-------
      ! spread
      !-------
      call stop_time ('compute_mean: write spread (GRIB)')
      fname = path_file (iopath, basename (ens_fname, '.m_ENS_'), '.spread')
!     call deallocate (enkf_spread, 'psr')
!     call deallocate (enkf_spread, 'pf')
      call deallocate (enkf_spread, dealloc_fields)
      call write_grib (enkf_spread, file=fname ,mode='w', &
                       ref=.false.                        )
    end if

  end subroutine compute_mean

!==============================================================================
  subroutine subtract_mean (enkf_ens,enkf_mean,par_trans,k_enkf)
  !---------------------------------------------
  ! subtract mean from forecast ensemble members
  ! (only variables in par_trans!)
  !---------------------------------------------
    type(t_atm)        ,pointer        :: enkf_ens(:)! ens to process
    type(t_atm)        ,intent(in)     :: enkf_mean  ! spread to compute
    character(len=*)   ,intent(in)     :: par_trans  ! parameters
    integer            ,intent(in)     :: k_enkf     ! number of ens members

    !---------------
    !local variables
    !---------------
    integer                     :: k        ! ensemble member index
    character(len=16)           :: pars(nm) ! names of variables to derive
    integer                     :: n        ! return code
    integer                     :: m        ! variable number

    call split (pars, par_trans, n)
    do m=1, size(enkf_ens(1)% m)
       if (      enkf_ens(1)% m(m)% i% alloc .and. &
           .not. enkf_ens(1)% m(m)% i% ref         ) then
          if(any(pars==enkf_ens(1)% m(m)% i% name)) then
             if (dace% lpio)                                          &
               write(6,*) 'subtract_mean: ', enkf_ens(1)% m(m)% i% name
             if (.not. enkf_mean% m(m)% i% alloc)                        &
               call finish('subtract_mean',                              &
                           'mean not present for '//enkf_mean%m(m)%i%name)
!$omp parallel do
             do k = 1, k_enkf
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
                enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
                enkf_ens(k)% m(m)% ptr(:,:,:,:) - &
                enkf_mean  % m(m)% ptr(:,:,:,:)
             end do
!$omp end parallel do
          end if
       end if
    end do
    do k = 1, k_enkf
       call print (enkf_ens(k), grid=(k==1),                       &
                   comment='forecast ensemble deviation '//char3(k))
    end do

  end subroutine subtract_mean
!==============================================================================
  subroutine add_mean (enkf_ens, enkf_mean, par_trans, k_enkf)
  !--------------------------------------
  ! add mean to ensemble members
  ! (only variables listed in par_trans!)
  !--------------------------------------
    type(t_atm)        ,pointer        :: enkf_ens(:)! ens to process
    type(t_atm)        ,intent(in)     :: enkf_mean  ! mean to compute
    character(len=*)   ,intent(in)     :: par_trans  ! parameters
    integer            ,intent(in)     :: k_enkf     ! number of ens members

    !---------------
    !local variables
    !---------------
    integer                     :: k        ! ensemble member index
    character(len=16)           :: pars(nm) ! names of variables to derive
    integer                     :: n        ! return code
    integer                     :: m        ! variable number

    call split (pars, par_trans, n)
    do m=1, size(enkf_ens(1)% m)
      if(any(pars==enkf_ens(1)% m(m)% i% name)) then
        if (dace% lpio) &
          write(6,*) ' add_mean: ', enkf_ens(1)% m(m)% i% name,  &
                                    enkf_ens(1)% m(m)% i% alloc, &
                              .not. enkf_ens(1)% m(m)% i% ref,   &
                                    enkf_mean  % m(m)% i% alloc
        if (      enkf_ens(1)% m(m)% i% alloc .and. &
            .not. enkf_ens(1)% m(m)% i% ref         ) then
          if (.not. enkf_mean% m(m)% i% alloc)                        &
            call finish('add_mean',                                   &
                        'mean not present for '//enkf_mean%m(m)%i%name)
!$omp parallel do
          do k = 1, k_enkf
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
            enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
            enkf_ens(k)% m(m)% ptr(:,:,:,:) + &
            enkf_mean  % m(m)% ptr(:,:,:,:)
          end do
!$omp end parallel do
        end if
      end if
    end do

  end subroutine add_mean
!==============================================================================
  subroutine add_noise (enkf_ens, enkf_additive, par_read_additive, k_enkf, scale_factor)
  !----------------------------------------------------
  ! add climatological perturbation to ensemble members
  ! (only variables listed in par_read_additive!)
  !----------------------------------------------------
    type(t_atm)        ,pointer        :: enkf_ens(:)       ! ens to process
    type(t_atm)        ,pointer        :: enkf_additive(:)  ! climatological perturbation
    character(len=*)   ,intent(in)     :: par_read_additive ! parameters
    integer            ,intent(in)     :: k_enkf            ! number of ens members
    real(wp)           ,intent(in)     :: scale_factor      ! scale factor

    !---------------
    !local variables
    !---------------
    integer                     :: k        ! ensemble member index
    character(len=16)           :: pars(nm) ! names of variables to derive
    integer                     :: n        ! return code
    integer                     :: m        ! variable number

    call split (pars, par_read_additive, n)
    do m=1, size(enkf_ens(1)% m)
      if(any(pars==enkf_ens(1)% m(m)% i% name)) then
        if (dace% lpio) &
          write(6,*) ' add_perturbation: ', enkf_ens(1)% m(m)% i% name,  &
                                            enkf_ens(1)% m(m)% i% alloc, &
                                            .not. enkf_ens(1)% m(m)% i% ref,   &
                                            enkf_additive(1)% m(m)% i% alloc
        if (      enkf_ens(1)% m(m)% i% alloc .and. &
            .not. enkf_ens(1)% m(m)% i% ref         ) then
          if (.not.enkf_additive(1)% m(m)% i% alloc)                         &
            call finish('add_noise',                                         &
                        'perturbation not present for '//enkf_additive(1)%m(m)%i%name)
!$omp parallel do
          do k = 1, k_enkf
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
            enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
            enkf_ens(k)% m(m)% ptr(:,:,:,:) + &
            scale_factor*enkf_additive(k)% m(m)% ptr(:,:,:,:)
          end do
!$omp end parallel do
        end if
      end if
    end do

  end subroutine add_noise
!==============================================================================
   subroutine apply_H (H_x, obs, enkf_ens, k_enkf, H_x_mean, enkf_mean, &
                       check, bg, bge)
    !------------------------------------------------------------------------
    ! apply H operator on ensemble [and optionally ensemble mean]
    ! the ensemble enkf_ens(:) is mapped to enkf_ens(:,1) for dummy time slot
    !------------------------------------------------------------------------
    type(t_vector) ,pointer                 :: H_x (:)     ! H applied
    type(t_obs_set),intent(inout)           :: obs         ! observation info
    type(t_atm)    ,intent(inout)           :: enkf_ens(:) ! ensemble
    integer        ,intent(in)              :: k_enkf      ! ensemble size
    type(t_vector) ,intent(inout) ,optional :: H_x_mean    ! H applied to mean
    type(t_atm)    ,intent(inout) ,optional :: enkf_mean   ! ensemble mean
    integer        ,intent(in)    ,optional :: check       ! check operator
    logical        ,intent(in)    ,optional :: bg          ! called   for bg/fg
    logical        ,intent(in)    ,optional :: bge         ! ensemble for bg/fg
    target                                  :: enkf_ens

      type(t_atm) ,pointer :: enkf_ens_t  (:,:) ! ensemble + time slot
      type(t_atm)          :: enkf_mean_t (1)   ! ensemble + time slot

      enkf_ens_t (1:k_enkf, 1:1) => enkf_ens (1:k_enkf)

      if (present (enkf_mean)) then
        call flat_copy (enkf_mean_t(1), enkf_mean)
        call apply_H2 (H_x, obs, enkf_ens_t, k_enkf, H_x_mean, enkf_mean_t, &
             check=check, bg=bg, bge=bge)
        call flat_copy (enkf_mean, enkf_mean_t(1))
      else
        call apply_H2 (H_x, obs, enkf_ens_t, k_enkf, check=check,bg=bg,bge=bge)
      endif

  end subroutine apply_H
!==============================================================================
  subroutine apply_H2 (H_x, obs, enkf_ens, k_enkf, H_x_mean, enkf_mean, &
                       check, bg, bge)
    !------------------------------------------------------------
    ! apply H operator on ensemble [and optionally ensemble mean]
    ! enkf_ens(:,:) holds the ensemble        (dim=1)
    !               with different time slots (dim=2)
    ! H_x     (:)   holds the ensemble
    ! enkf_mean (:) holds time slots
    !------------------------------------------------------------
    type(t_vector) ,pointer                :: H_x (:)       ! H applied
    type(t_obs_set),intent(inout)          :: obs           ! observation info
    type(t_atm)    ,intent(inout)          :: enkf_ens(:,:) ! ensemble
    integer        ,intent(in)             :: k_enkf        ! ensemble size
    type(t_vector) ,intent(inout),optional :: H_x_mean      ! H applied to mean
    type(t_atm)    ,intent(inout),optional :: enkf_mean(:)  ! ensemble mean
    integer        ,intent(in)   ,optional :: check         ! check operator
    logical        ,intent(in)   ,optional :: bg            ! called for bg/fg
    logical        ,intent(in)   ,optional :: bge         ! ensemble for bg/fg

    !---------------
    !local variables
    !---------------
    integer        :: k     ! ensemble member index
    type(t_vector) :: Hi_x  ! state interp.to obs.locations
    integer        :: c     ! check for applicability of operator
    integer        :: n     ! number of operators to be removed
    logical        :: full  ! force TSK_SETUP_FULL on mean

    call stop_time ('apply H operator')
    if (dace% lpio) then
      write(6,'(a)') '    set up interpolation operator'
      write(6,'()')
    endif
    c = 0; if (present (check)) c = check
    n = 0
    !----------------
    ! allocate result
    !----------------
    allocate (H_x (k_enkf))
    !-------------------------------------
    ! determine model columns required
    ! for preconditioning boxes on each PE
    !-------------------------------------
    call process_obs (TSK_SETUP_COLS, obs, enkf_ens(1,1))
    call test_obs    (obs% o,'TSK_SETUP_COLS',0)
    if (dace% lpio) then
      write(6,'(a)') '    Apply H operator'
      write(6,'()')
    endif
    !----------------------------------------------------------------
    ! Need TSK_SETUP_FULL once for different model grid or background
    !----------------------------------------------------------------
    full = .false.
    if (present (bg)) full = bg
    !-----------------------
    ! apply to ensemble mean
    !-----------------------
    if (present (enkf_mean)) then
      call construct (H_x_mean, obs% oi)
      call apply_H_m (H_x_mean, obs, enkf_mean, Hi_x, TSK_K, bg=bg, setup=full)
    endif
    !---------------------------
    ! loop over ensemble members
    !---------------------------
    do k=1, k_enkf
      call apply_H_m (H_x(k), obs, enkf_ens(k,:), Hi_x, TSK_Y, bg=bge)
      call check_operator (obs, c, n, H_atm=H_x(k))
      if (dace% lpio .and. n>0) write(6,*) n,'observations rejected'
    enddo
    call destruct (Hi_x)

  end subroutine apply_H2
!============================================================================

  subroutine apply_H_m (H_x, obs, atm, Hi_x, task, bg, setup)
  !------------------------------------------------------------
  ! apply H to single ensemble member using multiple time slots
  !------------------------------------------------------------
  type (t_atm)   ,intent(inout)          :: atm(:) ! member to apply H on
  type(t_obs_set),intent(inout)          :: obs    ! observation derived type
  type(t_vector) ,intent(inout)          :: H_x    ! H applied to atm
  type(t_vector) ,intent(inout),optional :: Hi_x   ! temporary
  integer        ,intent(in)   ,optional :: task   ! TSK_Y or TSK_K
  logical        ,intent(in)   ,optional :: bg     ! called for bg/fg
  logical        ,intent(in)   ,optional :: setup  ! force TSK_SETUP_FULL?
  target :: Hi_x
    !----------------
    !local variables)
    !----------------
    integer                :: t                           ! time slot
    type(t_cols),pointer   :: cbgb(:)                     ! background (columns)
    logical                :: apf,aph,agh,agf,arh         ! alloc. status to keep
    logical                :: atv,aps,ar2,arl             ! alloc. status to keep
    integer                :: ntri                        ! number of boxes
    integer                :: tsk                         ! forward or Jakobian
    type(t_vector),pointer :: tmp
    logical                :: full                        ! need TSK_SETUP_FULL

      if (dace% lpio) write(6,'(a,i4)') '    member:', atm(1)% member
      tsk = TSK_Y; if (present (task)) tsk = task
      if (present (Hi_x)) then
        tmp => Hi_x
      else
        allocate (tmp)
      endif
      full = .false.; if (present (setup)) full = setup
      !--------------------------------------------------------
      ! derive all quantities required by observation operators
      !--------------------------------------------------------
      atv = associated (atm(1)% tv)
      aps = associated (atm(1)% ps)
      apf = associated (atm(1)% pf)
      aph = associated (atm(1)% ph)
      agh = associated (atm(1)% geoh)
      agf = associated (atm(1)% geof)
      arh = associated (atm(1)% rh)
      ar2 = associated (atm(1)% rh2m)
      arl = associated (atm(1)% rh2m_land)
      do t=1,size(atm)
        call set_tv               (atm(t))
        if (.not.apf) call set_p  (atm(t))
        call set_geo              (atm(t), geof=.true.)
        call allocate             (atm(t),'rh')
        atm(t)% rh = rh_q (atm(t)% q, atm(t)% t, atm(t)% pf)
        if (.not. associated (atm(t)% rh2m) &
            .and. associated (atm(t)% t2m ) &
            .and. associated (atm(t)% td2m) ) call set_rh2m(atm(t)) ! 2m rel.hum.
        if (.not. associated (atm(t)% rh2m_land) &
            .and. associated (atm(t)% t2m_land ) &
            .and. associated (atm(t)% td2m_land) ) call set_rh2m_land (atm(t),rhmax=1._wp)
        if (.not.aps) call set_ps (atm(t))
        if (atm(t)% grid% model == MO_ICON) &
          call set_ph       (atm(t))
      end do
      !----------------------------------
      ! send model columns to hosting PEs
      !----------------------------------
      ntri = size (obs% o)
      allocate (cbgb (ntri))
      call get_cols (obs% o% mc, atm, cbgb)
      !------------------------
      ! store ln p instead of p
      !------------------------
      call set_cols_logp (cbgb)
      !------------------------------
      ! allocate temporary and result
      !------------------------------
      if (.not.associated (obs% ii) .or. full) then
        call process_obs (TSK_SETUP_FULL, obs, atm(1), cbgb, local=.true.)
        call test_obs    (obs% o,'TSK_SETUP_FULL',0)
      end if
      if (.not.associated (tmp% s)) &
        call construct (tmp, obs% ii,'Hi_x') ! temporary in interpolation space
      if (.not.associated (H_x% s)) &
        call construct (H_x, obs% oi,'H_x' ) ! obs.operator applied
      !---------------------------
      ! apply observation operator
      !---------------------------
      call interpolate (cbgb,  obs, psi = tmp, bg=bg)
      call interpolate_strat(obs)
      call destruct (obs% l% H)
      call process_obs (tsk, obs, atm(1), cbgb, xi=tmp, y=H_x, context=POC_FG)
      !---------
      ! clean up
      !---------
      call dealloc_cols (cbgb)
      deallocate        (cbgb)
      if (.not.present (Hi_x)) then
        call destruct (tmp)
        deallocate    (tmp)
      endif
      if (.not.aph) call deallocate   (atm,'ph')
      if (.not.apf) call deallocate   (atm,'pf')
      if (.not.aps) call deallocate   (atm,'ps')
      if (.not.arh) call deallocate   (atm,'rh')
      if (.not.atv) call deallocate   (atm,'tv')
      if (.not.agh) call deallocate   (atm,'geoh')
      if (.not.agf) call deallocate   (atm,'geof')
      if (.not.ar2) call deallocate   (atm,'rh2m')
      if (.not.arl) call deallocate   (atm,'rh2m_land')

    end subroutine apply_H_m

!==============================================================================
subroutine add_model_error_new(ens,enkf_grid,cbgb,a_m,k_enkf,m_flag,mf,mf0,moderr_fc)

   type(t_atm) ,intent(inout) :: ens (:)  ! ensemble
   type (t_grid),pointer      :: enkf_grid! atmosph. grid
   integer     ,intent(in)    :: k_enkf   ! number of ens members
   integer     ,intent(in)    :: moderr_fc! 1: apply mod. err. to fg, 0: to ana
   integer     ,intent(in)    :: m_flag   ! model error mode
   real(wp)    ,intent(in)    :: mf       ! B-factor
   real(wp)    ,intent(in)    :: mf0      ! B-factor (mean)
   type (t_cols) ,pointer     :: cbgb(:)  ! background, redistributed
   type (t_cols) ,pointer     :: a_m(:)   ! member, redistributed

   !---------------
   !local variables
   !---------------

   real(wp) ,pointer :: x(:,:,:)  ! cntrol variables in wavelet representation
   real(wp) ,pointer :: x0(:,:,:) ! cntrol variables in wavelet representation (for mean)
   integer           :: k         ! ensemble member index
   integer           :: pespost=1 ! stride (parallel post-multp.)
   type(t_atm)       :: fci       ! increment (model error)
   type(t_atm)       :: fin       ! final ensemble member

    call stop_time ('apply add model error')

    !-------------------------
    ! allocate local variables
    !-------------------------
    if (mf < 0._wp) return
    allocate (x (covm%nx, covm%ny, covm%hc2_part_size))
    allocate (x0 (covm%nx, covm%ny, covm%hc2_part_size))

    !---------------------------------
    ! generate random numbers for mean
    !---------------------------------
!   print *, 'm_flag = ', m_flag, mf,mf0

    if (moderr_fc == 1 .and. m_flag == 0) then
       !-----------------------------------
       ! transform back to (rotated)/C-grid
       !-----------------------------------
       select case(enkf_grid% gridtype)
       case(WMO6_ROTLL, WMO6_LATLON)
          do k=1, k_enkf
             if (enkf_grid% rot) then
                call uvtrafo(ens(k),0)
             end if
             if (enkf_grid% arakawa == 'C') then
                call A_to_C_grid (ens(k)  ,restore=.true. )
             end if
          end do
       case default
       end select
       return
    else if (moderr_fc == 0) then
       select case (enkf_grid% gridtype)
       case (WMO6_LATLON, WMO6_ROTLL)
          do k=1, k_enkf
             if (enkf_grid% arakawa == 'C') then
                call C_to_A_grid (ens(k), save=.true.)
             endif
             if (enkf_grid% rot) then
                call uvtrafo(ens(k),1) !transform wind from rotated coord.
             endif
          end do
       case default
       end select
    end if

    select case (m_flag)
    case(1,2,3)
       select case(m_flag)
       case(1)
       case(2,3)
          call random_gauss (x0, covm)
       end select
       !---------------------------------------
       ! loop over ensemble members to generate
       !---------------------------------------
       do k = 1, k_enkf
          !--------------------------------------
          ! set up grid columns for communication
          !--------------------------------------
          if (k==1) then
             allocate (a_m (size(cbgb)) )
             call setup_post (cbgb, ens(k), pespost)
          end if
          !---------------------------------------------
          ! generate normally distributed random numbers
          !---------------------------------------------
          if (dace% lpio) then
             write(6,'()')
             write(6,'(a,i3)') '    forecast ensemble member ',k
             write(6,'()')
          endif
          if (m_flag==1 .or. m_flag==2) then
             call random_gauss (x, covm)
          endif
          select case(m_flag)
          case(1)
             x = x * mf
          case(2)
             x = x * mf + x0 * mf0
          case(3)
             x = x0 * mf0
          end select
          !-----------------
          ! apply L (sqrt B)
          !-----------------
          call apply_L_m   (x, a_m, cbgb, lnewpl=.true.)
          call finish_post (fci, a_m, ens(k))
          call dealloc_cols (a_m)
          if (.not.associated(ens(k)% psr)) &
               call allocate (ens(k), 'psr')
          ens(k)% psr = ens(k)% ps
          if (.not.associated(ens(k)% rh)) then
             call set_p (ens(k))
             call allocate (ens(k), 'rh')
             ens(k)% rh = rh_q (ens(k)% q, ens(k)% t, ens(k)% pf)
          endif
          call finish_ana  (fin, fci, ens(k), .true., .false., 4, 0, .false.)
          if (associated (ens(k)% tsurf)) then
            call allocate (fin, 'tsurf')
            fin% tsurf = ens(k)% tsurf
          endif
          ens(k) = fin
          call destruct (fci)
          call destruct (fin)
          call deallocate (ens(k), 'pf')
          call deallocate (ens(k), 'ph')
          call deallocate (ens(k), 'rh')
          call dealloc_cols (cbgb)
          deallocate        (cbgb)
       end do
       !---------
       ! clean up
       !---------
       deallocate        (a_m)
       deallocate (x)
    end select
 end subroutine add_model_error_new
!==========================================================================
subroutine mean_spread_talag (obs, o, H_x, H_x0, talag, spread, k_enkf)
!------------------------------------------------------------------
! calculate mean, spread, rank histogram index in observation space
!------------------------------------------------------------------
type(t_obs_set) ,intent(in)     :: obs      ! observation meta data
type(t_vector)  ,intent(in)     :: o        ! observed values
type(t_vector)  ,intent(inout)  :: talag    ! rank histogram index
type(t_vector)  ,intent(inout)  :: spread   ! spread at obs. location
type(t_vector)  ,pointer        :: H_x (:)  ! observation operator applied
type(t_vector)  ,intent(inout)  :: H_x0     ! mean of H_x
integer         ,intent(in)     :: k_enkf   ! number of ens members

  !---------------
  !Local variables
  !---------------
  integer           :: ib             ! observation box index
  integer           :: i              ! observation index
  integer           :: k              ! ensemble member index
  real(wp)          :: hx_loc(k_enkf) ! local value of H_x
  real(wp)          :: hx0            ! local value of mean
  real(wp)          :: uu, vv, dd     ! mean wind components / direction
  real(wp)          :: d1, d2         ! wind direction bounds
  real(wp)          :: ov             ! observed value
  !-----------------------
  ! loop over observations
  !-----------------------
  do ib = 1, talag% n_s
    if (talag% s(ib)% pe == dace% pe) then
      do i=1, talag% s(ib)% n
        !---------------------------------------------------
        ! store ensemble, observation in temporary variables
        !---------------------------------------------------
        do k = 1, k_enkf
          hx_loc(k) = H_x(k)% s(ib)% x(i)
        enddo
        ov = o% s(ib)% x(i)
        !------------------------------------
        ! special handling for invalid values
        !------------------------------------
        if (any (hx_loc == rvind) .or. ov == rvind) then
          spread % s(ib)% x(i) = rvind
          H_x0   % s(ib)% x(i) = rvind
          talag  % s(ib)% x(i) = rvind
        else
          !----------------------------------------------------------
          ! special handling (reordering) for ensemble wind direction
          !----------------------------------------------------------
          if (obs% o(ib)% varno(i) == VN_DD) then
            uu = 0._wp
            vv = 0._wp
            do k = 1, k_enkf
              dd = d2r * hx_loc(k)
              uu = uu - sin (dd)
              vv = vv - cos (dd)
            end do
            if (uu == 0._wp .and. vv == 0._wp) then
              dd = 0._wp
            else
              dd = r2d * atan2 (-uu ,-vv)
              if (dd < 0._wp) dd = dd + 360._wp
            endif
            d1 = dd - 180._wp
            d2 = dd + 180._wp
            do k = 1, k_enkf
              if (hx_loc(k) <  d1) hx_loc(k) = hx_loc(k) + 360._wp
              if (hx_loc(k) >= d2) hx_loc(k) = hx_loc(k) - 360._wp
            enddo
            if (ov <  d1) ov = ov + 360._wp
            if (ov >= d2) ov = ov - 360._wp
          endif
          !---------------------------------------------
          ! calculate mean, spread, rank histogram index
          !---------------------------------------------
          hx0 = (1._wp/k_enkf) * sum (hx_loc)
          if (k_enkf > 1) then
            spread % s(ib)% x(i) = sqrt (sum((hx_loc-hx0)**2)*(1._wp/(k_enkf-1)))
          else
            spread % s(ib)% x(i) = 0._wp
          end if
          if (obs% o(ib)% varno(i) == VN_DD) then
            if (hx0 <    0._wp) hx0 = hx0 + 360._wp
            if (hx0 >= 360._wp) hx0 = hx0 - 360._wp
          endif
          H_x0   % s(ib)% x(i) = hx0
          talag  % s(ib)% x(i) = count (hx_loc < ov)
        endif
      enddo
    endif
  enddo

end subroutine mean_spread_talag
!==============================================================================
  subroutine gather_coarse (y, x, grid)
  !-------------------------------------------------------------------------
  ! gather fields allocated on the coarse grid (holding LETKF gain matrices)
  ! in contrast to the routines in mo_atm_transp these routines
  ! can handle grid decompositions with overlap (halo)
  ! there is an additional 5th dimension for multiple fields
  !-------------------------------------------------------------------------
  real(wp)      ,pointer    :: y (:,:,:,:,:) ! (OUT) gathered    field
  real(wp)      ,pointer    :: x (:,:,:,:,:) ! (IN)  distributed field
  type (t_grid) ,intent(in) :: grid          ! 'coarse' grid meta data

    integer :: pe             ! processor index variable
    integer :: lb (4), ub (4)
    integer :: lbl(5), ubl(5)
    integer :: lbpio(4)       ! lower bounds for all procs
    integer :: ubpio(4)       ! upper bounds for all procs

    lbl = lbound(x)
    ubl = ubound(x)
    lb  = max (grid% lbg, lbl(1:4))
    ub  = min (grid% ubg, ubl(1:4))

    lbpio = lb
    ubpio = ub

    if (.not.dace% lpio) then
       call p_send (lb, dace% pio, 1)
       call p_send (ub, dace% pio, 1)
       if (lb(1) <= ub(1) .and. lb(2) <= ub(2))                     &
       call p_send (x (lb(1):ub(1), lb(2):ub(2),:,:,:), dace% pio, 1)
    else
      y = 0._wp
      do pe = 0, dace% npe -1
        if (pe /= dace% pio) then
          call p_recv (lbpio, pe, 1)
          call p_recv (ubpio, pe, 1)
          if (lbpio(1) <= ubpio(1) .and. lbpio(2) <= ubpio(2))             &
          call p_recv (y(lbpio(1):ubpio(1), lbpio(2):ubpio(2),:,:,:), pe, 1)
        else
          y (lb(1):ub(1), lb(2):ub(2),:,:,:) = x(lb(1):ub(1), lb(2):ub(2),:,:,:)
        endif
      end do
    endif
  end subroutine gather_coarse
!------------------------------------------------------------------------------
  subroutine scatter_coarse (y, x, grid)
  !--------------------------------------------------------------------------
  ! scatter fields allocated on the coarse grid (holding LETKF gain matrices)
  ! in contrast to the routines in mo_atm_transp these routines
  ! can handle grid decompositions with overlap (halo)
  ! there is an additional 5th dimension for multiple fields
  !--------------------------------------------------------------------------
  real(wp)      ,pointer    :: y (:,:,:,:,:) ! (IN)  gathered    field
  real(wp)      ,pointer    :: x (:,:,:,:,:) ! (OUT) distributed field
  type (t_grid) ,intent(in) :: grid          ! 'coarse' grid meta data

    integer :: pe             ! processor index variable
    integer :: lb (4), ub (4)
    integer :: lbl(5), ubl(5)
    integer :: lbpio(4)       ! lower bounds for all procs
    integer :: ubpio(4)       ! upper bounds for all procs

    lbl=lbound(x)
    ubl=ubound(x)
    lb = max (grid% lbg, lbl(1:4))
    ub = min (grid% ubg, ubl(1:4))

    if (dace% lpio) then
      do pe = 0, dace% npe -1
        if (pe /= dace% pio) then
          call p_recv (lbpio, pe, 1)
          call p_recv (ubpio, pe, 1)
          if      (size (y (lbpio(1):ubpio(1),lbpio(2):ubpio(2),:,:,:)) >0 ) &
            call p_send (y (lbpio(1):ubpio(1),lbpio(2):ubpio(2),:,:,:), pe, 1)
        else
          x  (lb(1):ub(1), lb(2):ub(2),:,:,:) = &
           y (lb(1):ub(1), lb(2):ub(2),:,:,:)
        end if
      end do
    else
      x = 0._wp
      call p_send (lb, dace% pio, 1)
      call p_send (ub, dace% pio, 1)
      if      (size (x(lb(1):ub(1), lb(2):ub(2),:,:,:)) > 0)        &
        call p_recv (x(lb(1):ub(1), lb(2):ub(2),:,:,:), dace% pio, 1)
    endif
  end subroutine scatter_coarse
!------------------------------------------------------------------------------
  subroutine read_infl_factor (rho, grid, fname)
  !--------------------------------------
  ! read adaptive inflation factor
  ! this is currently a simple ASCII file
  !--------------------------------------
  real(wp)         ,pointer    :: rho (:,:,:,:,:) ! inflation factor
  type(t_grid)     ,intent(in) :: grid            ! 'coarse' grid meta data
  character(len=*) ,intent(in) :: fname           ! file name

    real(wp), pointer  :: rho_g (:,:,:,:,:) ! full field to be read on 1 PE
    integer            :: iunit             ! Fortran file unit number
    integer            :: ierr              ! I/O error return variable
    logical            :: exists            ! check presence of file
    real(wp)           :: test              ! try to read file beyond EOF
    character(len=256) :: path              ! full pathname

    nullify (rho_g)
    if (dace% lpio) then
      allocate (rho_g (grid% lbg(1):grid% ubg(1),&
                       grid% lbg(2):grid% ubg(2),&
                       grid% lbg(3):grid% ubg(3),&
                       grid% lbg(4):grid% ubg(4),&
                       1                        ))
      iunit = get_unit_number()

      !------------------------
      ! check existence of file
      !------------------------
      path = path_file (input, fname)
      inquire (file=path, exist=exists)
!     if (.not. exists) then
!       !+++++++++++++++++++++++++++++
!       ! fallback to output directory
!       !+++++++++++++++++++++++++++++
!       call message ('read_infl_factor',                                  &
!                    'cannot open adaptive inflation file: '//trim(path)// &
!                    '; trying output directory'                           )
!       path = path_file (output, fname)
!       inquire (file=path, exist=exists)
!     endif
      !----------
      ! open file
      !----------
      ierr = -1
      if (exists) open (iunit, file=path, action='read', iostat=ierr)
      if (ierr /= 0) then
        call finish ('read_infl_factor',                                 &
                     'cannot open adaptive inflation file: '//trim(path))
      endif
      !----------
      ! read data
      !----------
      read(iunit,'(F4.2)', iostat=ierr) rho_g
      if (ierr /= 0) then
        call finish ('read_infl_factor',                                         &
                     'error while reading file: '//trim(path)//' (wrong size ?)')
      endif
      !--------------
      ! check for EOF
      !--------------
      read(iunit,'(F4.2)', iostat=ierr) test
      if (ierr == 0) then
        call finish ('read_infl_factor',                                     &
                     'error while reading file: '//trim(path)//' (too long)')
      endif

      close(iunit)
      call return_unit_number (iunit)
    endif
    call scatter_coarse (rho_g, rho, grid)
    if (dace% lpio) deallocate (rho_g)
  end subroutine read_infl_factor

!==========================================================================
  subroutine relax_prior_pert (enkf_ens, enkf_fc, par_trans, k_enkf, rtpp_alpha)
  !--------------------------------------
  ! add mean to ensemble members
  ! (only variables listed in par_trans!)
  !--------------------------------------
    type(t_atm)        ,pointer        :: enkf_ens(:)! ens to process
    type(t_atm)        ,pointer        :: enkf_fc(:)! ens to process
    character(len=*)   ,intent(in)     :: par_trans  ! parameters
    integer            ,intent(in)     :: k_enkf     ! number of ens members
    real(wp)           ,intent(in)     :: rtpp_alpha  ! weight of forecast perturbations

    !---------------
    !local variables
    !---------------
    integer                     :: k       ! ensemble member index
    character(len=16)           :: pars(nm) ! names of variables to derive
    integer                     :: n        ! return code
    integer                     :: m        ! variable number

    call split (pars, par_trans, n)
    do m=1, size(enkf_ens(1)% m)
      if(any(pars==enkf_ens(1)% m(m)% i% name)) then
          do k = 1, k_enkf
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
            enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
            (1._wp - rtpp_alpha) * enkf_ens(k)% m(m)% ptr(:,:,:,:) + &
            rtpp_alpha * enkf_fc(k)% m(m)% ptr(:,:,:,:)
          end do
      end if
    end do

  end subroutine relax_prior_pert
!==============================================================================
   subroutine relax_prior_spread (enkf_ens, enkf_fc_spread, enkf_an_spread, &
                                  par_trans, k_enkf, rtps_alpha, rtps_type  )
  !--------------------------------------
  ! add mean to ensemble members
  ! (only variables listed in par_trans!)
  !--------------------------------------
    type(t_atm)        ,pointer        :: enkf_ens(:)    ! ens to process
    type(t_atm)        ,intent(inout)  :: enkf_fc_spread ! forecast spread
    type(t_atm)        ,intent(inout)  :: enkf_an_spread ! analysis spread
    character(len=*)   ,intent(in)     :: par_trans      ! parameters
    integer            ,intent(in)     :: k_enkf         ! number of ens members
    real(wp)           ,intent(in)     :: rtps_alpha     ! weight of forecast perturbations
    integer            ,intent(in)     :: rtps_type      ! type of rtps applied

    !---------------
    !local variables
    !---------------
    integer                     :: k        ! ensemble member index
    character(len=16)           :: pars(nm) ! names of variables to derive
    integer                     :: n        ! return code
    integer                     :: m        ! variable number
    real(wp) ,parameter         :: min_infl = 1.e-9_wp  !  minimum value in case of zero spread

    select case (rtps_type)
    case(0)
      !-------------
      ! default case
      !-------------
      call split (pars, par_trans, n)
      do m=1, size(enkf_ens(1)% m)
        if (any(pars==enkf_ens(1)% m(m)% i% name)) then
          do k = 1, k_enkf
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
            enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
            enkf_ens(k)% m(m)% ptr(:,:,:,:) * &
            min( 100._wp, rtps_alpha * ( enkf_fc_spread  % m(m)% ptr(:,:,:,:) &
                                       - enkf_an_spread  % m(m)% ptr(:,:,:,:) )  / &
            (min_infl + enkf_an_spread  % m(m)% ptr(:,:,:,:) )  + 1 )
          end do
        end if
      end do

    case (1)
      !------------------------------
      ! italian met. center version 1
      !------------------------------
      call split (pars, par_trans, n)
      do m=1, size(enkf_ens(1)% m)
        if(any(pars==enkf_ens(1)% m(m)% i% name)) then
          select case (m)
          case (W_SO)
            do k = 1, k_enkf
              enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
              enkf_ens(k)% m(m)% ptr(:,:,:,:) * &
              min( 100._wp, max(1._wp,rtps_alpha * ( enkf_fc_spread  % m(m)% ptr(:,:,:,:) &
                                         - enkf_an_spread  % m(m)% ptr(:,:,:,:) )  / &
              (min_infl + enkf_an_spread  % m(m)% ptr(:,:,:,:) )  + 1 ))
            end do
          case (T, U, V, Q)
            do k = 1, k_enkf
              enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
              enkf_ens(k)% m(m)% ptr(:,:,:,:) * &
              min( 100._wp, max(1._wp,rtps_alpha * ((0.25*(enkf_fc_spread % m(T)% ptr(:,:,:,:) &
                                               + enkf_fc_spread % m(U)% ptr(:,:,:,:) &
                                               + enkf_fc_spread % m(V)% ptr(:,:,:,:) &
                                               + enkf_fc_spread % m(Q)% ptr(:,:,:,:)) - &
                                           0.25*(enkf_an_spread % m(T)% ptr(:,:,:,:) &
                                               + enkf_an_spread % m(U)% ptr(:,:,:,:) &
                                               + enkf_an_spread % m(V)% ptr(:,:,:,:) &
                                               + enkf_an_spread % m(Q)% ptr(:,:,:,:))))  / &
                               (min_infl + 0.25*(enkf_an_spread % m(T)% ptr(:,:,:,:) &
                                               + enkf_an_spread % m(U)% ptr(:,:,:,:) &
                                               + enkf_an_spread % m(V)% ptr(:,:,:,:) &
                                               + enkf_an_spread % m(Q)% ptr(:,:,:,:)))  + 1 ))
            end do
          end select
        end if
      end do

    case (2)
      !------------------------------
      ! italian met. center version 2
      !------------------------------
      call split (pars, par_trans, n)
      do m=1, size(enkf_ens(1)% m)
        if(any(pars==enkf_ens(1)% m(m)% i% name)) then
          select case (m)
          case (W_SO)
            do k = 1, k_enkf
              enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
              enkf_ens(k)% m(m)% ptr(:,:,:,:) * &
              min( 100._wp, max(1._wp,sqrt (rtps_alpha * ( enkf_fc_spread  % m(m)% ptr(:,:,:,:)**2 &
                                         - enkf_an_spread  % m(m)% ptr(:,:,:,:)**2 )  / &
              (min_infl + enkf_an_spread  % m(m)% ptr(:,:,:,:)**2 )  + 1) ))
            end do
          case (T, U, V, Q)
            do k = 1, k_enkf
              enkf_ens(k)% m(m)% ptr(:,:,:,:) = &
              enkf_ens(k)% m(m)% ptr(:,:,:,:) * &
              min( 100._wp, max(1._wp,sqrt (rtps_alpha * ((0.25*(enkf_fc_spread % m(T)% ptr(:,:,:,:)**2 &
                                               + enkf_fc_spread % m(U)% ptr(:,:,:,:)**2 &
                                               + enkf_fc_spread % m(V)% ptr(:,:,:,:)**2 &
                                               + enkf_fc_spread % m(Q)% ptr(:,:,:,:)**2) - &
                                           0.25*(enkf_an_spread % m(T)% ptr(:,:,:,:)**2 &
                                               + enkf_an_spread % m(U)% ptr(:,:,:,:)**2 &
                                               + enkf_an_spread % m(V)% ptr(:,:,:,:)**2 &
                                               + enkf_an_spread % m(Q)% ptr(:,:,:,:)**2)))  / &
                               (min_infl + 0.25*(enkf_an_spread % m(T)% ptr(:,:,:,:)**2 &
                                               + enkf_an_spread % m(U)% ptr(:,:,:,:)**2 &
                                               + enkf_an_spread % m(V)% ptr(:,:,:,:)**2 &
                                               + enkf_an_spread % m(Q)% ptr(:,:,:,:)**2))  + 1) ))
            end do
          end select
        end if
      end do

    case default
      call finish ("relax_prior_spread","unknown value of 'rtps_type'")
    end select

  end subroutine relax_prior_spread
!==========================================================================

  elemental subroutine sat_adj (ana)
    type (t_atm) ,intent(inout) :: ana
    !----------------------
    ! saturation adjustment
    !----------------------
    real(wp), allocatable  :: qv_sat(:,:,:,:)

    !-------------------
    ! allocate temporary
    !-------------------
    allocate (qv_sat(ana% grid% lb(1):ana% grid% ub(1), &
                     ana% grid% lb(2):ana% grid% ub(2), &
                     ana% grid% ub(3),                  &
                     ana% grid% lb(4):ana% grid% ub(4) ))
    !initialize because of cray optimization compiler bug
    qv_sat = 0._wp
!---------------------------------------------------------
!!!!!!!!!!! change qv and qcl/qci depending on temperature
!---------------------------------------------------------
!!$    do k=1, k_ens
!!$       !determine p_sat_qv -> qv_sat
!!$       qv_sat = q_e(es_t(enkf_an_ens(k)% t),ana% grid% p0 + enkf_an_ens(k)% pp)
!!$       !add qcl and qci in case of subsaturation
!!$       where ((enkf_an_ens(k)% q - qv_sat) < 0._wp)
!!$          enkf_an_ens(k)% q = enkf_an_ens(k)% q + enkf_an_ens(k)% qcl + enkf_an_ens(k)% qci
!!$          enkf_an_ens(k)% qcl = 0._wp
!!$          enkf_an_ens(k)% qci = 0._wp
!!$       end where
!!$       !correct/reduce qcl/qci in case of supersaturation
!!$       qv_sat = enkf_an_ens(k)% q - qv_sat
!!$       where (qv_sat > 0._wp)
!!$          where (enkf_an_ens(k)% t >= b3 )!T>=273.16 K
!!$             enkf_an_ens(k)% qcl = enkf_an_ens(k)% qcl + qv_sat
!!$          elsewhere                   !T<273.16 K
!!$             enkf_an_ens(k)% qci = enkf_an_ens(k)% qci + qv_sat
!!$          end where
!!$          enkf_an_ens(k)% q = enkf_an_ens(k)% q - qv_sat
!!$       end where
!!$    end do

    !---------------------------------------------------
    ! change only qv and qcl, do nothing if t < 233.15 K
    !---------------------------------------------------
    !---------------------------------
    ! check if temperature >= 233.15 K
    !---------------------------------
    where (ana% t >= 233.15_wp)
      !-----------------------------
      ! determine p_sat_qv -> qv_sat
      !-----------------------------
      qv_sat = q_e(esw_t(ana% t),ana% pf)
      !---------------------------------------------
      ! add qcl in case of subsaturation
      ! where ((ana% q - qv_sat) < 0._wp)
      !---------------------------------------------
      ana% q = ana% q + ana% qcl
      ana% qcl = 0._wp
      !----------------------------------------------
      ! end where
      ! correct/reduce qcl in case of supersaturation
      !----------------------------------------------
      qv_sat = ana% q - qv_sat
      where (qv_sat > 0._wp)
        !---------------------------------------------------
        ! ana% qcl = ana% qcl + qv_sat
        !---------------------------------------------------
        ana% qcl = qv_sat
        ana% q   = ana% q - qv_sat
      end where
    end where

    deallocate (qv_sat)

  end subroutine sat_adj
!==========================================================================
  elemental subroutine apply_hum_bounds (atm, rhmax, rhmin, qmin)
    type(t_atm) ,intent(inout)          :: atm     ! Atmospheric state
    real(wp)    ,intent(in)   ,optional :: rhmax   ! Upper limit on rel.hum.
    real(wp)    ,intent(in)   ,optional :: rhmin   ! Lower limit on rel.hum.
    real(wp)    ,intent(in)   ,optional :: qmin    ! Lower limit on spec.hum.

    real(wp) :: rh_max                             ! Local copy of rhmax
    real(wp) :: rh_min                             ! Local copy of rhmin
    real(wp) :: q_min                              ! Local copy of qmin

    rh_max = 1.1_wp; if (present (rhmax)) rh_max = rhmax
    rh_min = 0.0_wp; if (present (rhmin)) rh_min = rhmin
    q_min  = 0.0_wp; if (present (qmin))  q_min  = qmin

    if   (associated (atm% qci)) then
      if (associated (atm% qcl)) then
        where (atm% qci < 0._wp)
                      atm% qcl = atm% qcl + atm% qci
        endwhere
      endif
                      atm% qci = max (atm% qci, 0.0_wp)
    endif

    if (associated (atm% qcl)) then
      if (associated (atm% q)) then
        where (atm% qcl < 0._wp)
                    atm% q   = atm% q + atm% qcl
        endwhere
      endif
                    atm% qcl = max (atm% qcl, 0.0_wp)
    endif

    if (associated (atm% q  )) &
                    atm% q   = max (atm% q  , q_min )
    if (associated (atm% qr  )) &
                    atm% qr  = max (atm% qr , 0.0_wp)
    if (associated (atm% qs  )) &
                    atm% qs  = max (atm% qs , 0.0_wp)
    if (associated (atm% qg  )) &
                    atm% qg  = max (atm% qg , 0.0_wp)
    if (associated (atm% qh  )) &
                    atm% qh  = max (atm% qh , 0.0_wp)

    if (associated (atm% rh  )) then
                    atm% rh   = max (atm% rh  , rh_min)
                    atm% rh   = min (atm% rh  , rh_max)
    endif
    if (associated (atm% rh2m)) then
                    atm% rh2m = max (atm% rh2m, rh_min)
                    atm% rh2m = min (atm% rh2m, rh_max)
    endif

    if (associated (atm% nccloud  )) &
                    atm% nccloud  = max (atm% nccloud , 0.0_wp)
    if (associated (atm% ncrain   )) &
                    atm% ncrain   = max (atm% ncrain ,  0.0_wp)
    if (associated (atm% ncice    )) &
                    atm% ncice    = max (atm% ncice  ,  0.0_wp)
    if (associated (atm% ncsnow   )) &
                    atm% ncsnow   = max (atm% ncsnow ,  0.0_wp)
    if (associated (atm% ncgraupel)) &
                    atm% ncgraupel= max (atm% ncgraupel,0.0_wp)
    if (associated (atm% nchail   )) &
                    atm% nchail   = max (atm% nchail ,  0.0_wp)

  end subroutine apply_hum_bounds
!==========================================================================
  elemental subroutine apply_tke_bounds (atm)
    type (t_atm) ,intent(inout) :: atm
    atm% tke = max (atm% tke, 0.00004_wp)
  end subroutine apply_tke_bounds
!==========================================================================
! GME specific auxiliary routines:
!==============================================================================
  subroutine set_d_ps (datm_dps, enkf_fc_mean, enkf_fc)
  type(t_atm) ,intent(out)   :: datm_dps       ! impact of ps on atm
  type(t_atm) ,intent(inout) :: enkf_fc_mean   ! forecast ensemble mean
  type(t_atm) ,intent(inout) :: enkf_fc (:)    ! forecast ensemble
  !-------------------------------------------------------------------------
  ! 1) Determine the impact of surface pressure variations (in the hybrid
  !    vertical coordinate system) on atmospheric quantities interpolated to
  !    a fixed grid. Allthough this could be done analytically (by automatic
  !    differentiation) this is done here by a finite difference approach.
  !
  !    The impact of surface pressure variations is calculated for the
  !    ensemble mean only. Processing each member individually might be more
  !    accurate but requires more resources.
  !
  ! 2) Modify the forecast ensemble, so that the effect of surface pressure
  !    variations due to a shift of the hybrid vertical coordinate system is
  !    not included in the correlations derived from the ensemble.
  !
  ! This routine should be applied in conjunction with subroutine apply_d_ps
  ! GME-specific
  !-------------------------------------------------------------------------
    type(t_atm) :: tmp(2) ! temporary for finite differencing
    integer     :: k      ! level index
    integer     :: m      ! atmospheric field index
    integer     :: k_ens  ! ensemble member index
    !-----------------------------
    ! return if not on hybrid grid
    !-----------------------------
    if (enkf_fc_mean% grid% vct /= VCT_P_HYB) return
!   if (fix_plev <= 0)                        return
    !------------------------------------------
    ! determine gradient by finite differencing
    !------------------------------------------
    if (dace% lpio) write(6,'(/"  set_d_ps:"/)')
    call construct (tmp(1)   ,template=enkf_fc_mean)
    call construct (tmp(2)   ,template=enkf_fc_mean)
    call construct (datm_dps ,template=enkf_fc_mean)
    call allocate  (tmp, 'psr')
    tmp(2)% psr = enkf_fc_mean% ps
    tmp(1)% psr = enkf_fc_mean% ps
    enkf_fc_mean% ps = tmp(1)% psr + 50._wp  ! + 0.5 hPa
    call vert_intp (enkf_fc_mean, tmp(2))
    enkf_fc_mean% ps = tmp(1)% psr - 50._wp  ! - 0.5 hPa
    call vert_intp (enkf_fc_mean, tmp(1))
    enkf_fc_mean% ps = tmp(1)% psr
    datm_dps    = (tmp(2) - tmp(1)) / 100._wp
    call print (datm_dps, comment='set_d_ps: datm_dps')
    call deallocate(datm_dps, 'ps')
    call destruct  (tmp)
    !-----------------------------
    ! modify the forecast ensemble
    !-----------------------------
    if (dace% lpio) write(6,'(/"  set_d_ps: modify forecast ensemble"/)')
    do k_ens = 1, size (enkf_fc)
      do m = 1, size (datm_dps% m)
        if (.not. ALLOCATED  (datm_dps      % m(m)% ptr)  ) cycle
        if (.not. ALLOCATED  (enkf_fc(k_ens)% m(m)% ptr)  ) cycle
        if (                  datm_dps      % m(m)% i% ref) cycle
        if (dace% lpio .and. k_ens == 1)                        &
          write(6,'("  set_d_ps: ",a)') datm_dps% m(m)% i% name
        do k = lbound(datm_dps% m(m)% ptr,3), ubound(datm_dps% m(m)% ptr,3)
          enkf_fc(k_ens)% m(m)% ptr(:,:,k,:) = &
          enkf_fc(k_ens)% m(m)% ptr(:,:,k,:) + &
          enkf_fc(k_ens)%       ps (:,:,1,:) * &
          datm_dps      % m(m)% ptr(:,:,k,:)
        end do
      end do
    end do
    if (dace% lpio) write(6,'(/"  set_d_ps: return"/)')
  end subroutine set_d_ps
!==============================================================================
  subroutine apply_d_ps (enkf_an, datm_dps)
  type(t_atm) ,intent(inout) :: enkf_an (:)    ! analysis ensemble mean
  type(t_atm) ,intent(inout) :: datm_dps       ! impact of ps on atm
  !-----------------------------------------------------------------------
  ! Account for the impact of surface pressure variations (due to the
  ! shift in the vertical hybrid coordinate system) on atmospheric
  ! quantities for each ensemble member.
  !
  ! This routine should be applied in conjunction with subroutine set_d_ps
  ! GME-specific
  !-----------------------------------------------------------------------
    integer     :: k      ! level index
    integer     :: m      ! atmospheric field index
    integer     :: k_ens  ! ensemble member index
    !-----------------------------
    ! return if not on hybrid grid
    !-----------------------------
!   if (fix_plev <= 0)                    return
    if (datm_dps% grid% vct /= VCT_P_HYB) return
    if (dace% lpio) write(6,'(/"  apply_d_ps:"/)')
    !-----------------------------
    ! modify the analysis ensemble
    !-----------------------------
    if (dace% lpio) write(6,'(/"  apply_d_ps: modify analysis ensemble"/)')
    do k_ens = 1, size (enkf_an)
      do m = 1, size (datm_dps% m)
        if (.not. ALLOCATED  (datm_dps      % m(m)% ptr)  ) cycle
        if (.not. ALLOCATED  (enkf_an(k_ens)% m(m)% ptr)  ) cycle
        if (                  datm_dps      % m(m)% i% ref) cycle
        if (dace% lpio .and. k_ens == 1) &
          write(6,'("  apply_d_ps: ",a)') datm_dps% m(m)% i% name
        do k = lbound(datm_dps% m(m)% ptr,3), ubound(datm_dps% m(m)% ptr,3)
          enkf_an(k_ens)% m(m)% ptr(:,:,k,:) = &
          enkf_an(k_ens)% m(m)% ptr(:,:,k,:) - &
          enkf_an(k_ens)%       ps (:,:,1,:) * &
          datm_dps      % m(m)% ptr(:,:,k,:)
        end do
      end do
    end do
    !--------------------
    ! deallocate gradient
    !--------------------
    call destruct (datm_dps)
    if (dace% lpio) write(6,'(/"  apply_d_ps: return"/)')
  end subroutine apply_d_ps
!=============================================================================
  subroutine to_common_grid (enkf_fc, enkf_fc_mean)
  type(t_atm) ,intent(inout) :: enkf_fc (:)    ! first k_enkf ens members
  type(t_atm) ,intent(in)    :: enkf_fc_mean   ! forecast ensemble mean
  !--------------------------------------------------------------
  ! interpolate to common vertical grid (p-levs of forecast mean)
  ! (GME-specific:currently hybrid pressure coordinates only)
  !--------------------------------------------------------------
    integer      :: k    ! ensemble member index
    type (t_atm) :: tmp  ! temporary atmospheric state
    !-----------------------------
    ! return if not on hybrid grid
    !-----------------------------
    if (enkf_fc_mean% grid% vct /= VCT_P_HYB) return
!   if (fix_plev >= 0)                        return
    !--------------------
    ! construct temporary
    !--------------------
    call construct (tmp, enkf_fc_mean% grid)
    !---------------------------
    ! loop over ensemble members
    !---------------------------
    do k = 1, size (enkf_fc)
      !-----------------------
      ! set reference pressure
      !-----------------------
      call allocate (tmp, 'psr'); tmp% psr = enkf_fc_mean% ps
      !------------
      ! interpolate
      !------------
      call vert_intp (enkf_fc(k), tmp, inplace=.true.)
      !----------
      ! copy back
      !----------
      enkf_fc(k) = tmp
      call deallocate (enkf_fc(k), 'pf')
      call deallocate (enkf_fc(k), 'ph')
      call deallocate (enkf_fc(k), 'psr')
    end do
    !---------------------
    ! deallocate temporary
    !---------------------
    call destruct (tmp)

  end subroutine to_common_grid
!==============================================================================
  subroutine from_common_grid (enkf_an, enkf_fc_mean)
  type(t_atm) ,intent(inout) :: enkf_an (:)    ! first k_enkf ens members
  type(t_atm) ,intent(in)    :: enkf_fc_mean   ! forecast ensemble mean
  !--------------------------------------------------------------
  ! interpolate to common vertical grid (p-levs of forecast mean)
  ! (GME-specific:currently hybrid pressure coordinates only)
  !--------------------------------------------------------------
    integer      :: k    ! ensemble member index
    type (t_atm) :: tmp  ! temporary atmospheric state
    !-----------------------------
    ! return if not on hybrid grid
    !-----------------------------
    if (enkf_fc_mean% grid% vct /= VCT_P_HYB) return
!   if (fix_plev >= 0)                        return
    !--------------------
    ! construct temporary
    !--------------------
    call construct (tmp, enkf_fc_mean% grid)
    !---------------------------
    ! loop over ensemble members
    !---------------------------
    do k = 1, size (enkf_an)
      !-----------------------
      ! set reference pressure
      !-----------------------
      call allocate (enkf_an(k), 'psr')
      call allocate (tmp,        'psr')
      enkf_an(k)% psr = enkf_fc_mean% ps
      tmp%        psr = enkf_an(k)% ps
      !------------
      ! interpolate
      !------------
      call vert_intp (enkf_an(k), tmp, inplace=.true.)
      !----------
      ! copy back
      !----------
      enkf_an(k) = tmp
      call deallocate (enkf_an(k), 'pf')
      call deallocate (enkf_an(k), 'ph')
      call deallocate (enkf_an(k), 'psr')
    end do
    !---------------------
    ! deallocate temporary
    !---------------------
    call destruct (tmp)
  end subroutine from_common_grid
!==============================================================================
  subroutine adjust_soil (ens, adj, adj_wso)
    !-------------------------------------------------
    ! adjust ensemble soil moisture (w_so) to analysis
    !-------------------------------------------------
    type(t_atm) ,intent(inout) :: ens  (:)   ! ensemble state
    type(t_atm) ,intent(in)    :: adj        ! adjustment
    real(wp)    ,intent(in)    :: adj_wso(:) ! nudging factor

    integer  :: i,j,d                          ! horizontal index
    integer  :: k                              ! ensemble index
    integer  :: l                              ! layer index
    integer  :: n                              ! ensemble size
    integer  :: lb(4), ub(4)                   ! local bounds
    real(wp) :: f   (maxls)                    ! adjustment factor
    real(wp) :: old (ens(1)% grid% shape(1), & ! w_so before adjustment
                     ens(1)% grid% shape(2), &
                     ens(1)% grid% shape(4)  )
    real(wp) :: upp (ens(1)% grid% lb(1):ens(1)% grid% ub(1), &
                     ens(1)% grid% lb(2):ens(1)% grid% ub(2), &
                     ens(1)% grid% lb(4):ens(1)% grid% ub(4)  )
    real(wp) :: low (ens(1)% grid% lb(1):ens(1)% grid% ub(1), &
                     ens(1)% grid% lb(2):ens(1)% grid% ub(2), &
                     ens(1)% grid% lb(4):ens(1)% grid% ub(4)  )
    real(wp), pointer :: soiltyp (:,:,:)

    !------------------------------------------------------------
    ! set adjustment factor, 'adj_wso' is related to 1 day period
    !------------------------------------------------------------
    f = 1._wp - adj_wso
    f = max (f, tiny (1._wp))          ! safeguard for adj_wso == 1
    f = f ** (fc_hours / 24._wp)
    f = 1._wp - f

    !--------------------------------------------------------
    ! set upper (field capacity), lower (wilting point) bound
    !--------------------------------------------------------
    if (.not. associated (ens(1)% grid% soiltyp))                           &
      call finish ('adjust_soil','soil type not present in invariant fields')

    lb = ens(1)% grid% lb
    ub = ens(1)% grid% ub
    soiltyp (lb(1):, lb(2):, lb(4): ) => ens(1)% grid% soiltyp (lb(1):ub(1), &
                                                                lb(2):ub(2), &
                                                             1, lb(4):ub(4)  )

    do d = lb(4), ub(4)
    do j = lb(2), ub(2)
    do i = lb(1), ub(1)
      if (soiltyp (i,j,d) == 9999._wp) cycle       ! skip missing value
      upp (i,j,d) = Cfcap (int (soiltyp (i,j,d)))
      low (i,j,d) = Cpwp  (int (soiltyp (i,j,d)))
    end do
    end do
    end do

    n = size(ens)
    do k = 1, n                                    ! ensemble loop
      do l = 1, ens(1)% grid% ns                   ! layer    loop
        where (ens(k)% t_so (:,:,l,:) >  t0c     & ! not for frozen soil
              .and. soiltyp (:,:,  :) /= 9999._wp) ! not missing value
          !---------------------------
          ! adjust, remember old value
          !---------------------------
          old                    = ens(k)% w_so (:,:,l,:)
          ens(k)% w_so (:,:,l,:) = ens(k)% w_so (:,:,l,:) &
                             + f(l) * adj% w_so (:,:,l,:)
          !------------------------
          ! keep w_so within bounds
          !------------------------
          ens(k)  % w_so (:,:,l,:) = min (                  &
            ens(k)% w_so (:,:,l,:),  max (old, upp * st (l)))
          ens(k)  % w_so (:,:,l,:) = max (                  &
            ens(k)% w_so (:,:,l,:),  min (old, low * st (l)))
        endwhere
      end do
    end do
  end subroutine adjust_soil
!==============================================================================
  subroutine inflate_wso_spread (x, infl_fac, stdv_limit)
    type(t_atm) ,intent(inout) :: x         (:) ! ensemble state
    real(wp)    ,intent(in)    :: infl_fac  (:) ! inflation factor
    real(wp)    ,intent(in)    :: stdv_limit(:) ! ensemble stdv. limiter

    integer  :: l                         ! ensemble size
    integer  :: ne                        ! number of relevant ensemble members
    integer  :: i, j, k, d                ! indices
    integer  :: m                         ! soil level indices
    integer  :: st                        ! soil type
    integer  :: lb(4), ub(4)              ! local bounds
    logical  :: snow    (size(x))         ! snow covered pixel
    logical  :: invalid (size(x))         ! snow covered or frozen
    logical  :: valid_l (n_sl)            ! valid level
    real(wp) :: c       (n_st)            ! soil capacity
    real(wp) :: lower   (n_st)            ! lower limit on smi
    real(wp) :: upper   (n_st)            ! upper limit on smi
    real(wp) :: stmm    (n_sl)            ! soil thickness in mm
    real(wp) :: wso_0   (n_sl,n_st)       ! wso at smi=0
    real(wp) :: w       (n_sl)            ! default inflation
    real(wp) :: f       (n_sl)            ! inflation factor
    real(wp) :: spread_lim(n_sl)          ! smi spread limit
    real(wp) :: wso                       ! soil humidity
    real(wp) :: wso_mean(n_sl)            ! soil humidity mean
    real(wp) :: wso_max                   ! max. soil humidity of ensemble
    real(wp) :: wso_min                   ! min. soil humidity of ensemble
    real(wp) :: mean, spread              ! mean, spread
    real(wp) :: smi_emean, smi_spread     ! ensemble mean, spread
    real(wp) :: smi_emax,  smi_emin       ! ensemble max, min
    real(wp) :: fac, tmp, tmp1, tmp2      ! temporaries

    !----------------------------------
    ! check presence of required fields
    !----------------------------------
    if (.not. associated (x(1)% w_so )) then
      if (dace% lpio) then
        write (6,*) 'inflate_wso_spread: W_SO not present!'
        write (0,*) 'inflate_wso_spread: W_SO not present!'
      endif
      return
    endif

    if (size (infl_fac) < n_sl-1) then
       call finish ("inflate_wso_spread","size of array infl_fac too small")
    end if
    if (size (stdv_limit) < n_sl-1) then
       call finish ("inflate_wso_spread","size of array stdv_limit too small")
    end if

    m = x(1)% grid% ns
    if (m < n_sl-1) then
      if (dace% lpio) then
        write (6,*) 'inflate_wso_spread: found', m, '<', n_sl, 'soil levels!'
        write (0,*) 'inflate_wso_spread: found', m, '<', n_sl, 'soil levels!'
      end if
      call finish ("inflate_wso_spread","incompatible soil model?")
    end if

    w                =   1._wp  ! default: no inflation or deflation
    spread_lim       = 999._wp  ! default: no limit
    spread_lim(n_sl) =   0._wp  ! do not touch climatological layer

    do m = 1, n_sl-1
       if (infl_fac  (m) > 1._wp) w(m)          = infl_fac  (m)
       if (stdv_limit(m) > 0._wp) spread_lim(m) = stdv_limit(m)
    end do

    c    = cfcap - cpwp                   ! soil capacity for smi=1
    stmm = soil_thickness (1:) * 1000._wp ! thickness to mm
    !------------------------------------------------------------------------
    ! smi sanity upper and lower limits at pore volume and air dryness point;
    ! wso_0: soil moisture at smi=0
    !------------------------------------------------------------------------
    lower = - HUGE (0._wp)
    upper =   HUGE (0._wp)
    wso_0 =   HUGE (0._wp)
    do st = 1, n_st
       select case (st)
       case (:0, ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE, n_st+1:)
          cycle
       case default
          upper(st) = (cporv(st)-cpwp(st)) / c(st)
          lower(st) = (cadp (st)-cpwp(st)) / c(st)
          do m = 1, n_sl - 1
             wso_0(m,st) = stmm(m) * cpwp(st)
          end do
       end select
    end do

    l  = size (x)               ! ensemble size
    lb = x(1)% grid% lb
    ub = x(1)% grid% ub
    !---------------------
    ! horizontal grid loop
    !---------------------
    do d = lb(4), ub(4)
    do j = lb(2), ub(2)
    do i = lb(1), ub(1)

      !----------
      ! land only
      !----------
      if (x(1)%grid% lsm(i,j,1,d) <= 0.5_wp) cycle
      st = nint (x(1)% grid% soiltyp (i,j,1,d))
      select case (st)
      case (:0, ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE, n_st+1:)
         cycle
      end select
      !---------------------
      ! check for snow cover
      !---------------------
      snow = .false.
      do k = 1, l
        if (associated (x(k)% h_snow)) snow(k) = (x(k)% h_snow(i,j,1,d) > 0._wp)
      end do
      !-------------------------------
      ! check if spread is above limit
      !-------------------------------
      invalid  = snow
      valid_l  = .true.
      f        = w
      wso_mean = 0._wp
      !-----------------
      ! loop over levels
      !-----------------
      do m = 1, n_sl - 1
        !-----------------------------------------------
        ! calculate mean and spread of non-frozen points
        !-----------------------------------------------
        mean    = 0._wp
        spread  = 0._wp
        ne      = 0
        wso_min =   HUGE (0._wp)
        wso_max = - HUGE (0._wp)
        do k = 1, l
          if (x(k)% t_so (i,j,m,d) < t0c) invalid(k) = .true.
          if (invalid(k)) cycle
          wso     = x(k)% w_so (i,j,m,d)
          mean    = mean   + wso
          spread  = spread + wso ** 2
          wso_min = min (wso, wso_min)
          wso_max = max (wso, wso_max)
          ne      = ne   + 1
        end do
        if (ne == 0) cycle
        mean        =            mean   / ne
        spread      = sqrt (max (spread / ne - mean ** 2, 0._wp))
        wso_mean(m) = mean
        !---------------
        ! convert to smi
        !---------------
        fac = 1._wp / (stmm(m) * c(st))
        smi_emin   = (wso_min - wso_0(m,st)) * fac
        smi_emax   = (wso_max - wso_0(m,st)) * fac
        smi_emean  = (mean    - wso_0(m,st)) * fac
        smi_spread =  spread                 * fac

        if (spread_lim(m) <= 0._wp .or. smi_spread <  1.e-12_wp    &
                                   .or. smi_spread >= spread_lim(m)) then
          valid_l(m) = .false.
        else
          !------------------------------------
          ! determine maximum allowed inflation
          !------------------------------------
          tmp1 = max (smi_emax  - smi_emean, 1.e-12_wp)
          tmp2 = max (smi_emean - smi_emin,  1.e-12_wp)
          select case (m)
          case (1:2)
            ! Layers 1-2: allow range from air dryness point to pore volume
            tmp = min ((upper(st) - smi_emean) / tmp1, &
                       (smi_emean - lower(st)) / tmp2  )
          case (3:)
            ! Layers 3-7: limit range to 0 (wilt.point) <= smi <= 1 (field cap.)
            tmp = min ((1._wp     - smi_emean) / tmp1, &
                       (smi_emean - 0._wp    ) / tmp2  )
          end select
          !---------------------------------
          ! derive maximum desired inflation
          !---------------------------------
          tmp = min (tmp, spread_lim(m) / smi_spread, w(m))

          if (tmp <= 1._wp) then
            valid_l(m) = .false.
          else
            f(m) = tmp
          end if
        endif
      end do ! m
      !---------------------
      ! ensemble member loop
      !---------------------
      do k = 1, l
        if (snow(k)) cycle
        !-----------------
        ! loop over levels
        !-----------------
        do m = 1, n_sl - 1
          !------------------------------
          ! no action below frozen layers
          ! skip if limits exhausted
          !------------------------------
          if (x(k)% t_so (i,j,m,d) < t0c) exit
          if (.not. valid_l(m))           cycle
          !------------------------------------------------------------
          ! scale soil moisture differences to mean by inflation factor
          !------------------------------------------------------------
          x(k)%       w_so (i,j,m,d) = wso_mean(m)  +    &
               (x(k)% w_so (i,j,m,d) - wso_mean(m)) * f(m)
        end do
      end do

    end do ! i
    end do ! j
    end do ! d

  end subroutine inflate_wso_spread
!==============================================================================
end module mo_t_enkf
