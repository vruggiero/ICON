!
!+ Interface from 3D-Var to radiance bias-correction module
!
module mo_radbias_3dv
!
! Description:
!   Interface from 3D-Var to radiance bias-correction module
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2009/12/09 Andreas Rhodin
!  Interface from 3D-Var to radiance bias-correction module
! V1_9         2010/04/20 Andreas Rhodin
!  3dvar-interface to bias-correction module
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  bug fix (online bias correction using IWV predictor)
!  option to calculate bias correction statistics on obs-ana
! V1_13        2011/11/01 Mashrab Kuvatov
!  changes for IASI, HIRS, AMSUB
! V1_19        2012-04-16 Robin Faulwetter
!  Improved read_fdbk to enable rad_biascor to deal with multiple grids
! V1_20        2012-06-18 Andreas Rhodin
!  fix memory leak
! V1_22        2013-02-13 Robin Faulwetter
!  implementation of vectorized K-mode
!  changes for HIRS cloud check (O.Stiller)
! V1_23        2013-03-26 Robin Faulwetter
!  Improve bias correction modules: t_decay = 0. : infinite memory
!                                   t_decay < 0. : static bias correction.
! V1_26        2013/06/27 Harald Anlauf
!  Write diagnostics files to output directories, not the current one
! V1_27        2013-11-08 Andreas Rhodin
!  implementation of Variational Bias Correction (VarBC)
! V1_28        2014/02/26 Andreas Rhodin
!  VarBC: fix for zero number of observations
!  namelist variable 'precon': account for VarBC in preconditioner
! V1_31        2014-08-21 Robin Faulwetter
!  Unify mo_rad with COSMO. Improve mo_rttov_ifc. New write_rttov_prof
! V1_37        2014-12-23 Robin Faulwetter
!  Cleanup
! V1_44        2015-09-30 Harald Anlauf
!  Cleanup (mdlsfc is integer, not real)
! V1_47        2016-06-06 Robin Faulwetter
!  Many improvements for radiances.
! V1_48        2016-10-06 Andreas Rhodin
!  Old 1dvar-style static bias correction is not enabled by default any more
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured/unified cloud detection for radiances.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

!-------------
! Modules used
!-------------
use mo_kind,         only: sp,                  &! single  precision kind
                           wp,                  &! working precision kind
                           i1                    ! 1 byte integer kind
use mo_exception,    only: finish                ! abort routine
use mo_time,         only: cyyyymmddhhmm,       &! convert time to string
                           time_cyyyymmddhhmm,  &! convert string to time
                           operator(-),         &! subtract times
                           days                  ! convert time to days
use mo_mpi_dace,     only: dace,                &! MPI group info
                           p_bcast,             &! broadcast routine
                           p_bcast_ptr,         &! broadcast pointer
                           p_sum,               &! sum over processor elements
                           p_or                  ! .or. over PEs
use mo_biasc_io,     only: bc_filename,         &! derive file name
                           bc_path,             &! common path (optional)
                           bc_paths,            &! full pathnames
                           bc_obstyp,           &! observation type in files
                           nbcf,                &! number of files in list
                           fallback              ! use template if no input file
use mo_rad,          only: t_radv,              &
                           t_rad_set,           &
                           rinvalid,            &! invalid value
                           set_indx,            &! get bias correction set index
                           chan_indx,           &! get channel index
                           print_rad_set,       &!
                           rad_set,             &!
                           n_set
use mo_radbiascor,   only: bccoef,              &! coefficients
                           t_bcor_coef,         &! type to hold coefficients
                           t_stats_data,        &! type to hold statistics
                           assignment (=),      &! t_rad_set = t_rad_set
                           read_bcor_file,      &! read coefficients
                           write_bcor_file,     &! write coefficients
                           construct,           &! allocate   t_stats_data% ..
                           destruct,            &! deallocate t_stats_data% ..
                           calc_bcor_coefs,     &! calculate bias corr.coeff.
                           calc_bcor_coefs_ac,  &! bias corr.coeff. from ac.stat
                           calc_bcor_t_stats,   &! calculate bias corr. stats
                           calc_predictors,     &! calculate predictors
                           merge_bcor_stats,    &! merge bias corr. stats
                           scale_bcor_stats,    &! change size of statistics
                           apply_bias_corr,     &! apply bias correction
                           atm_thick_vers,      &! testing
                           min_entries,         &! min. number of entries
                           offset_max_age,      &! max. age of stat. in offset calc.
                           offset_t_decay,      &! t_decay of offset
                           P_THICK_TR,          &! thickness weighted with transmission delta in layer
                           P_THICK_TRF,         &! thickness weighted with top-of-layer transmission
                           P_TRANS_T,           &! T integrated over transmission interval
                           P_TRANS_TS,          &! T integrated over transmission interval
                           P_TRANS_P,           &! pressure of transmission value
                           P_TRANS_DP,          &! pressure difference between transmission values
                           chan_dep_pred
use mo_namelist,     only: position_nml,        &! position namelist
                           nnml,                &! namelist Fortran unit number
                           POSITIONED            ! ok    code from position_nml
use mo_fortran_units,only: get_unit_number,     &! reserve a unit number
                           return_unit_number    ! release the unit number
use mo_run_params,   only: input,               &! path name of input  files
                           output,              &! path name of output files
                           path_file,           &! concatenate: path/file.sufx
                           run_time,            &! time of run of this program
                           ana_time,            &! analysis time
                           nex,                 &! experiment id
                           flag_biasc_rad        ! flag for online biascorr.
use mo_obs_set,      only: t_obs_set             ! observation data derived type
use mo_t_obs,        only: t_spot,              &! report/fov  data derived type
                           t_obs,               &! observation data derived type
                           ldeb, usd, dpref      ! Debug selected spot(s)
use mo_dec_matrix,   only: t_vector              ! vector data type
use mo_tovs_prof,    only: fill_rad,            &! fill t_radv from t_obs_set
                           valid                 ! check whether spot is valid
use mo_instrid,      only: rttov_instr,         &! RTTOV from WMO instrument id
                           instr_rttov           ! WMO from RTTOV instrument id
use mo_t_use,        only: STAT_PASSIVE,        &! not used, only monitored
                           STAT_REJECTED,       &! not used, suspicious quality
                           STAT_DEFAULT,        &! invalid value flag
                           STAT_PAS_REJ,        &! not used, suspicious quality, passive
                           CHK_BIASCOR,         &! flag for no bias correction
                           CHK_CLOUD,           &! flag for cloudy pixel/channel
                           CHK_INSDAT,          &! flag for insufficient data
                           decr_use,            &! decrease the state of a datum
                           stats,               &! status value table
                           chks => chk,         &! status flags table
                           change_bits,         &! convert 3dvar -> feedback
                           reverse_code,        &! convert code back
                           reverse_bits          ! convert bit-field back
use mo_fdbk,         only: t_fdbk,              &! feedback file data type
                           setup_fdbk,          &! set up feedback-file data structure
                           open_fdbk_read      ,&! open feedback file to read
                           open_fdbk_write     ,&! open feedback file to write
                           read_meta           ,&! read verification meta data
                           get_veri_index      ,&! return indices of verification runs
                           get_veri            ,&! read verification entry
                           add_history,         &! add history entry
                           close_fdbk,          &! close feedback file
                           cleanup_fdbk          ! deallocate components
use mo_fdbk_tables,  only: VT_FIRSTGUESS,       &! first guess flag
                           VT_ANALYSIS,         &! analysis    flag
                           VE_DETERM,           &! detrministic run flag
                           OT_RAD                ! Radiances report type ID
use mo_rad,          only: t_radv                ! Type for radiance data
use netcdf,          only: nf90_open,      &
                           nf90_close,     &
                           nf90_inq_varid, &
                           nf90_put_var,   &
                           NF90_NOWRITE,   &
                           NF90_NOERR
use mo_t_netcdf,     only: ncid,           &
                           strnc,          &
                           stanc,          &
                           counc,          &
                           get_dim,        &
                           get_attr,       &
                           get_var,        &
                           chk
implicit none

!----------------
! Public entities
!----------------
private
public :: radbias_init     ! initialize module: read namelist, biascor.coeff.
public :: radbias_bfg      ! apply bias correction, called before fg-check
public :: radbias_afg      ! update bias correction coefs., after fg-check
public :: read_fdbk        ! read 3dvar feedback file
public :: write_tbcor_fdbk ! restore bias correction in feedback file
public :: biascor_mode     ! bias correction mode
public :: no_bias          ! for testing
public :: precon           ! account for VarBC in preconditioner
public :: BC_NOBC,        &! no bias correction
          BC_1DVAR,       &! keep correction from 1D-Var
          BC_3DVAR,       &! apply  bias corr. in 3D-Var
          BC_UP_FG,       &! update bias corr. in 3D-Var
          BC_UP_AN,       &! update bias corr. from analysis
          BC_VARBC         ! variational bias correction

!---------------
! biascor_modes:
!---------------
integer ,parameter :: BC_NOBC  =  0 ! flag value: no bias correction
integer ,parameter :: BC_1DVAR =  1 !             keep correction from 1D-Var
integer ,parameter :: BC_3DVAR =  2 !             apply  bias corr. in 3D-Var
integer ,parameter :: BC_UP_FG =  4 !             update bias corr. in 3D-Var
integer ,parameter :: BC_UP_AN =  8 !             update bias corr. from ana.
integer ,parameter :: BC_VARBC = 16 !             variational bias correction

!--------------------
! Namelist RADBIASCOR
!--------------------
character(len=256) :: coeff_ifile = 'bias.out' ! coefficient input  file
character(len=256) :: coeff_ofile = 'bias.new' ! coefficient output file
integer            :: biascor_mode= BC_1DVAR   ! bias correction mode
logical            :: update_ana  = .true.     ! use updated coefs. in analysis
integer            :: format      = 0          ! 0: use old naming convention
                                               ! 1: write 1 file
                                               ! 2: write 1 file per satellite
                                               ! 3: write 1 file per instrument
logical            :: no_pred     = .false.    ! TESTING: deactivate predictors
logical            :: no_bias     = .false.    ! TESTING: deact. initial biasc.
logical            :: precon      = .true.     ! account for VarBC in precond.

namelist /RADBIASCOR/ coeff_ifile, coeff_ofile, format, biascor_mode, &
                      update_ana, min_entries, no_pred, no_bias, precon,&
                      atm_thick_vers, offset_max_age, offset_t_decay
!-----------------
! Module variables
!-----------------
  !
  ! 1DVAR bias correction coefficients for all satellites
  !

!==============================================================================
contains
!==============================================================================

  subroutine radbias_bfg (obs, y, x, mask_rs)
    !------------------------------------------------------------------------
    ! Radiance bias correction routine to be called before first guess check.
    ! Applies bias correction to observations.
    !------------------------------------------------------------------------
    type (t_obs_set) ,intent(inout)        :: obs  ! observation
    type (t_vector)  ,intent(in)           :: y    ! background, observation   space
    type (t_vector)  ,intent(in), optional :: x    ! background, interpolation space
    logical          ,intent(in), optional :: mask_rs(:) ! mask rad_set, i.e. do not process all data

    character(len=*), parameter :: proc = 'radbias_bfg'
    integer                :: ib      ! observation 'box'      index
    integer                :: is      ! observation report/fov index
    integer                :: i_rs    ! index in rad_set
    type (t_spot) ,pointer :: s       ! pointer to report/fov data
    real(sp)      ,pointer :: obc(:)  ! pointer to bias corrected obs.
    real(sp)      ,pointer :: bc (:)  ! pointer to bias correction value
    logical,      allocatable :: mask_bc(:)
    type(t_radv), allocatable :: rad    (:)
    integer                   :: nbc

    if (biascor_mode == 0) then
      !--------------------------------------
      ! remove bias correction applied so far
      !--------------------------------------
      do ib = 1,size(obs% o)                    ! loop over observation 'boxes'
        do is=1,obs% o(ib)% n_spot              ! loop over reports/fovs
          s => obs% o(ib)%spot(is)              !
          if (s% hd% obstype /= OT_RAD) cycle   ! only treat radiances
          obc => obs% o(ib)% body(s%o%i+1:s%o%i+s%o%n)% o
          bc  => obs% o(ib)% body(s%o%i+1:s%o%i+s%o%n)% bc
          obc = obc - bc                        ! reset corrected value
          bc  = 0._sp                           ! reset bias correction
        end do
      end do
    else if (iand (biascor_mode, BC_1DVAR) == BC_1DVAR) then
      !--------------------------------------------------------------
      ! keep bias correction from 1D-Var feedback-file, nothing to do
      !--------------------------------------------------------------
    else if (iand (biascor_mode, BC_3DVAR) == BC_3DVAR) then
      !-----------------------------------
      ! apply bias-correction within 3dvar
      !-----------------------------------
      if (associated(bccoef)) then
        nbc = size(bccoef)
      else
        nbc = 0
      end if
      if (nbc == 0) return
      allocate(mask_bc(nbc),rad(nbc))
      if (present(mask_rs)) then
        mask_bc = .false.
        do is = 1, nbc
          i_rs = set_indx(rad_set(1:n_set), satid=bccoef(is)%i%satid, grid=bccoef(is)%i%grid)
          if (i_rs > 0) mask_bc(is) = mask_rs(i_rs)
        end do
        if (count(mask_bc) /= count(mask_rs)) call finish(proc,&
             'Did not find biascorrection for all datasets')
      else
        mask_bc = .true.
      end if
      call fill_rad(rad, bccoef(:)% i, obs, x, y, cp_opts=.true., lbc=.true., mask_rs=mask_bc)
      call calc_bcor_coefs_ac   (bccoef, mask_rs=mask_bc)
      call calc_predictors (rad, bccoef, mask_rs=mask_bc)
      call apply_bias_corr (rad, bccoef, mask_rs=mask_bc)
      !---------------------------
      ! feed back changes to 3dvar
      !---------------------------
      call restore_bcor    (rad, obs, mask_rs=mask_bc)
      call destruct        (rad)
    else
      call finish ("radbias_bfg","invalid value for 'biascor_mode'")
    endif

  end subroutine radbias_bfg

!------------------------------------------------------------------------------

  subroutine restore_bcor (rad, obs, cloud, mask_rs)
    !------------------------------------------------------------------
    ! restore bias correction and bias corrected brightness temperature
    ! optionally restore cloud flag
    ! in 3dvar derived type
    !------------------------------------------------------------------
    type (t_radv)     ,intent(inout)        :: rad(:)  ! radiance data
    type (t_obs_set)  ,intent(inout)        :: obs     ! observation data (3dvar type)
    logical           ,intent(in), optional :: cloud   ! restore cloud flag
    logical,           intent(in), optional :: mask_rs(:)
    integer :: i
    integer :: j, l, k, ib, is, ip0, ipn
    logical :: lc, lp

    lc = .false.; if (present(cloud)) lc = cloud

    do     i  = 1, size (rad)
      if (present(mask_rs)) then
        if (.not.mask_rs(i)) CYCLE
      end if
      do   l  = 1, rad(i)% n_rec
        ib  = rad(i)% i_box   (l)
        is  = rad(i)% i_reprt (l)
        lp  =      associated (obs% o(ib)% bcpred) &
             .and. associated (rad(i)%     pred  )
        ip0 = (obs% o(ib)% spot(is)% bcp% i)
        ipn = (obs% o(ib)% spot(is)% bcp% n)
        do j = 1, rad(i)% i% n_chan
          k  =    rad(i)% i_body (j,l)
          if (k>0) then
            !-------------------------------
            ! check for valid channel number
            !-------------------------------
            if (nint(obs% o(ib)% olev (k)) /= rad(i)% i% chan(j)) then
              print*,'restore_bcor invalid',i,l,j,k,ib,nint(obs% o(ib)% olev (k)),rad(i)% i% chan(j)
              print*,'file invalid',rad(i)% i%satid,rad(i)% i%grid,rad(i)% i%n_chan
              call finish('restore_bcor','invalid channel number')
            end if
            !---------------------------------------------------
            ! restore bias correction and brightness temperature
            !---------------------------------------------------
            if (rad(i)% bcor_(j,l) /= rinvalid) then
              rad(i)% bt_bcor(j,l) = rad(i)% bt_obs(j,l) - rad(i)% bcor_(j,l)
              obs% o(ib)% body(k)% o  =   rad(i)% bt_bcor(j,l)
              obs% o(ib)% body(k)% bc = - rad(i)% bcor_  (j,l)
            else
              rad(i)% bt_bcor(j,l) = rad(i)% bt_obs(j,l)
              obs% o(ib)% body(k)% o  = rad(i)% bt_obs(j,l)
              obs% o(ib)% body(k)% bc = 0.
              call decr_use (obs% o(ib)% body(k)% use, &
                             check = CHK_BIASCOR        )
            endif
            if (ldeb(obs% o(ib)% spot(is))) then
              write(usd,*) dpref,'bcor',obs% o(ib)% spot(is)%hd%id,&
                   rad(i)% i% chan(j),- rad(i)% bcor_(j,l)
            end if
            !------------------------------
            ! optionally restore predictors
            !------------------------------
            if (lp) then
              if (size(rad(i)% pred,1) /= ipn) then
                write(0,*) 'restore_bcor: set, size (rad%pred), ipn =',&
                            i, size(rad(i)% pred,1),ipn
                call finish ('restore_bcor','size (rad%pred) /= ipn')
              endif
              if (size(rad(i)% pred,3) > 1) &
                   call finish ('restore_bcor','P_THICK_TR* predictors not implmented for VarBC so far')
              obs% o(ib)% bcpred (ip0+1:ip0+ipn) = rad(i)% pred (:,l,1)
            endif
            !------------------------------
            ! optionally restore cloud flag
            !------------------------------
            if (lc) then
              select case (rad(i)% cloudy(j,l))
              case (1)
                call decr_use (obs% o(ib)% body(k)% use, check= CHK_CLOUD, &
                                                         lflag= .true.     )
              case (2)
                obs%   o(ib)% body(k)% use% flags = ibclr (  &
                  obs% o(ib)% body(k)% use% flags, CHK_CLOUD )
              case (3)
                if ( obs% o(ib)% body(k)% use% state >  STAT_PAS_REJ .and. &
                     obs% o(ib)% body(k)% use% state /= STAT_REJECTED) then
                  call decr_use (obs% o(ib)% body(k)% use, check= CHK_INSDAT,  &
                                                           state= STAT_REJECTED)
                end if
              end select
            endif
          endif
        end do
      end do
    end do

  end subroutine restore_bcor

!==============================================================================

  subroutine radbias_afg (obs, pass, y, x)
  type (t_obs_set) ,intent(inout) :: obs  ! observation
  integer          ,intent(in)    :: pass ! 1 for fg; 2 for analysis
  type (t_vector)  ,intent(in)    :: y    ! background, observation space
  type (t_vector)  ,intent(in)    :: x    ! background, interpolation space
  optional                        :: x
  !-----------------------------------------------------------------------
  ! Radiance bias correction routine to be called after first guess check
  ! or after analysis.
  ! Updates bias correction coefficients.
  ! Write updated correction coefficient file.
  ! Use updated correction in subsequent analysis.
  !-----------------------------------------------------------------------

    !----------------
    ! local variables
    !----------------
    integer                  :: iset     ! sat./instr. set index
    integer                  :: nset     ! number of sat./instr. sets
    integer                  :: iu       ! I/O unit number
    integer                  :: sat_wmo  ! wmo satellite  id
    integer                  :: ins_wmo  ! wmo instrument id

    type (t_stats_data) ,pointer     :: st       ! pointer to statistics data
    type (t_bcor_coef)  ,pointer     :: b        ! pointer to bias corr.coeff.
    type (t_stats_data) ,pointer     :: stat (:) ! bias statistics
    type (t_radv)       ,pointer     :: rad  (:) ! radiance data

    !---------------------------------
    ! check if update shall be applied
    !---------------------------------
    select case (pass)
    case (1)
      if (iand (biascor_mode, BC_UP_FG) /= BC_UP_FG) return
    case (2)
      if (iand (biascor_mode, BC_UP_AN) /= BC_UP_AN) return
    case default
      call finish("radbias_afg","invalid value for 'pass'")
    end select
    if (.not. associated (bccoef)) &
      call finish ('radbias_afg','bias correction not available.')

    !-------------------------------------
    ! allocate memory to gather statistics
    !-------------------------------------

    nset = size (bccoef)
    allocate (stat (nset))
    allocate (rad  (nset))
    do iset = 1, nset
      call construct (stat   (iset),            &
                      bccoef (iset)%    n_pred, &
                      bccoef (iset)% i% n_fov,  &
                      bccoef (iset)% i% n_chan  )
    end do

    !-----------------------------
    ! derive correction statistics
    !-----------------------------
    call fill_rad          (rad,  bccoef% i,  obs, x, y, lerr=.true., cp_opts=.true., lbc=.true.)
    call calc_predictors   (rad,  bccoef)
    call calc_bcor_t_stats (stat, bccoef, rad)

    !-------------------------------------------
    ! loop over sets of satellites / instruments
    !-------------------------------------------
    do iset = 1, nset
      st => stat   (iset)
      b  => bccoef (iset)
      !----------------------------
      ! sum over processor elements
      !----------------------------
      st% btdn       = p_sum (st% btdn)
      st% btdscan    = p_sum (st% btdscan)
      st% d2         = p_sum (st% d2)
      st% pred_i     = p_sum (st% pred_i)
      st% btd_pred_i = p_sum (st% btd_pred_i)
      st% pred_ij    = p_sum (st% pred_ij)
      st% n          =   sum (st% btdn, dim=1)
      !--------------------------
      ! merge with old statistics
      !--------------------------
!     b% exp = nex
      call merge_bcor_stats (b, st,  &
           cyyyymmddhhmm (ana_time), &
           cyyyymmddhhmm (run_time), &
           days=0._wp                )
      !------------------------------------
      ! update bias correction coefficients
      !------------------------------------
      call calc_bcor_coefs (b% sta, b)
    end do

    !--------------------------------------------
    ! write bias coefficient and statistics files
    !--------------------------------------------
    if (dace% lpio) then
      iu = get_unit_number()
      call   write_bcor_file (bccoef, 6,  meta=.true.)
      select case (format)
      case (0, 1)
         if (format==1) &
              coeff_ofile = path_file (output,                       &
              path_file (bc_path, bc_filename (OT_RAD)))
         call write_bcor_file (bccoef, iu, file=coeff_ofile)
      case (2, 3)
         do iset = 1, nset
            sat_wmo =              bccoef(iset)% i% satid
            ins_wmo = instr_rttov (bccoef(iset)% i% instr(1), sat_wmo)
            if (format==2) coeff_ofile = path_file (output,            &
                 path_file (bc_path,           &
                 bc_filename (OT_RAD, sat_wmo)))
            if (format==3) coeff_ofile = path_file (output,            &
                 path_file (bc_path,                  &
                 bc_filename(OT_RAD, sat_wmo,ins_wmo)))
            call write_bcor_file (bccoef(iset:iset), iu, file=coeff_ofile)
         end do
      case default
         call finish("radbias_afg","invalid value for 'format'")
      end select
      call return_unit_number (iu)
    endif

    !--------------------------------------------------
    ! apply bias correction to be used in this analysis
    !--------------------------------------------------
    if (update_ana                         .and. &
        pass == 1                          .and. &
        iand (biascor_mode, BC_3DVAR) == BC_3DVAR) then
      call apply_bias_corr (rad, bccoef)
      call restore_bcor    (rad, obs)
    endif

    !---------
    ! clean up
    !---------
    call destruct (stat)
    deallocate    (stat)
    call destruct (rad)
    deallocate    (rad)

  end subroutine radbias_afg
!------------------------------------------------------------------------------
  subroutine radbias_init (obs)
  type (t_obs) ,intent(in) :: obs  ! observations (global)
  !------------------------------------------------------------------
  ! Initialize this module:
  ! read namelist /rad_biascor/
  ! read bias-correction coefficient files
  ! optionally: read bias-correction templates from namelist
  !             extend bias-correction coefficient list by templates
  !-----------------------------------------------------------------
    integer                    :: ierr  ! error return parameter
    integer                    :: iu    ! I/O unit number
    integer                    :: i, j  ! bias-correction set index
    integer                    :: l, k  ! loop indices
    integer                    :: npr   ! number of special predictors
    integer                    :: ic, icb, ii ! channel/instrument indices
    integer                    :: ntmpl ! number of templates read
    integer                    :: is    ! spot index
    type(t_spot)      ,pointer :: s     ! pointer to report/fov data
    integer                    :: satid ! satellite id
    integer                    :: grid  ! grid
    type(t_bcor_coef) ,pointer :: tmp(:)! temporary bias correction set
    type(t_bcor_coef) ,pointer :: tpl(:)! templates read from namelist
    type(t_bcor_coef) ,pointer :: bc
    logical           ,pointer :: tpu(:)! used template flag
    type(t_rad_set),   pointer :: rs
    real(wp)                   :: x

    character(len=120)         :: message

    !===========================
    ! read namelist /RADBIASCOR/
    !===========================
    !--------------------------------
    ! set namelist variables defaults
    !--------------------------------
    update_ana     = .true.
    select case (flag_biasc_rad)
    case (-1)
      !------------------------------------------
      ! use old 1dvar static bias correction file
      !------------------------------------------
      coeff_ifile  = 'bias.out' ! coefficient input  file name
      coeff_ofile  = 'bias.new' ! coefficient output file name
      biascor_mode = 0          ! derive bias correction in 3dvar
      format       = 0          ! naming convention for file write
    case (0:1)
      !--------------------------------------------------------
      ! online bias correction, use cycled bias correction file
      !--------------------------------------------------------
      coeff_ifile  = ''
      coeff_ofile  = ''
      biascor_mode = BC_3DVAR + BC_UP_FG
      format       = 1
    end select
    !---------------------------
    ! read namelist /RADBIASCOR/
    !---------------------------
    if (dace% lpio) then
      call position_nml ('RADBIASCOR', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=RADBIASCOR ,iostat=ierr)
        if (ierr/=0) call finish ('nml_run_flags',&
                                  'ERROR in namelist /RADBIASCOR/')
#else
        read (nnml ,nml=RADBIASCOR)
#endif
      end select
      !------------------
      ! adjust path names
      !------------------
      coeff_ifile = path_file (input,  coeff_ifile)
      coeff_ofile = path_file (output, coeff_ofile)
      !---------------------------
      ! adjust biascor_mod for VBC
      !---------------------------
      if (iand (BC_VARBC, biascor_mode) /= 0) then
        biascor_mode = BC_3DVAR + BC_UP_AN + BC_VARBC
        update_ana   = .false.
      endif
      !---------------------------
      ! printout adjusted namelist
      !---------------------------
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      write(6,'(a)')  '  Namelist /RADBIASCOR/ (adjusted settings)'
      write(6,'( )')
      write(6,'(2a)')         '    coeff_ifile    = ',trim(coeff_ifile)
      write(6,'(2a)')         '    coeff_ofile    = ',trim(coeff_ofile)
      write(6,'(a,i2)')       '    biascor_mode   = ',biascor_mode
      write(6,'(a,l1)')       '    update_ana     = ',update_ana
      write(6,'(a,i2)')       '    format         = ',format
      write(6,'(a,f6.1)')     '    min_entries    = ',min_entries
      write(6,'(a,l1)')       '    no_pred        = ',no_pred
      write(6,'(a,l1)')       '    no_bias        = ',no_bias
      write(6,'(a,l1)')       '    precon         = ',precon
      write(6,'(a,i1)')       '    atm_thick_vers = ',atm_thick_vers
      write(6,'(a,f7.3)')     '    offset_max_age = ',offset_max_age
      write(6,'(a,f7.3)')     '    offset_t_decay = ',offset_t_decay
      write(6,'( )')
    endif
    !-------------------
    ! broadcast namelist
    !-------------------
    call p_bcast (coeff_ifile,   dace% pio)
    call p_bcast (coeff_ofile,   dace% pio)
    call p_bcast (biascor_mode,  dace% pio)
    call p_bcast (update_ana,    dace% pio)
    call p_bcast (format,        dace% pio)
    call p_bcast (min_entries,   dace% pio)
    call p_bcast (offset_max_age,dace% pio)
    call p_bcast (offset_t_decay,dace% pio)
    call p_bcast (no_pred,       dace% pio)
    call p_bcast (no_bias,       dace% pio)
    call p_bcast (precon,        dace% pio)
    call p_bcast (atm_thick_vers,dace% pio)
    !=======================================
    ! read bias-correction coefficient files
    ! (currently on all PEs)
    !=======================================
    if (iand (biascor_mode, BC_3DVAR+BC_UP_FG+BC_UP_AN+BC_VARBC) /= 0) then
      if (dace% lpio) then
        write(6,'(a)') repeat('-',79)
        write(6,'( )')
        write(6,'(a)')  '  Read bias correction coefficient files: '//trim(coeff_ifile)
        write(6,'( )')
      endif
      iu = get_unit_number()

      if (format == 0) then
        !------------------------
        ! old obsolete convention
        !------------------------
        nullify             (bccoef)
        call read_bcor_file (bccoef, iu, ierr, coeff_ifile)
        if (ierr /= 0) &
          call finish('radbias_init','error while reading '//trim(coeff_ifile))
      else
        !---------------
        ! new convention
        !---------------
        allocate (bccoef(0))
        ierr = 1
        do i = 1, nbcf
          if (bc_obstyp(i) == OT_RAD) then
            if (dace% lpio) write(6,'(a)')  '      '//trim(bc_paths(i))
            call read_bcor_file (bccoef, iu, ierr, bc_paths(i))
            if (ierr /= 0) call finish('radbias_init',                        &
                                     'error while reading '//trim(bc_paths(i)))
          endif
        end do
        if (ierr /= 0 .and. .not.fallback) &
          call finish('radbias_init','no bias correction read')
      endif
      bccoef% exp = nex
      !--------------------
      ! diagnostic printout
      !--------------------
      if (dace% lpio) then
        call write_bcor_file (bccoef, 6,  meta=.true.)
        call write_bcor_file (bccoef, iu, meta=.true.,            &
                              file=path_file (output, 'bias.tmpl'))
      endif
      call return_unit_number (iu)
      !----------------------------------------------
      ! fallback option: read templates from namelist
      !----------------------------------------------
      allocate (tpl(0))
      if (fallback) then
        if (dace% lpio) then
          write(6,'(a)') repeat('-',79)
          write(6,'( )')
          write(6,'(a)')  '  Read bias correction coefficient templates:'
          write(6,'( )')
          rewind (nnml)
        endif
        call read_bcor_file (tpl, nnml, ierr)
        do i = 1, size(tpl)
          if (tpl(i)% sti% alloc) call destruct (tpl(i)% sti)
          if (tpl(i)% sta% alloc) call destruct (tpl(i)% sta)
        end do
      endif
      ntmpl = size(tpl)
      allocate (tpu (ntmpl))
      tpu = .false.
      !-----------------------------
      ! check for valid coefficients
      !-----------------------------

      do is=1,obs% n_spot                    ! loop over reports/fovs
        s => obs% spot(is)                   !
        if (s% hd% obstype /= OT_RAD) cycle  ! only treat radiances
        !--------------------------------------------
        ! choose combination of satellite/target grid
        !--------------------------------------------
        satid = s% hd% satid                  ! satellite id
        grid  = s% hd% grid_id                ! grid id
        i     = set_indx (bccoef% i, satid=satid, grid=grid)
        if (i < 1) then
          if (.not.fallback) then
            write(message,'("no suitable bias correction set for satellite ",&
                 &I3," and grid ",I3)') satid, grid
            write(0,'()')
            write(0,'(a)') message
            do i = 1, size (bccoef)
              if (bccoef(i)% i% satid < 0) cycle
              write(0,'("                                     set ",i2," : ",&
                 &I3," and grid ",I3)') i, bccoef(i)% i% satid, bccoef(i)% i% grid
            end do
            call finish('radbias_init', trim(message))
          else
            !------------------------------------------
            ! select templates for missing coefficients
            !------------------------------------------
            i = set_indx (tpl% i, satid=satid, grid=grid)
            if (i < 1) then
              write(message,'("no suitable bias correction template for &
                 &satellite ",I3," and grid ",I3)') satid, grid
              call finish('radbias_init', trim(message))
            else
              tpu (i) = .true.
            endif
          endif
        endif
      end do
      !----------------------------------------
      ! Prepare transmission factor calculation
      !----------------------------------------
      do j = 1, size(bccoef)
        bc => bccoef(j)
        npr = 0
        do k = 1, size(chan_dep_pred)
          npr = npr + count(bc%pred(1:bc%n_pred) == chan_dep_pred(k))
        end do
        if (npr > 0) then
          bc%chan_dep_pred = .true.
          i = set_indx(rad_set(1:n_set), satid=bc%i%satid, grid=bc%i%grid)
          if (i > 0) then
            rs => rad_set(i)
            allocate(rs%bc%p_tr(npr,2), rs%bc%ib_ch(rs%n_chan), rs%bc%type(npr))
            rs%bc%p_tr  = -1._wp
            rs%bc%ib_ch = 0
            rs%bc%type  = 0
            allocate(bc%i%bc%p_tr(npr,2), bc%i%bc%ib_ch(bc%i%n_chan), bc%i%bc%type(npr))
            bc%i%bc%p_tr  = -1._wp
            bc%i%bc%ib_ch = 0
            bc%i%bc%type  = 0
            l = 0
            do k = 1, bc%n_pred
              if (any(bc%pred(k) == chan_dep_pred)) then
                l = l + 1
                rs%bc%p_tr(l,1) = bc% pred_p1(k)
                rs%bc%p_tr(l,2) = bc% pred_p2(k)
                if (rs%bc%p_tr(l,2) > rs%bc%p_tr(l,1)) then
                  x               = rs%bc%p_tr(l,2)
                  rs%bc%p_tr(l,2) = rs%bc%p_tr(l,1)
                  rs%bc%p_tr(l,1) = x
                end if
                select case(bc%pred(k))
                case(P_THICK_TR)
                  rs%bc%type(l) = 1
                case(P_THICK_TRF)
                  rs%bc%type(l) = 2
                case(P_TRANS_T)
                  rs%bc%type(l) = 3
                case(P_TRANS_TS)
                  rs%bc%type(l) = 4
                case(P_TRANS_P)
                  rs%bc%type(l) = 5
                case(P_TRANS_DP)
                  rs%bc%type(l) = 6
                end select
                ii = 1
                do ic = 1, rs%n_chan
                  if (ic > rs%o_ch_i(ii) + rs%n_ch_i(ii)) ii = ii+1
                  icb = chan_indx(rs%instr(ii), rs%chan(ic), bc%i)
                  if (icb > 0) then
                    if (bc%pred_use(k,icb) >= 0) &
                         rs%bc%ib_ch  (ic ) = ibset(rs%bc%ib_ch  (ic ), l)
                  else
                    write(0,*) 'satid=',bc%i%satid,'grid=',bc%i%grid,'ii=',ii,&
                         'instr=',rs%instr(ii), 'chan=',rs%chan(ic)
                    call finish('radbias_init','no corresponding channel')
                  end if
                end do
                do icb = 1, bc%i%n_chan
                  if (bc%pred_use(k,icb) >= 0) &
                       bc%i%bc%ib_ch(icb) = ibset(bc%i%bc%ib_ch(icb), l)
                end do
              end if
            end do
            rs%bc%n_tr   = l
            bc%i%bc%n_tr = l
            bc%i%bc%p_tr = rs%bc%p_tr
            bc%i%bc%type = rs%bc%type
          end if
        end if
      end do

      !-------------------
      ! rescale statistics
      !-------------------
      do j = 1, size(bccoef)
        call scale_bcor_stats (bccoef(j),                            &
          cyyyymmddhhmm (ana_time),                                  &
          cyyyymmddhhmm (run_time),                                  &
          days (ana_time - time_cyyyymmddhhmm (bccoef(j)% last_date)))
        if (no_pred) then
          where (bccoef(j)% pred_use == 1) bccoef(j)% pred_use = -1
        endif
        if (no_bias .and. bccoef(j)% sta% alloc) then
          bccoef(j)% sta% btdscan    = 0._wp
          bccoef(j)% sta% d2         = 0._wp
          bccoef(j)% sta% btd_pred_i = 0._wp
        endif
      end do
      !-------------------------
      ! extend list by templates
      !-------------------------
      tpu = p_or (tpu)
      i = count (tpu)
      if (i>0) then
        j = size (bccoef)
        tmp => bccoef
        allocate (bccoef (j+i))
        bccoef (1:j) = tmp
        deallocate (tmp)
        do i = 1, ntmpl
          if (tpu(i)) then
            j = j + 1
            bccoef (j) = tpl (i)
          endif
        end do
      endif
      deallocate (tpl)
      deallocate (tpu)
    else
      allocate (bccoef(0))
    endif

  end subroutine radbias_init

!==============================================================================
  subroutine read_fdbk (rad, file, set, pfile, lana)
  !----------------------------------------------------
  ! read 3dvar feedback file, store data in type t_radv
  !----------------------------------------------------
  type (t_radv)    ,intent(out) :: rad(:) ! derived type to fill
  character(len=*) ,intent(in)  :: file   ! 3dvar feedback file name
  type (t_rad_set) ,intent(in)  :: set(:) ! sets of instruments/channels
  character(len=*) ,intent(in)  :: pfile  ! 3dvar profiles file name
  logical          ,intent(in)  :: lana   ! read analysis (default: firstguess)
  optional                      :: pfile, lana
  target                        :: rad

    integer  ,allocatable :: i_body  (:)  ! Index of the record begining
    integer  ,allocatable :: l_body  (:)  ! Length of the record
    integer  ,allocatable :: obstype (:)  ! Observation type
    integer  ,allocatable :: ident   (:)  ! Satellite IDs
    integer  ,allocatable :: grid    (:)  ! grid (instrument ID)
    real(wp) ,allocatable :: lon     (:)  ! Longitude (degree)
    real(wp) ,allocatable :: lat     (:)  ! Latitude  (degree)
    integer  ,allocatable :: mdlsfc  (:)  ! model surface flag
    integer  ,allocatable :: phase   (:)  ! FOVs
    integer  ,allocatable :: state   (:)  ! 3dVar-state of the data
    integer  ,allocatable :: flags   (:)  ! 3dVar flags
    real(wp) ,allocatable :: stzen   (:)  ! Satellite zenith angle
    real(wp) ,allocatable :: bc_obs  (:)  ! Bias corrected observations
    real(wp) ,allocatable :: bcor    (:)  ! Bias
    real(wp) ,allocatable :: fg      (:)  ! First guess
    real(wp) ,allocatable :: level   (:)  ! Channel# + f_ins * instrument#
    integer  ,allocatable :: chan    (:)  ! Channel number
    integer  ,allocatable :: instr   (:)  ! Instrument
    integer               :: idx (1)      ! indices (from get_veri_index)
    integer               :: nidx         ! number of indices found
    integer               :: i,is         ! indices
    integer               :: j, k, l      ! indices
    integer  ,allocatable :: i_set  (:)   ! index to set
    integer  ,allocatable :: i_prof (:)   ! profile index
    integer  ,allocatable :: i_reor (:)   ! reordering for profile files
    integer  ,allocatable :: id_mon (:)   ! observation id in mon-file
    integer  ,allocatable :: id_prof(:)   ! observation id in prof-file
    real(wp) ,allocatable :: surf   (:)   ! surface variable
    real(wp) ,allocatable :: prof (:,:)   ! profile variable
    integer               :: n_set (size(set)) ! number of entries in set
    type(t_radv) ,pointer :: r            ! pointer into rad(:)
    type(t_fdbk)          :: fdbk         ! feedback file derived type
    logical               :: lprof        ! profile file present
    integer               :: n_prof       ! number of profiles   in prof-file
    integer               :: n_lev        ! number of levels     in prof-file
    integer               :: n_param      ! number of parameters in prof-file
    logical               :: la           ! read analysis (default: firstguess)

    !-----------------------------
    ! check for optional arguments
    !-----------------------------
    lprof = .false. ;if (present(pfile)) lprof = (pfile /= '')
    la    = .false. ;if (present(lana))  la    = lana

    !-------------------------------------------
    ! Open 3DVAR monitoring file, read meta data
    !-------------------------------------------
    call setup_fdbk(fdbk% nc)
    call open_fdbk_read(fdbk, file)
    if(fdbk% nc% error /= NF90_NOERR) then
       write(6,*) 'Error. Could not open file: ', file
       call finish('read_fdbk','cannot open: '//trim(file))
    end if
    call read_meta(fdbk)

    !--------------------------
    ! allocate temporary arrays
    !--------------------------
    allocate (i_body (fdbk% n_hdr ))
    allocate (l_body (fdbk% n_hdr ))
    allocate (obstype(fdbk% n_hdr ))
    allocate (ident  (fdbk% n_hdr ))
    allocate (grid   (fdbk% n_hdr ))
    allocate (phase  (fdbk% n_hdr ))
    allocate (stzen  (fdbk% n_hdr ))
    allocate (i_set  (fdbk% n_hdr ))
    allocate (i_prof (fdbk% n_hdr ))
    allocate (i_reor (fdbk% n_hdr ))
    allocate (id_mon (fdbk% n_hdr ))
    allocate (mdlsfc (fdbk% n_hdr ))
    allocate (lon    (fdbk% n_hdr ))
    allocate (lat    (fdbk% n_hdr ))
    allocate (state  (fdbk% n_body))
    allocate (flags  (fdbk% n_body))
    allocate (bc_obs (fdbk% n_body))
    allocate (bcor   (fdbk% n_body))
    allocate (level  (fdbk% n_body))
    allocate (chan   (fdbk% n_body))
    allocate (instr  (fdbk% n_body))
    allocate (fg     (fdbk% n_body))

    !---------------------
    ! read monitoring data
    !---------------------
    ncid  = fdbk%nc% ncid
    stanc = 1
    counc = fdbk% n_hdr
    strnc = 1
    call get_var (i_body,  'i_body')
    call get_var (l_body,  'l_body')
    call get_var (obstype, 'obstype')
    call get_var (ident,   'ident')
    call get_var (grid,    'instype') ; grid(:) = rttov_instr(grid(:), ident(:))
    call get_var (phase,   'phase')
    call get_var (stzen,   'sat_zenit')
    call get_var (lon,     'lon')
    call get_var (lat,     'lat')
    call get_var (mdlsfc,  'mdlsfc')
    call get_var (id_mon,  'obs_id')
    stanc = 1
    counc = fdbk% n_body
    strnc = 1
    call get_var (state,   'state')
    call get_var (flags,   'flags')
    call get_var (bc_obs,  'obs' , fill = rinvalid)
    call get_var (bcor,    'bcor', fill = 0._wp)
    call get_var (level,   'level')
    call get_var (instr,   'level_sig')

    !-----------------------------
    ! read first guess or analysis
    !-----------------------------
    if (la) then
      call get_veri_index (idx, nidx, fdbk,          &
                           run_type   = VT_ANALYSIS, &
                           ens_member = VE_DETERM    )
      if (nidx /= 1) call finish('read_fdbk','analysis not found')
    else
      call get_veri_index (idx, nidx, fdbk,            &
                           run_type   = VT_FIRSTGUESS, &
                           ens_member = VE_DETERM      )
      if (nidx /= 1) call finish('read_fdbk','first guess not found')
    endif

    fg = get_veri (fdbk, idx(1))

    !---------------------
    ! conversions required
    !---------------------
    chan   = nint (level)             ! channel number
    call reverse_code (state, stats)  ! 3dvar status
    call reverse_bits (flags, chks)   ! 3dvar flags

    !------------------------
    ! count data, set indices
    !------------------------
    n_set = 0
    do i = 1, fdbk% n_hdr
      if (obstype(i) /= OT_RAD) cycle
      if (grid(i) < 0) then
        grid(i) = rttov_instr(instr (i_body (i)), ident(i)) !2nd try: instr.of 1st chan
        write(6,'("*** WARNING: no valid grid-id for header ",I8,", use grid=",I2.2)') &
             &i, grid(i)
      end if
      is = set_indx (set, satid=ident(i), grid=grid(i))
      if (is < 1) then                 ! 3rd try: mapping as in 1dvar
        select case (ident(i))
        case (784)                     ! AQUA satellite
          is = set_indx (set, satid=ident(i), grid=3)       ! AMSU-A mapping
        case default
          is = set_indx (set, satid=ident(i), grid=0)       ! HIRS mapping
        end select
        write(6,'("*** WARNING: no valid grid-id for header ",I8," (satid ",I3.3,&
             &", grid ",I2.2,"), use 1dvar-convention.")') i,ident(i),grid(i)
      endif
      if (is < 1) then
        write (6,'("read_fdbk: satid/grid ",I3.3,"/",I2.2," not in set:",50(1x,I3.3,"/",I2.2))') &
             ident(i),grid(i),(set(j)% satid, set(j)%grid, j=1,size(set))
        call finish ('read_fdbk','satid not in set')
      endif
      i_set  (i)  = is
      n_set  (is) = n_set (is) + 1
      i_prof (i)  = n_set (is)
    end do

    !-------------------------------------------
    ! set meta data in rad (derived type t_radv)
    !-------------------------------------------
    rad% model_date = 100 * fdbk% refdate + fdbk% reftime / 100  ! yyyymmddhh
    do is = 1, size (set)
      r             => rad   (is)
      r% i          =  set   (is)
      r% filename   =  file
      r% n_rec      =  n_set (is)
      !--------------------
      ! allocate components
      !--------------------
      allocate (r% fov     (             r% n_rec))
      allocate (r% stzen   (             r% n_rec))
      allocate (r% i_reprt (             r% n_rec))
      allocate (r% dlon    (             r% n_rec))
      allocate (r% dlat    (             r% n_rec))
      allocate (r% mdlsfc  (             r% n_rec))
      allocate (r% obsnum  (             r% n_rec))
      allocate (r% i_body  (r%i% n_chan, r% n_rec))
      allocate (r% bt_obs  (r%i% n_chan, r% n_rec))
      allocate (r% bt_bcor (r%i% n_chan, r% n_rec))
      allocate (r% bcor_   (r%i% n_chan, r% n_rec))
      allocate (r% bt_fg   (r%i% n_chan, r% n_rec))
      allocate (r% not_rej (r%i% n_chan, r% n_rec))
      allocate (r% cloudy  (r%i% n_chan, r% n_rec))
      allocate (r% state   (r%i% n_chan, r% n_rec))
      allocate (r% flags   (r%i% n_chan, r% n_rec))
      !---------------------------
      ! preset with invalid values
      !---------------------------
      r% bt_obs  = rinvalid
      r% bt_bcor = rinvalid
      r% bcor_   = rinvalid
      r% bt_fg   = rinvalid
      r% not_rej = .false.
      r% cloudy  = 0
      r% state   = STAT_DEFAULT
      r% flags   = 0
      r% i_reprt = 0
      r% i_body  = 0
      r% dlon    = rinvalid
      r% dlat    = rinvalid
      r% mdlsfc  = 0
    end do

    !--------------------
    ! store data in 'rad'
    !--------------------
    do i = 1, fdbk% n_hdr
      !-------
      ! header
      !-------
      if (obstype(i) /= OT_RAD) cycle
      is             =  i_set  (i)
      j              =  i_prof (i)
      r              => rad    (is)
      r% fov     (j) =  phase  (i)
      r% stzen   (j) =  stzen  (i)
      r% i_reprt (j) =          i
      r% dlon    (j) =  lon    (i)
      r% dlat    (j) =  lat    (i)
      r% mdlsfc  (j) =  mdlsfc (i)
      r% obsnum  (j) =  id_mon (i)
      !-----
      ! body
      !-----
      do k = i_body(i), i_body(i) + l_body(i) - 1
        l   = chan_indx (             instr(k),  chan(k), r% i)
        if (l<0) then   ! new convention: WMO sensor id
          l = chan_indx (rttov_instr (instr(k), ident(i)), chan(k), r% i)
        endif
        if (l<0) then
           call print_rad_set(r%i)
           write (6,*)'read_fdbk: satid ',r%i%satid,' instr ',rttov_instr(instr(k), ident(i)),&
                &instr(k),' channel ', chan(k),' not in set.'
           call finish ('read_fdbk','channel/instrument not in set')
        end if
        r% bt_fg   (l,j) =            fg     (k)
        r% bt_obs  (l,j) =            bc_obs (k) - bcor (k)
        r% bt_bcor (l,j) =            bc_obs (k)
        r% bcor_   (l,j) =          - bcor   (k)
        r% not_rej (l,j) = valid (int(state  (k), i1))
        r% state   (l,j) =            state  (k)
        r% flags   (l,j) =            flags  (k)
        r% i_body  (l,j) =                    k
        if (               btest (    flags  (k), CHK_CLOUD)) then
          r% cloudy  (l,j) = 1
        else
          r% cloudy  (l,j) = 2
        endif
      end do
    end do

    !----------------------------
    ! deallocate temporary arrays
    !----------------------------
    deallocate (i_body )
    deallocate (l_body )
    deallocate (obstype)
    deallocate (ident  )
    deallocate (grid   )
    deallocate (phase  )
    deallocate (stzen  )
    deallocate (lon    )
    deallocate (lat    )
    deallocate (mdlsfc )
    deallocate (state  )
    deallocate (flags  )
    deallocate (bc_obs )
    deallocate (bcor   )
    deallocate (level  )
    deallocate (chan   )
    deallocate (instr  )
    deallocate (fg     )

    !----------------------------
    ! Close 3DVAR monitoring file
    !----------------------------
    call close_fdbk   (fdbk)
    call cleanup_fdbk (fdbk)

    !=========================
    ! Optionally read profiles
    !=========================
    if (lprof) then
      !-----------------------------------------
      ! Open 3DVAR profiles file, read meta data
      !-----------------------------------------
      call chk (nf90_open (pfile, nf90_NOWRITE, ncid))
      call get_dim (n_prof,    "profile"  )
      call get_dim (n_lev,     "level"    )
      call get_dim (n_param,   "parameter")
      !----------------
      ! allocate arrays
      !----------------
      do is = 1, size (set)
        r        => rad   (is)
        r% n_lev =  n_lev
        allocate (r%  p    (n_lev, 1       )); r%   p    = rinvalid
        allocate (r%  t_fg (n_lev, r% n_rec)); r%   t_fg = rinvalid
        allocate (r%  q_fg (n_lev, r% n_rec)); r%   q_fg = rinvalid
        allocate (r% ps_fg (       r% n_rec)); r%  ps_fg = rinvalid
        allocate (r% ts_fg (       r% n_rec)); r%  ts_fg = rinvalid
        allocate (r% v10_abs_fg(   r% n_rec)); r% v10_abs_fg = rinvalid
      end do
      allocate (id_prof (       n_prof))
      allocate (surf    (       n_prof))
      allocate (prof    (n_lev, n_prof))
      !-------------------
      ! set up index array
      !-------------------
      counc = (/n_prof,1,1/)
      stanc = (/1,     1,1/)
      call get_var (id_prof,  'id')
      if (n_prof /= size (id_mon))                   &
        call finish ('read_fdbk','inconsistent files')
      j = 1
      do i = 1, size (id_mon)
        if (j > n_prof) j = 1
        if (id_prof (j) /= id_mon (i)) j = sum(minloc(abs(id_prof-id_mon(i))))
        if (id_prof (j) /= id_mon (i)) then
          write (0,*) 'missing id in profiles file:',id_mon (i)
          call finish ('read_fdbk','missing id in profiles file')
        endif
        i_reor (i) = j
        j = j + 1
      end do
      !----------------------------
      ! fill arrays: surface fields
      !----------------------------
      counc = (/1,n_prof,1/)
      stanc = (/1,     1,1/)
      call get_var (surf,  'ts')
      do i = 1, size (id_mon)
        rad(i_set(i))% ts_fg (i_prof(i)) = surf(i_reor(i))
      end do
      call get_var (surf,  'ps')
      do i = 1, size (id_mon)
        rad(i_set(i))% ps_fg (i_prof(i)) = surf(i_reor(i))
      end do
      call get_var (surf,  'v10m', fill = rinvalid)
      do i = 1, size (id_mon)
        rad(i_set(i))% v10_abs_fg(i_prof(i)) = surf(i_reor(i))
      end do
      !---------------------------------
      ! fill arrays: pressure coordinate
      !---------------------------------
      counc = (/n_lev,1,1/)
      stanc = (/    1,1,1/)
      call get_var (prof(:,1),  'level')
      do is = 1, size (set)
        rad(is)% n_lev = n_lev
        rad(is)% p(:,1)     = prof(:,1)
      end do
      !----------------------
      ! fill arrays: profiles
      !----------------------
      counc = (/n_lev,1,n_prof/)
      stanc = (/    1,1,     1/)
      call get_var (prof,  't')
      do i = 1, size (id_mon)
        rad(i_set(i))% t_fg(:,i_prof(i)) = prof(:,i_reor(i))
      end do
      call get_var (prof,  'q')
      do i = 1, size (id_mon)
        rad(i_set(i))% q_fg(:,i_prof(i)) = prof(:,i_reor(i))
      end do
      !--------------------------
      ! Close 3DVAR profiles file
      !--------------------------
      call chk (nf90_close (ncid))
      !----------------------------
      ! deallocate temporary arrays
      !----------------------------
      deallocate (id_prof)
      deallocate (surf   )
      deallocate (prof   )
    endif
    deallocate (i_set  )
    deallocate (i_prof )
    deallocate (i_reor )
    deallocate (id_mon )

  end subroutine read_fdbk

!------------------------------------------------------------------------------

  subroutine write_tbcor_fdbk (rad, r_tb, r_flags)
  !------------------------------------------------------------------
  ! restore bias correction and bias corrected brightness temperature
  ! in 3dvar-feedback file
  !------------------------------------------------------------------
  type(t_radv), intent(in) ,target :: rad (:)
  logical                ,optional :: r_tb    ! restore brightn.temp.,correction
  logical                ,optional :: r_flags ! restore (cloud) flags

    type(t_fdbk)              :: fdbk
    real(wp)     ,allocatable :: bcor  (:)
    real(wp)     ,allocatable :: obs   (:)
    integer      ,allocatable :: flags (:)
    integer                   :: i,j,k,l,m
    type(t_radv) ,pointer     :: r
    integer                   :: ncstat, vid_obs, vid_bcor, vid_flags
    real(wp)                  :: fillvalue     ! Invalid values
    logical                   :: l_tb, l_flags

    !----------------------------
    ! process optional parameters
    !----------------------------
    l_tb    = .true.  ;if (present(r_tb))    l_tb    = r_tb
    l_flags = .false. ;if (present(r_flags)) l_flags = r_flags
    !-------------------------------------------
    ! Open 3DVAR monitoring file, read meta data
    !-------------------------------------------
    if (size(rad)==0) return
    call open_fdbk_write (fdbk, rad(1)% filename)
    if (fdbk% nc% error /= NF90_NOERR) then
       write(6,*) 'Error. Could not open file: ', rad(1)% filename
       call finish('read_fdbk','cannot open: '//trim(rad(1)% filename))
    end if
    call read_meta(fdbk)

    !---------------------
    ! allocate temporaries
    !---------------------
    allocate (bcor  (fdbk% n_body))
    allocate (obs   (fdbk% n_body))
    allocate (flags (fdbk% n_body))
    !---------------------------------------------
    ! gather information from different satellites
    !---------------------------------------------
    do m = 1, size(rad)
      r => rad (m)
      do j = 1, r% n_rec
        i = r% i_reprt (j)
        if (i <= 0) cycle
        do l = 1, r% i% n_chan
          k = r% i_body (l,j)
          if (k > 0) then
            if (r% bcor_(l,j) == rinvalid) then
              bcor (k) = rinvalid
              if (l_flags) flags (k) = ibset (r% flags(l,j), CHK_BIASCOR)
            else
              bcor (k) = - r% bcor_   (l,j)
              if (l_flags) flags (k) = ibclr (r% flags(l,j), CHK_BIASCOR)
            endif
            obs    (k) =   r% bt_bcor (l,j)
          endif
        end do
      end do
    end do

    !------------------------------------
    ! re-write bias corrected observation
    !------------------------------------
    if (l_tb) then
      ncid = fdbk% nc% ncid
      ncstat = nf90_inq_varid (fdbk% nc% ncid, 'obs', vid_obs)
      if (ncstat /= NF90_NOERR) &
        call finish("write_tbcor_fdbk","inq_varid('obs') failed")
      call get_attr (fillvalue, 'obs', '_FillValue')
      where (obs == rinvalid) obs = fillvalue
      ncstat = nf90_put_var (fdbk% nc% ncid, vid_obs, obs)
      if (ncstat /= NF90_NOERR) &
        call finish("write_tbcor_fdbk","put_var('obs') failed")
      !-------------------------
      ! re-write bias correction
      !-------------------------
      ncstat = nf90_inq_varid (fdbk% nc% ncid, 'bcor', vid_bcor)
      if (ncstat /= NF90_NOERR) &
        call finish("write_tbcor_fdbk","inq_varid('bcor') failed")
      call get_attr (fillvalue, 'bcor', '_FillValue')
      where (bcor == rinvalid) bcor = fillvalue
      ncstat = nf90_put_var (fdbk% nc% ncid, vid_bcor, bcor)
    endif
    !-------------
    ! update flags
    !-------------
    if (l_flags) then
      call change_bits (flags, chks% code)   ! 3dvar flags -> feedback
      ncstat = nf90_inq_varid (fdbk% nc% ncid, 'flags', vid_flags)
      if (ncstat /= NF90_NOERR) &
        call finish("write_tbcor_fdbk","inq_varid('flags') failed")
      ncstat = nf90_put_var (fdbk% nc% ncid, vid_flags, flags)
    endif
    !---------------
    ! update history
    !---------------
    call add_history (fdbk, 'apply', 'bias', 'correction, flags')
    if (fdbk%nc%error /= NF90_NOERR) &
        call finish("write_tbcor_fdbk","add_history failed")
    !----------------------------
    ! Close 3DVAR monitoring file
    !----------------------------
    call close_fdbk   (fdbk)
    call cleanup_fdbk (fdbk)

  end subroutine write_tbcor_fdbk

!------------------------------------------------------------------------------
! Broadcast type t_bcor_coef
!---------------------------
#define VECTOR
#define DERIVED type(t_bcor_coef),dimension(:)
#define p_bcast_DERIVED p_bcast_bcor_coef
#undef  MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
end module mo_radbias_3dv
