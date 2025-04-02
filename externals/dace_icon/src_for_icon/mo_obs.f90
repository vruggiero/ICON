!
!+ provide generic operators on observation data type T_OBS.
!
MODULE mo_obs
!
! Description:
!   This module defines the generic operators (PRINT, SETUP, PROCESS) on the
!   observation data type T_OBS. Subroutine process_obs is the main entry
!   point. It calls the observation type specific routines from the
!   respective modules.
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
!  Changes for verification mode; cleanup
! V1_5         2009/05/25 Andreas Rhodin
!  Changed interfaces to subroutines
! V1_7         2009/08/24 Harald Anlauf
!  Changes for feedback file read; reactivate AMV from AVHRR (DK_AMV_NOAA)
! V1_9         2010/04/20 Andreas Rhodin
!  TSK_SHRINK in subroutines process: pass parameter 'state' to shrink_report
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  set_rules: add dbkz=10385 to list of SYNOPs
! V1_13        2011/11/01 Alexander Cress
!  changes for windprofilers, aircrafts, IASI
! V1_19        2012-04-16 Andreas Rhodin
!  nominal_height: moved from module mo_letkf; write diagnostics to feedback file
! V1_20        2012-06-18 Harald Anlauf
!  cleanup
! V1_22        2013-02-13 Andreas Rhodin
!  changes for wind profiler, STD and COSMO observation operators
!  nominal_height: pressure level for GPSRO impact parameter
!  implementation of vectorized K-mode (Robin Faulwetter)
! V1_23        2013-03-26 Andreas Rhodin
!  remove references to removed module mo_obs_odb
! V1_26        2013/06/27 Andreas Rhodin
!  subroutine nominal_height: check for validity of the Jacobi matrix
! V1_27        2013-11-08 Andreas Rhodin
!  dont write humidity weight (RAD) to pcc; fix bug for ICON+IFSstratosphere
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF
! V1_29        2014/04/02 Andreas Rhodin
!  changes for MODE-S and STD/ZTD observations
! V1_37        2014-12-23 Andreas Rhodin
!  t_obs_set: new component si; Add new BUFR codes to default rules tables
! V1_40        2015-02-27 Harald Anlauf
!  changes for new SHIP BUFR reports (A.Cress)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_43        2015-08-19 Andreas Rhodin
!  MEC: call observation module initialisation only if required
! V1_44        2015-09-30 Harald Anlauf
!  Preparations for Jason-2
! V1_45        2015-12-15 Andreas Rhodin
!  nominal_height: option to use actual background, use p instead of ln(p)
! V1_47        2016-06-06 Andreas Rhodin
!  disentangle /RULES/ default settings for ICON/COSMO/MEC
! V1_48        2016-10-06 Andreas Rhodin
!  changes for passive monitoring
! V1_50        2017-01-09 Andreas Rhodin
!  minor change
! V1_51        2017-02-24 Andreas Rhodin
!  process_obs: make parameter 'atm' optional
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM  2000  original code
! Andreas Rhodin  DWD    2001  changes for PSAS (common observation data type)
! Oliver Schmid   DWD    2005  new obs data type
! Gerhard Paul    DWD    2008  NetCDF input of observations
! Harald Anlauf   DWD    2008  minor changes and fixes
!==============================================================================
  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp, sp, i8        ! working precision kind parameter
  use mo_exception,  only: finish            ! abort routine
#ifndef NOBUFR
  use mo_run_params, only: p_readbufr        ! PE used to read BUFR data
#endif
  use mo_usstd,      only: p_h_usstd         ! US 1976 Standard Atmosphere
  !-----------------------------------------------
  ! observation data type and predefined constants
  !-----------------------------------------------
  use mo_dwd_tables,                        &! Datenbank-Kennziffer:
                     only: DK_SATOB_2,      &! Satellitenbeob, SATOB, Section 2
                           DK_SATOB_3,      &! Satellitenbeob, SATOB, Section 3
                           DK_AMV_EUMETSAT, &!
                           DK_AMV_GOES,     &!
                           DK_AMV_MODIS,    &!
                           DK_AMV_MODIS_OLD,&!
                           DK_AMV_MODIS_A,  &!
                           DK_AMV_MTV,      &!
                           DK_AMV_NOAA,     &!
                           DK_AMV_SENTINEL, &!
                           DK_AMDAR,        &!
                           DK_AIREP,        &!
                           DK_ACARS_USA,    &!
                           DK_ACARS_EU,     &!
                           DK_ACARS_CH,     &!
                           DK_ACARS,        &!
                           DK_ACARS_SINGLE, &!
                           DK_MODES_OLD,    &!
                           DK_MODES          !
  use mo_obs_set,    only: t_obs_set,       &!
                           t_obs_block,     &!
                           obs_block         !
  use mo_obs_tables, only: read_nml_report, &! read namelist /REPORT/
                           rept_char         ! report characteristics table
  use mo_t_obs,      only: t_obs,           &! observation data type
                           t_spot,          &!   component of t_obs
                           release_mem,     &! release unused memory
                           invalid,         &! invalid observation value
                           !-------------------------------------------
                           ! identifier for observation operator module
                           !-------------------------------------------
                           TEMP,            &! TEMP  soundings
                           SYNOP,           &! SYNOP observations
                           SCATT,           &! Scatterometer observations
                           TOVS,            &! TOVS
                           GPSRO,           &! GNSS Radio Occultations
                           GPSGB,           &! GNSS Ground Based
                           COSMO,           &! COSMO conventional operators
                           RADAR,           &! RADAR operator
                           AMV,             &! AMV satellite winds
                           AIREP,           &! Aircraft reports
                           SATEM,           &! IASI satem reports
                           WLIDAR,          &! Wind Lidar observations
                           SOIL,            &! ASCAT soil moisture operator
                           !-----------------------------------------
                           ! values of task argument to setup/process
                           !-----------------------------------------
                           TSK_INIT,TSK_READ,TSK_SET_CHR,TSK_SETUP_COLS,    &
                           TSK_SETUP_FULL, TSK_SETUP_FUL0, TSK_SHRINK,      &
                           TSK_R,TSK_K,TSK_Y,TSK_YH,TSK_H,TSK_SETUP_OP,     &
                           !-----------------------
                           ! observation type flags
                           !-----------------------
                           OBS_TV,                                          &
                           OBS_RH,                                          &
                           OBS_DUM,                                         &
                           !----------------------------
                           ! characteristics of operator
                           !----------------------------
                           CHR_NONL,                                        &
                           !--------------------
                           ! interpolation types
                           !--------------------
                           ITY_MCOLS,&
                           !-------------------------
                           ! byte sizes of data types
                           !-------------------------
                           obs_bytes, spot_bytes, body_bytes, int_bytes, &
                           wp_bytes,                                     &
                           set_byte_size, &! routine to setup sizes
                           !------------------------
                           ! namelist /OBSERVATIONS/
                           !------------------------
                           read_obs_nml,  &!
                           read_bufr,     &!
                           read_netcdf,   &!
                           !----------------------------
                           ! context of process_obs call
                           !----------------------------
                           po_context,    &!
                           po_ilns,       &!
                           po_lmon,       &!
                           POC_NN
  use mo_fdbk_tables,only: OT_GPSRO,      &! radio occultation  report id.
                           OT_GPSGB,      &! GNN ground based (STD/ZTD)
                           OT_RAD,        &! Radiances          report id.
                           VN_P            ! pressure 'varno'
  use mo_thinning,   only: read_nml_thin   ! read  namelist /THINNING/
  use mo_obs_rules,  only: new_rule,      &! set observation processing rules
                           t_set,         &! data type definition
                           PRC_NONE        ! special processing flag
  use mo_t_use,      only: STAT_DISMISS,  &! constant: don't use this parameter
                           STAT_ACTIVE_1   !          assimilate this parameter
  !---------------------
  ! matrices and vectors
  !---------------------
  use mo_dec_matrix, only: t_vector,      &! vector
                           t_vector_segm, &! vector segment
                           CSC,           &! compressed sparse column repr.
                           BDIAG,         &! block diagonal flag value
                           FULL,          &! full n x m representation flag
                           mp,            &! real kind for matrix coeffs.
                           operator (*),  &! matrix * vector
                           operator (-),  &! vector - vector
                           assignment(=), &! vector = vector
                           sub_block,     &! return submatrix
                           construct,     &!
                           destruct,      &!
                           allocate_block,&!
                           reallocate,    &!
                           delete_storage  !
  use mo_allocate,   only: enter_function,&!
                           leave_function  !
  !-----------------------------------
  ! observation error covariance model
  !-----------------------------------
  use mo_obs_err,      only: init_obs_err
  !-------------------------------------
  ! variational quality control defaults
  !-------------------------------------
  use mo_vqc,          only: svqc,         &! VQC threshold
                             vqc_form       ! formulation of VQC: shape of J_o
  !------------------------------------------------------------------
  ! data types for atmospheric (gridpoint) state and 'physical' space
  !------------------------------------------------------------------
  use mo_atm_state,    only: t_atm,            &! atmospheric state data type
                             p_sum_atm          ! sum over processors
  use mo_time,         only: t_time,           &! date & time derived type
                             operator(-),      &! time difference
                             iyyyymmdd,        &! conversion routines
                             ihhmm,            &!
                             seconds            !
  use mo_t_col,        only: t_cols             ! data type for model columns
  !-------------------------------
  ! specific observation operators
  !-------------------------------
  use mo_temp,         only: process_temp, print_temp
  use mo_occ_1d,       only: process_occ
  use mo_std,          only: process_std,      &!
              std_method =>  pl_method          ! method for plev estimation
  use mo_satem,        only: process_satem
  use mo_synop,        only: process_synop
  use mo_scatt,        only: process_scatt
  use mo_tovs,         only: process_tovs,     &! RADiance specific routine
                             rttov_mult_prof,  &! multiple profiles to RTTOV
                             process_tovs_mult,&! mult.p.  specific routine
                             read_fdbk_rad    ,&! post-process fdbk-read for TOVS
                             print_tovs_levs,  &! level info for TOVS
                             destruct_tovs_oe, &! destruct the t_range_fparse structures in rad_set
                             pl_method,        &! method for plev estimation
                             pl_e_bg_t,        &! temperature weight
                             pl_e_bg_rh,       &! rel.hum.    weight
                             pl_log,           &! use log(p) for mean, stdev
                             iatm_tovs_all
  use mo_rad,          only: rad_set, n_set     ! descriptions of radiance datasets
  use mo_amv,          only: process_amv        ! AMV        specific routine
  use mo_airep,        only: process_airep      ! Aircraft   specific routine
  use mo_cosmo_conv,   only: process_cosmo_conv ! COSMO      specific routine
  use mo_radar,        only: process_radar      ! RADAR      specific routine
  use mo_wlidar,       only: process_wlidar     ! Wind Lidar specific routine
  use mo_soil_obs,     only: process_soil       ! ASCAT soil specific routine
  use mo_occ,          only: set_occ,          &! set occ atmospheric data
                             destruct_occ,     &! deallocate occ atmosph. data
                             setup_occ_msis,   &! prepare transfer to MSIS
                             t_occ,            &! GPSRO specific table
                             t_ray,            &! GPSRO specific table
                             load_occ           ! load  GPSRO specific tables
  use mo_std,          only: set_std            ! initialize auxiliary data
  use mo_obs_bufr_dwd, only: read_obs_bufr
  use mo_obs_netcdf,   only: read_obs_netcdf
  use mo_obs_nml,      only: read_nml_obs
  use mo_fdbk_in,      only: read_feedback
  use ICO_grid,        only: gf                 ! occ atmospheric data
  !-----------------------
  ! communication routines
  !-----------------------
  use mo_mpi_dace,     only: dace,             &! MPI group info
                             p_sum,            &! parallel sum
                             p_max,            &! parallel max
                             p_ior,             &! parallel or
                             p_bcast,          &! generic MPI broadcast routine
                             p_bcast_ptr        ! generic MPI pointer bcast
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  !-----------------
  ! basic operations
  !-----------------
  public :: process_obs     ! generic observation processing routine
  public :: print_obs
  public :: nominal_height  ! Derive the nominal height (hPa) of observations.
  public :: set_time_slice  ! set meta data for temporal interpolation
  public :: PRT_COM, PRT_ALL, PRT_HEAD, PRT_TEMP
  public :: set_rules_icon  ! /rules/ defaults for global data assimilation
  public :: set_rules_cosmo ! /rules/ defaults for KENDA
  public :: set_dec_info    ! set dimension meta data of vectors and matrices
  !----------------------------------
  ! definition of print routine flags
  !----------------------------------
  integer, parameter :: PRT_COM  = 1 ! print components common to all observat.
  integer, parameter :: PRT_ALL  = 2 ! print all observations
  integer, parameter :: PRT_HEAD = 4 ! print header of all observations
  integer, parameter :: PRT_TEMP = 8 ! print TEMP observations
  !----------
  ! debugging
  !----------
  integer :: obs_debug = 0
!------------------------------------------------------------------------------
  !===========
  ! Interfaces
  !===========
  interface process_obs
    module procedure process_obs0 ! one      box
    module procedure process_obs1 ! multiple boxes
  end interface process_obs
!==============================================================================
contains
!==============================================================================
  subroutine print_summary (obs)
  type (t_obs) ,intent (in)           :: obs

    type (t_spot) ,pointer :: s
    integer                :: i

    character(len=*),parameter :: f1 = &
'(/"   id type bty sty statid   yyyymmddhhmm synp     lat     lon lev obs qcf"/)'
    character(len=*),parameter :: f2 = &
       '(i5,i5,2i4,1x,a8,i9.8,i4.4,i5.4,2f8.2,i4,i4,i4)'

    write (6,'(a)') repeat('-',79)
    write (6,f1)
    do i=1,obs% n_spot
      s => obs% spot(i)
      write (6,f2) s%id, s%hd%modtype, s%hd%buf_type, s%hd%buf_subtype, &
        s%statid, iyyyymmdd(s%actual_time),ihhmm(s%actual_time),        &
        ihhmm(s%hd%time), s%col%c%dlat, s%col%c%dlon, s%col% nlev, s%o%n
    end do
    write (6,'()')

  end subroutine print_summary
!==============================================================================
  subroutine print_obs (obs, flags, force)
  type (t_obs) ,intent (in)           :: obs
  integer      ,intent (in) ,optional :: flags
  logical      ,intent (in) ,optional :: force
    integer :: i
    integer :: fl
    logical :: fo

    fl = PRT_COM;    if (present(flags)) fl = flags
    fo = dace% lpio; if (present(force)) fo = force
    if (fo) then
      print *,'type (t_obs)'
      if(0/=fl) call print_head
      print *
      do i=1,obs% n_spot
        select case (obs% spot(i)% hd% modtype)
        case (TEMP)
          if(0/=iand(fl,PRT_ALL+PRT_TEMP+PRT_HEAD)) call print_spot
          if(0/=iand(fl,PRT_ALL+PRT_TEMP)) call print_temp (obs% spot(i), obs)
        case default
          if(0/=iand(fl,PRT_HEAD+PRT_ALL)) call print_spot
        end select
      end do
        if(0/=fl .and. fl/=PRT_COM) call print_head
        call print_mem
      print *,'end type (t_obs)'
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine print_head
      print *,'n_spot:',obs% n_spot
      print *,'n_obs :',obs% n_obs
      print *,'n_int :',obs% n_int
      print *,'n_par :',obs% n_par
      if(associated(obs%varno)) print *       ,&
        'varno',minval(obs% varno(:obs%n_obs)),&
          '...',maxval(obs% varno(:obs%n_obs))
      if(associated(obs%body )) print *         ,&
        'obs  ',minval(obs% body (:obs%n_obs)%o),&
          '...',maxval(obs% body (:obs%n_obs)%o)
      if(associated(obs%bger )) print *       ,&
        'bger ',minval(obs% bger (:obs%n_obs)),&
          '...',maxval(obs% bger (:obs%n_obs))
      if(associated(obs%t_int)) print *       ,&
        't_int',minval(obs% t_int(:obs%n_int)),&
          '...',maxval(obs% t_int(:obs%n_int))
      if(associated(obs%lev  )) print *       ,&
        'lev  ',minval(obs% lev  (:obs%n_int)),&
          '...',maxval(obs% lev  (:obs%n_int))
      if(associated(obs%bgeri)) print *       ,&
        'bgeri',minval(obs% bgeri(:obs%n_int)),&
          '...',maxval(obs% bgeri(:obs%n_int))
      if(associated(obs%bgi  )) print *       ,&
        'bgi  ',minval(obs% bgi  (:obs%n_int)),&
          '...',maxval(obs% bgi  (:obs%n_int))
      if(associated(obs%par  )) print *       ,&
        'par  ',minval(obs% par  (:obs%n_par)),&
          '...',maxval(obs% par  (:obs%n_par))
    end subroutine print_head
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine print_spot
      print *,'type (t_spot)'
      print *,'      i,id,modtype',i,obs% spot(i)% id,obs% spot(i)% hd% modtype
      print *,'      i',obs% spot(i)% o%i,&
                        obs% spot(i)% i%i,obs% spot(i)% p%i
      print *,'      n',obs% spot(i)% o%n,&
                        obs% spot(i)% i%n,obs% spot(i)% p%n
      if(associated(obs%body)) print*,'obs'                                  ,&
        minval(obs%body(obs%spot(i)%o%i+1:obs%spot(i)%o%i+obs%spot(i)%o%n)%o),&
        '...',                                                                &
        maxval(obs%body(obs%spot(i)%o%i+1:obs%spot(i)%o%i+obs%spot(i)%o%n)%o)
    end subroutine print_spot
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine print_mem
      integer :: nbytes
      character(len=*),parameter:: head=&
        '("component           :     words     bytes")'
      character(len=*),parameter:: form='(a20,":",2i10)'
      print *,'memory usage:'
      call set_byte_size
      nbytes = 0
      write (*,head)
      write (*,form) 'common', 1, obs_bytes
      nbytes = nbytes + obs_bytes
      write (*,form) 'metainfo', obs% n_spot, obs% n_spot*spot_bytes
      nbytes = nbytes + obs% n_spot*spot_bytes
      write (*,form) 'body', obs% n_obs, obs% n_obs*body_bytes
      nbytes = nbytes + obs% n_obs*body_bytes
      if (associated (obs%varno)) then
        write (*,form) 'observation types', obs% n_obs, obs% n_obs*int_bytes
        nbytes = nbytes + obs% n_obs*int_bytes
      endif
      if (associated (obs%varno)) then
        write (*,form) 'observations', obs% n_obs, obs% n_obs*wp_bytes
        nbytes = nbytes + obs% n_obs*wp_bytes
      endif
      if (associated (obs%t_int)) then
        write (*,form) 'psas types', obs% n_int, obs% n_int*int_bytes
        nbytes = nbytes + obs% n_int*int_bytes
      endif
      if (associated (obs% lev)) then
        write (*,form) 'psas levels', obs% n_int, obs% n_int*wp_bytes
        nbytes = nbytes + obs% n_int*wp_bytes
      endif
      write (*,form) 'total', 0, nbytes
    end subroutine print_mem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine print_obs
!==============================================================================
  subroutine process_obs1 (task, obs, atm, cols, xi, y, Jo, Jo_atm, local, &
                           state, context, ilns, lmon)
  !---------------------------------------------------------
  ! generic observation processing routine
  ! to be called for the complete set of observation 'boxes'
  !---------------------------------------------------------
  integer       ,intent(in)              :: task   ! what to do
  type(t_obs_set),intent(inout) ,target  :: obs    ! observations
  type(t_atm)   ,intent(in)    ,optional :: atm    ! atmospheric state
  type(t_cols)  ,intent(in)    ,optional :: cols(:)! model columns
  type(t_vector),intent(in)    ,optional :: xi     ! background (interp.space)
  type(t_vector),intent(inout) ,optional :: y      ! result     (observ.space)
  real(wp)      ,intent(inout) ,optional :: Jo     ! obs. cost funct. Jo
  type(t_atm)   ,intent(inout) ,optional :: Jo_atm ! gradient:d Jo/d atm
  logical       ,intent(in)    ,optional :: local  ! only process boxes locally
  integer       ,intent(in)    ,optional :: state  ! status flag
  integer       ,intent(in)    ,optional :: context! context of this call
  integer       ,intent(in)    ,optional :: ilns   ! call within in linesearch
  logical       ,intent(in)    ,optional :: lmon   ! monitoring call?

    logical                       :: llocal
    integer                       :: ib, io
    integer                       :: i
    type (t_obs_block)            :: ob
    type (t_spot)        ,pointer :: o
    integer                       :: n
    integer                       :: n_obs
    integer                       :: n_spot
    integer                       :: tsk
    type (t_vector_segm) ,pointer :: ys

    if (task==0) return

!   call enter_function

    !-------------------------------------------------------
    ! on subroutine entry: process some (optional) arguments
    !-------------------------------------------------------
    llocal =  .false.; if(present(local)) llocal   = local
    tsk    = task
    ! context of this call
    if (present(context)) then
      po_context = context
    else
      po_context = POC_NN
    end if
    if (present(ilns)) then
      po_ilns = ilns
    else
      po_ilns = -99
    end if
    if (present(lmon)) then
      po_lmon = lmon
    else
      po_lmon = .false.
    end if

    !----------------------------------------------------------------------
    ! for some specific tasks perform some initialisation of derived types:
    ! TSK_SETUP_FULL: oi,ii,di: description of vector spaces
    !----------------------------------------------------------------------
    if (iand (TSK_SETUP_FULL ,tsk) /= 0) then
      if (.not.associated(obs% oi)) then
        n_obs  = sum (obs% o% n_obs)
        n_spot = sum (obs% o% n_spot)
        allocate       (obs% si)
        allocate       (obs% oi)
        allocate       (obs% ii)
        allocate       (obs% di)
        call construct (obs% si, n_spot, size(obs% o))
        call construct (obs% oi, n_obs,  size(obs% o))
        call construct (obs% ii, 0,      size(obs% o))
        call construct (obs% di, 0,      size(obs% o))
        if (.not.associated(obs% R% b)) then
          call construct (obs% R, obs% oi, 'R', BDIAG)
        endif
      endif
      !----------------------------------------------------------
      ! set some quantities for vertical extension with IFS state
      !----------------------------------------------------------
      if (.not.present (atm)) call finish ('process_obs1','atm not present')
      !-------------------------------------------
      ! Derive pressure level for transfer to MSIS
      !-------------------------------------------
      call setup_occ_msis   (atm)
    endif
    if (iand(TSK_SETUP_COLS, tsk) /= 0) then
      iatm_tovs_all = 0
    end if
    !----------------------------------------------------------
    ! TSK_K: initialisation for linearized observation operator
    !----------------------------------------------------------
    if (iand (TSK_K ,tsk) /= 0) then
      if (.not.associated(obs% l% H% b)) then
        call construct   (obs% l% H ,obs% oi ,'H' ,BDIAG ,ci=obs% ii)
        call construct   (obs% l% x ,obs% ii ,'xi')
        call construct   (obs% l% y ,obs% oi ,'yi')
      endif
      obs% l% x = xi
      obs% l% y = 0._wp
      if (associated (obs% b% xb% s)) then
        if (.not.associated (obs% b% xl_xb% s))                &
          call construct    (obs% b% xl_xb, obs% ii ,'xl_xb')
        if (.not.associated (obs% b% ul_ub% s))                &
          call construct    (obs% b% ul_ub, obs% oi ,'ul_ub')
        obs% b% xl_xb = obs% l% x - obs% b% xb
      endif
    endif
    !-------------------------------------------------
    ! TSK_Y: default initialize dismissed observations
    !-------------------------------------------------
    if (iand (TSK_Y ,tsk) /= 0) then
      if (present (y)) y = invalid
    end if

    !-------------------------------------------
    ! TSK_INIT: we don't need the argument 'obs'
    !-------------------------------------------
    if (iand (TSK_INIT ,tsk) /= 0) then
      call process_obs (task)
    else
      !-------------------------------------------------------
      ! loop over observation 'boxes', call respective routine
      !-------------------------------------------------------
      nullify (ys)
      do ib=1,size(obs% o)
        if (llocal .and. obs% o(ib)% pe /= dace% pe) cycle
        if (present (y)) ys => y% s(ib)
        call obs_block (ob,obs,ib)
        if (present (xi)) then
          if (present (cols)) then
            call process_obs (task, ob, atm, cols(ib), xi% s(ib), ys, &
                              Jo, Jo_atm, state=state)
          else
            call process_obs (task, ob, atm, xi=xi%s(ib), y=ys, &
                              Jo=Jo, Jo_atm=Jo_atm, state=state)
          endif
        else
          if (present (cols)) then
            call process_obs (task, ob, atm, cols(ib), y=ys, &
                              Jo=Jo, Jo_atm=Jo_atm, state=state)
          else
            call process_obs (task, ob, atm, y=ys, &
                              Jo=Jo, Jo_atm=Jo_atm, state=state)
          endif
        endif
      end do
    endif

    !----------------------------------------------------
    ! TSK_Y, TSK_K : special handling for rttov_mult_prof
    !----------------------------------------------------
    if ((rttov_mult_prof .and. (IAND(TSK_Y+TSK_K, tsk) /= 0)) .or. &
        IAND(TSK_SETUP_OP, tsk) /= 0) then
      call process_tovs_mult (task, obs, atm=atm, xi=xi, y=y)
    end if

    !-------------------------
    ! TSK_SETUP_COL: broadcast
    !-------------------------
    if (iand(TSK_SETUP_COLS ,tsk) /= 0) then
      n = 0
      do ib=1,size(obs% o)
        ob% o => obs% o(ib)
        if (dace% pe == ob% o% pe) ob% o% n_spt = sum (ob% o% spot(:)% n_spt)
        if (ob% o% bcast) then
          call p_bcast (ob% o% spot(:)% n_spt, ob% o% pe)
          do io=1,ob% o% n_spot
            o => ob% o% spot (io)
            if (o% n_spt > 0 .and. o% int_type==ITY_MCOLS) then
              !-----------------------------------------------------------
              ! Change in # of model columns used by observation operator?
              !-----------------------------------------------------------
              if (associated (o% imcol)) then
                 if (size (o% imcol) /= o% n_spt) deallocate (o% imcol)
              end if

              if (.not.associated(o% imcol)) allocate(o% imcol(o%n_spt))
              if (n==0) n = size(transfer(o% imcol(1),(/' '/)))
              call p_bcast_derivedtype (o% imcol, n * o% n_spt, &
                                        ob% o% pe, (dace% comm) )
            endif
          end do
          ob% o% n_spt = sum (ob% o% spot(:)% n_spt)
        endif
        ob% o% mc% pe = ob% o% pe
      end do
      iatm_tovs_all = p_ior(iatm_tovs_all)
    end if

    !------------------------------------------
    ! TSK_SETUP_FULL: broadcast, set dimensions
    !------------------------------------------
    if (iand(TSK_SETUP_FULL ,tsk) /= 0) then
      do ib=1,size(obs% o)
        ob% o => obs% o(ib)
        call p_bcast     (ob% o% n_int,         ob% o% pe)
        if (.not.ob% o% bcast) cycle
        call p_bcast     (ob% o% spot(:)% i% i, ob% o% pe)
        call p_bcast     (ob% o% spot(:)% i% n, ob% o% pe)
        call p_bcast_ptr (ob% o% t_int,         ob% o% pe)
        call p_bcast_ptr (ob% o% lev,           ob% o% pe)
      end do
      call set_dec_info (obs)
      ! Must be called after TSK_SETUP_FULL...
      call print_tovs_levs(obs)
    end if

    !---------------------------
    ! TSK_K: set ul_ub = H xl_xb
    !---------------------------
    if (iand (TSK_K ,tsk) /= 0) then
      if (associated (obs% b% xl_xb% s)) obs%b% ul_ub = obs%l% H * obs%b% xl_xb
    endif

    !------------------------------------------------------------------------
    ! TSK_SHRINK: make dimension meta data of vectors and matrices consistent
    !------------------------------------------------------------------------
    if (iand (TSK_SHRINK ,tsk) /= 0) then
      call release_mem (obs% o)
      if (associated (obs% oi)) call set_dec_info (obs)
    endif

    call enter_function
    if (present(xi)) call delete_storage (xi)
    if (present(y))  call delete_storage (y)
    call leave_function

  end subroutine process_obs1
!------------------------------------------------------------------------------
  subroutine set_dec_info (obs)
  !------------------------------------------------
  ! set dimension meta data of vectors and matrices
  ! consistent with observations
  !------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs
    integer :: n_obs, n_int, n_spot, id, ib
    if (.not. associated (obs% si)) allocate (obs% si)
    if (.not. associated (obs% oi)) allocate (obs% oi)
    if (.not. associated (obs% ii)) allocate (obs% ii)
    if (.not. associated (obs% di)) allocate (obs% di)
    call destruct (obs% si)
    call destruct (obs% oi)
    call destruct (obs% ii)
    call destruct (obs% di)
    n_spot = sum (obs% o% n_spot)
    n_obs  = sum (obs% o% n_obs)
    n_int  = sum (obs% o% n_int)
    call construct (obs% si, n_spot, size(obs% o))
    call construct (obs% oi, n_obs,  size(obs% o))
    call construct (obs% ii, n_int,  size(obs% o))
    call construct (obs% di, 0,      size(obs% o))
    do ib=1,size(obs% o)
      call construct (obs% si% b(ib), obs% o(ib)% n_spot, obs% o(ib)% pe)
      call construct (obs% oi% b(ib), obs% o(ib)% n_obs,  obs% o(ib)% pe)
      call construct (obs% ii% b(ib), obs% o(ib)% n_int,  obs% o(ib)% pe)
      id = 0
      if (associated(obs% o(ib)% t_int))       &
        id = count(obs% o(ib)% t_int == OBS_DUM)
      call construct (obs% di% b(ib), id, obs% o(ib)% pe)
    end do
    obs% di% n = sum (obs% di% b(:)% n)
  end subroutine set_dec_info
!------------------------------------------------------------------------------
  subroutine process_obs0 (task, obs, atm, cols, xi, y, Jo, Jo_atm, state)
  integer            ,intent(in)              :: task   ! what to do
  type(t_obs_block)  ,intent(inout) ,optional :: obs    ! observations
  type(t_atm)        ,intent(in)    ,optional :: atm    ! atmospheric state
  type(t_cols)       ,intent(in)    ,optional :: cols   ! model columns
  type(t_vector_segm),intent(in)    ,optional :: xi     ! background(int.space)
  type(t_vector_segm),pointer       ,optional :: y      ! result (observ.space)
  real(wp)           ,intent(inout) ,optional :: Jo     ! obs. cost funct. Jo
  type(t_atm)        ,intent(inout) ,optional :: Jo_atm ! gradient:d Jo/d atm
  integer            ,intent(in)    ,optional :: state  ! status flag
    !----------------
    ! local variables
    !----------------
    integer                     :: i             ! loop indices
    integer                     :: n, no         ! array sizes for setup
    integer                     :: flag          ! for setup
    integer            ,pointer :: itmp (:)      ! temporary for setup
    real(wp)           ,pointer ::  tmp (:)      ! temporary for setup
    type(t_spot)       ,pointer :: s
    integer                     :: tsk           ! local copy of task
    type(t_vector_segm),pointer :: py            ! pointer to y
    type(t_obs)                 :: op1 (1)
    type(t_atm)        ,pointer :: patm
    type(t_atm)        ,target  :: dummy
    target                      :: atm

    !--------------------------------------------------------------------------
    ! TASKS:
    !
    ! TSK_INIT       M    initialize modules:
    !                     Set and allocate module variables.
    !                     Read files and namelists.
    !
    ! TSK_READ       M    read observations
    !
    ! TSK_SET_CHR    S    set characteristics of the observation
    !
    ! TSK_SETUP_COLS SC   setup columns:
    !                     Determine size of interpolation space
    !                     Determine type of interpolated quantities
    !                     Determine model columns required by the operator.
    !
    ! TSK_SETUP_FULL S B  setup description of PSAS-space:
    !                     Allocate interpolation space
    !                     Determine levels of interpolated quantities
    !
    ! TSK_R          ???  Set R (observational error covariance).
    !
    ! TSK_K          SB
    !                     Set H (tangent-linear observation operator).
    !
    ! TSK_XK              broadcast H
    !
    ! TSK_H          GB   run tangent linear operator
    !--------------------------------------------------------------------------
    nullify (py); if(present(y)) py => y

    nullify (dummy% grid)
    patm   => dummy;   if(present(atm))   patm     => atm

    tsk = task
    !=========
    ! printout
    !=========
    if (obs_debug > 0 .and. dace% lpio) then
      write(6,'(a)') repeat('=',79)
      write(6,'(a)')
      write(6,'(a)') 'process_obs : model fields present (atm):'
      if(associated (patm% t ))        write(6,'(a)')'  temperature'
      if(associated (patm% q ))        write(6,'(a)')'  specific humidity'
      if(associated (patm% pf))        write(6,'(a)')'  pressure (full levels)'
      if(associated (patm% ph))        write(6,'(a)')'  pressure (half levels)'
      if(associated (patm% ps))        write(6,'(a)')'  pressure (surface)'
      write(6,'(a)') 'process_obs : arguments present:'
      if(present(obs))                write(6,'(a)')'  obs'
      if(present(Jo))                 write(6,'(a)')'  Jo'
      if(present(Jo_atm))             write(6,'(a)')'  Jo_atm'
      if(present(cols))               write(6,'(a)')'  cols'
      write(6,'(a,i8)') 'process_obs : flags present:',tsk
      if(iand(tsk,TSK_INIT      )/=0) write(6,'(a)')'  TSK_INIT'
      if(iand(tsk,TSK_READ      )/=0) write(6,'(a)')'  TSK_READ'
      if(iand(tsk,TSK_SET_CHR   )/=0) write(6,'(a)')'  TSK_SET_CHR'
      if(iand(tsk,TSK_SETUP_COLS)/=0) write(6,'(a)')'  TSK_SETUP_COLS'
      if(iand(tsk,TSK_SETUP_FULL)/=0) write(6,'(a)')'  TSK_SETUP_FULL'
      if(iand(tsk,TSK_SETUP_FUL0)/=0) write(6,'(a)')'  TSK_SETUP_FUL0'
      if(iand(tsk,TSK_SETUP_OP  )/=0) write(6,'(a)')'  TSK_SETUP_OP'
      if(iand(tsk,TSK_SHRINK    )/=0) write(6,'(a)')'  TSK_SHRINK'
      if(iand(tsk,TSK_H         )/=0) write(6,'(a)')'  TSK_H'
      if(iand(tsk,TSK_K         )/=0) write(6,'(a)')'  TSK_K'
      if(iand(tsk,TSK_R         )/=0) write(6,'(a)')'  TSK_R'
    endif
    !===================
    ! Initialise modules
    !===================
    if (iand(tsk,TSK_INIT)/=0) then
      call read_obs_nml      ! namelist /OBSERVATIONS/
      call read_nml_report   ! namelist /REPORT/
      call read_nml_thin     ! namelist /THINNING/ (thinning parameters)
      flag = TSK_INIT
      call process_general
      call init_obs_err
      tsk = tsk - TSK_INIT
    endif
    !==================
    ! read observations
    !==================
    if (iand(tsk,TSK_READ)/=0) then
      if (dace% lpio) then
        write(6,*)
#if defined (_DACE_) && !defined (NOBUFR)
        write(6,*) 'reading BUFR on PE',p_readbufr
        write(6,*)
#endif
      endif
      !----------------------------------------
      ! artificial observations (from namelist)
      !----------------------------------------
      call read_nml_obs (obs% o)
      !-----------------------------------
      ! read observations from NetCDF file
      !-----------------------------------
      if (read_netcdf) call read_obs_netcdf (obs% o)
      !----------------------------
      ! read observations from BUFR
      !----------------------------
      if (read_bufr) call read_obs_bufr (obs% o)
      !-------------------------------------
      ! read observations from feedback file
      !-------------------------------------
      op1 = obs% o
      call read_feedback (pass=2, obs=op1, chk=.true., spec=.true.)
      call read_fdbk_rad('post',op1)
      obs% o = op1(1)
      !------------------------
      ! read other observations
      !------------------------
      flag = TSK_READ
      call process_general
      !-----------------------------
      ! release memory not used
      ! print memory usage
      !-----------------------------
      call release_mem   (obs% o)
!     call print_obs     (obs% o)
!     call print_summary (obs% o)
      tsk = tsk - TSK_READ
    endif
    if (tsk == 0) return
    !========================================
    ! allocate memory for feedback quantities
    ! set to zero
    !========================================
    if (iand (TSK_SET_CHR ,tsk) /= 0) then
      !========================================
      ! call specific setup routines
      ! to set observation characteristics
      !========================================
      flag = TSK_SET_CHR
      call setup_specific
    endif

    if (iand(TSK_SETUP_COLS ,tsk) /= 0) then
      !========================================
      ! call specific setup routines
      ! to determine required model columns
      ! broadcast spot% imcol
      !========================================
      if (dace% pe==obs% o%pe) then
        if (obs_debug>0) print *,'  task = TSK_SETUP_COLS'
        flag = TSK_SETUP_COLS
        allocate (obs% o% mc% idx (patm% grid% lbg(1):patm% grid% ubg(1), &
                                   patm% grid% lbg(2):patm% grid% ubg(2), &
                                   patm% grid% lbg(4):patm% grid% ubg(4), &
                                   obs% o% n_time                      ))

        allocate (obs% o% mc% c (2 * obs% o% n_spot))
        obs% o% mc% idx = 0
        obs% o% mc% n   = 0
        call setup_specific
        deallocate (obs% o% mc% idx)
      endif

    endif

    if (iand (TSK_SETUP_FULL ,tsk) /= 0) then

      if (obs_debug>0) print *,'  task = TSK_SETUP_FUL0'
      flag = TSK_SETUP_FUL0
      call setup_specific

      if (obs_debug>0) print *,'  task = TSK_SETUP_FULL'
      !==========================================
      ! allocate memory for interpolated space
      ! and call specific setup routines again
      ! broadcast interpolation space information
      !==========================================
      flag = TSK_SETUP_FULL
      n    = sum (obs% o% spot(1:obs% o% n_spot)% i% n)
      no = 0; if (associated(obs% o% t_int)) no = size (obs% o% t_int)
      itmp => obs% o% t_int
      tmp  => obs% o% lev
      allocate (obs% o% t_int(n))
      allocate (obs% o% lev  (n))
      obs% o% n_int = n
      n = 0
      do i=1,obs% o% n_spot
        if (no>0) &
         obs%o%t_int(                n+1:n                +obs%o%spot(i)%i% n)&
              = itmp(obs%o%spot(i)%i%i+1:obs%o%spot(i)%i%i+obs%o%spot(i)%i% n)
         obs% o%lev (                n+1:n                +obs%o%spot(i)%i% n)&
              =  tmp(obs%o%spot(i)%i%i+1:obs%o%spot(i)%i%i+obs%o%spot(i)%i% n)
        obs% o% spot(i)% i% i = n
        n = n + obs% o% spot(i)% i% n
      end do
      if (associated(itmp)) deallocate (itmp)
      if (associated( tmp)) deallocate ( tmp)

      !--------------------------------------------
      ! allocate space in R (obs. err. cov. matrix)
      !--------------------------------------------
      if (obs% o% pe == dace% pe) then
        if (associated(obs% R)) then
          if (.not.associated(obs% R% packed)) then
            obs% R% repr = CSC
            obs% R% n = 0
            obs% R% nonzero = 0
            obs% R% pe = dace% pe
            do i=1,obs% o% n_spot
              obs% R% n = obs% R% n + obs% o% spot(i)% o% n
              obs% R% nonzero = obs% R% nonzero + obs% o% spot(i)% nr
            end do
            obs% R% m = obs% R% n
            allocate (obs% R% packed (obs% R% nonzero))
            allocate (obs% R% ja  (obs% R% nonzero))
            allocate (obs% R% ia  (obs% R% n+1))
            obs% R% ia = 1
            obs% R% packed = -9._mp
            obs% R% ja     = -9
          endif
        endif
      endif
      call setup_specific
      if (obs% o% pe == dace% pe .and. associated(obs% R)) then
        obs% R% nonzero = obs% R% ia(obs% R% n+1) - 1
      endif
    endif

    !===============================
    ! fill in R (obs.err.cov.matrix)
    !===============================
    if (iand (TSK_R ,tsk) /= 0) then
      if (obs_debug>0) print *,'  task = TSK_R'
      flag = TSK_R
      if (obs% o% pe == dace% pe) then
        call setup_specific
        obs% R% nonzero = obs% R% ia(obs% R% n+1) - 1
        !----------------------------
        ! check for variables not set
        !----------------------------
!xxxxxxxxxxxx
!        if(any(obs% R% packed(1:obs% R% nonzero)<=0._wp)) &
!          call finish ('process_obs0,TSK_R','R <= 0.')
!xxxxxxxxxxxx
        if(any(obs% R% ja    (1:obs% R% nonzero)<=0    )) &
          call finish ('process_obs0,TSK_R','ja <= 0')
      endif
    endif

    !===========================
    ! remove unused observations
    !===========================
    if (iand(TSK_SHRINK ,tsk) /= 0) then
      flag = TSK_SHRINK
      if (.not.present(state)) &
        call finish("process_obs0 (TSK_SHRINK)","argument 'state' is required")
      call setup_specific
    endif

    !===================================
    ! initialise GPSRO/GB data structure
    !===================================
    if (iand (TSK_Y+TSK_K+TSK_YH ,tsk) /= 0) then
      if (obs% o% pe == dace% pe) then
        call set_occ (gf, patm% grid, cols, obs% o, xi)
        call set_std (patm% time)
      end if
    endif

    !===================================
    ! set up H (tangent linear operator)
    !===================================
    if (iand (TSK_K ,tsk) /= 0) then
      if (obs% o% pe == dace% pe) then
        flag = TSK_K
        !--------------------------------------------
        ! first time: allocate H with sufficient size
        !--------------------------------------------
        if (.not.associated(obs% H% packed)) then
          obs% H% repr = CSC
          obs% H% m = 0
          obs% H% n = 0
          obs% H% nonzero = 0
          obs% H% pe = dace% pe
          do i=1,obs% o% n_spot
            obs% H% m = obs% H% m + obs% o% spot(i)% o% n
            obs% H% n = obs% H% n + obs% o% spot(i)% i% n
            obs% H% nonzero = obs% H% nonzero + &
                              obs% o% spot(i)% o% n * obs% o% spot(i)% i% n
          end do
          call allocate_block (obs% H, repr=CSC, ns=obs% H% nonzero)
          obs% H% ia =  1
          obs% H% ja = -9
        endif
        !-----------------------------------------------------
        ! store non-zero coefficients in sparse representation
        !-----------------------------------------------------
        call setup_specific
        obs% H% nonzero = obs% H% ia(obs% H% n+1) - 1
        !-----------------------------------------
        ! resize storage for sparse representation
        !-----------------------------------------
        if (obs% H% nonzero < size (obs% H% packed)) &
          call reallocate (obs% H, ns=obs% H% nonzero)
      end if
    end if

    !=======================
    ! run nonlinear operator
    !=======================
    if (iand (TSK_Y ,tsk) /= 0) then
      if (associated (py)) then
        if (obs% o% pe == dace% pe) then
          do i=1,obs% o% n_spot
            if (obs% o% spot(i)% use% state <= STAT_DISMISS) cycle
            flag = TSK_Y
            call call_specific (i, y=y)
          end do
        end if
      end if
    endif

    !=================================
    ! run linear or nonlinear operator
    !=================================
    if (iand (TSK_YH+TSK_K ,tsk) /= 0) then
      if (associated (py)) then
        if (obs% o% pe == dace% pe) then
          !----------------
          ! linear operator
          !----------------
          y% x = obs% H * (xi% x - obs% xi% x) + obs% yi% x
          !------------------------------------
          ! if required run nonlinear operators
          !------------------------------------
          do i=1,obs% o% n_spot
            s => obs% o% spot(i)
            if (iand(s% char,CHR_NONL)/=0) then
              if (iand (TSK_YH ,tsk)/=0) then
                flag = TSK_Y
                call call_specific (i, y=y)
              endif
            endif
          end do
        end if
      end if
    endif

    !===================================
    ! deallocate GPSRO/GB data structure
    !===================================
    if (iand (TSK_Y+TSK_K+TSK_YH ,tsk) /= 0) then
      if (obs% o% pe == dace% pe) call destruct_occ
    endif

    !==================================================
    ! destruct TOVS obserror parameterization structure
    !==================================================
    if (iand (TSK_R ,tsk) /= 0) then
      if (obs% o% pe == dace% pe) call destruct_tovs_oe
    endif


!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------
    subroutine setup_specific
    !-----------------------------
    ! call specific setup routines
    !-----------------------------
      integer :: i
      if (associated (py)) then
        do i=1,obs% o% n_spot
          call call_specific (i, y=py)
        end do
      else
        do i=1,obs% o% n_spot
          call call_specific (i)
        end do
      end if
    end subroutine setup_specific
!------------------------------------------------------------------------------
    subroutine call_specific (i,y)
    integer             ,intent(in)              :: i
    type(t_vector_segm) ,intent(inout) ,optional :: y
    !-----------------------------
    ! call specific setup routines
    !-----------------------------
      select case (obs% o% spot(i)% hd% modtype)
      case (TEMP)
        call process_temp       (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (SYNOP)
        call process_synop      (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (SCATT)
        call process_scatt      (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (TOVS)
        call process_tovs       (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (SATEM)
        call process_satem      (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (GPSRO)
        call process_occ        (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (GPSGB)
        call process_std        (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (COSMO)
        call process_cosmo_conv (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (RADAR)
        call process_radar      (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (AMV)
        call process_amv        (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (AIREP)
        call process_airep      (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (WLIDAR)
        call process_wlidar     (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case (SOIL)
        call process_soil       (flag, obs%o%spot(i), obs, patm, cols, xi, y, &
                                 Jo, Jo_atm, state=state)
      case default
!       if(dace% lpio) write(0,*)&
!         'setup_obs: unknown module type',obs% o% spot(i)% hd% modtype
!       call finish ('call_specific','unknown module type')
      end select
    end subroutine call_specific
!------------------------------------------------------------------------------
    subroutine process_general
      if (any(rept_char% mod == TEMP  )) call process_temp       (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == SYNOP )) call process_synop      (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == SCATT )) call process_scatt      (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == TOVS  )) call process_tovs       (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == SATEM )) call process_satem      (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == GPSRO )) call process_occ        (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == GPSGB )) call process_std        (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == COSMO )) call process_cosmo_conv (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == RADAR )) call process_radar      (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == AMV   )) call process_amv        (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == AIREP )) call process_airep      (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == WLIDAR)) call process_wlidar     (flag ,obs=obs ,atm=patm)
      if (any(rept_char% mod == SOIL))   call process_soil       (flag ,obs=obs ,atm=patm)
    end subroutine process_general
!------------------------------------------------------------------------------
  end subroutine process_obs0
!==============================================================================
  subroutine nominal_height (obs, e_fi)
  type(t_obs_set) ,intent(inout)        :: obs  ! observation data
  type(t_vector)  ,intent(in) ,optional :: e_fi ! actual background error
  !-------------------------------------------------------------------
  ! Derive the nominal height of observations.
  ! In general (for in-situ observations) the vertical coordinate is
  ! defined in pressure coordinates and this height is returned.
  ! For other observations (currently radiances and radiooccultations)
  ! the vertical coordinate is derived from the maximum absolute value
  ! of H (the respective column of the Jacoby matrix)
  !-------------------------------------------------------------------
    integer                   :: ib         ! 'box'  index
    integer                   :: is         ! 'spot' index
    integer                   :: j          ! observ.(channel) index
    integer                   :: k          ! level index
    integer                   :: i          ! max.Jacobian level index
    integer                   :: i0, i1, i2 ! indices
    integer                   :: m          ! index
    integer                   :: method     ! copy of pl_method
    real(wp)                  :: e_bg_t     ! temperature weight
    real(wp)                  :: e_bg_rh    ! rel.hum.    weight
    logical                   :: p_log      ! copy of pl_log
    type(t_obs)  ,pointer     :: bi         ! 'box'  pointer
    type(t_spot) ,pointer     :: si         ! 'spot' pointer
    type(t_occ)               :: occ        ! GPSRO specific table
    type(t_ray)  ,pointer     :: ray(:)     ! GPSRO specific table
    real(wp)     ,allocatable :: H(:,:)     ! Jacoby matrix
    real(wp)     ,allocatable :: z(:)       ! height: p or ln(p)
    real(wp)                  :: w          ! weight for mean Jacobian
    real(wp)                  :: sw         ! sum of weights
    real(wp)                  :: wz         ! weight * Jacobian
    real(wp)                  :: swz        ! sum of weights*Jacobian
    real(wp)                  :: swz2       ! sum of w * Jacobian^2
    real(wp)                  :: wmax       ! max weighted Jacobian
    real(wp)                  :: wq         ! humidity weight

    !------------------------------
    ! loop over observation 'boxes'
    !------------------------------
    do ib = 1, size(obs% o)
      bi => obs% o(ib)
      if (bi% pe /= dace% pe) cycle
      !-----------------------------------------------------
      ! loop over single reports or satellite fields of view
      !-----------------------------------------------------
      do is = 1, bi% n_spot
        si => bi% spot(is)    ! 'spot' pointer
        i0 = si% o% i
        i1 = i0 + 1           ! indices to ..
        i2 = i0 + si% o% n    ! observation arrays
        select case (si% hd% obstype)
        !-------------
        ! in situ data
        !-------------
        case default
          !----------------------------------
          ! pass explicit pressure coordinate
          !----------------------------------
          where (bi% body(i1:i2)% lev_typ <  0   ) &
                 bi% body(i1:i2)% lev_typ =  VN_P
          where (bi% body(i1:i2)% lev_typ == VN_P) &
                 bi% body(i1:i2)% plev    =  bi% olev (i1:i2)
        !-------------------
        ! radio occultations
        !-------------------
        case (OT_GPSRO)
          !-------------------------------------------------
          ! from geometrical height of impact parameter (km)
          !-------------------------------------------------
          call load_occ (obs% o(ib), si, occ, ray)
          where (ray% geo% h > -1._wp)
            bi% body(i1:i2)% plev = p_h_usstd (ray% geo% h * 1000._wp)
          elsewhere
            bi% body(i1:i2)% plev = -1._wp
          endwhere
          deallocate (ray)
        !----------
        ! radiances
        !----------
        case (OT_RAD, OT_GPSGB)
          select case (si% hd% obstype)
          case (OT_RAD)
            method  = pl_method
            e_bg_t  = pl_e_bg_t
            e_bg_rh = pl_e_bg_rh
            p_log   = pl_log
          case (OT_GPSGB)
            if (std_method == 0) bi% body(i1:i2)% plev = si% ps_bg - 20000._wp
            if (std_method <= 0)                cycle
            if (std_method == 20)               cycle
!           if (si% use% state <= STAT_DISMISS) cycle
            method  = std_method
            e_bg_t  = 0._wp
            e_bg_rh = 1._wp
            p_log   = .false.
          end select
          if (method==4 .and. .not. present(e_fi))      &
            call finish('nominal_height','e_fi not present')
          !--------------------------------
          ! get Jacoby (sensitivity) matrix
          !--------------------------------
          if (.not. associated (obs% l% H% b)) &
            call finish('nominal_height','Jacobi matrix not set')
          allocate (H (si% o% n, si% i% n))
          allocate (z (si% i% n))
          H = sub_block (obs% l% H% b(ib,ib), si% o% i, si% o% n,     &
                                              si% i% i, si% i% n, FULL)
          if (p_log) then
            z =      bi% lev (si% i% i+1 : si% i% i + si% i% n)
          else
            z = exp (bi% lev (si% i% i+1 : si% i% i + si% i% n))
          endif
          !--------------------------
          ! loop over channel indices
          !--------------------------
          do j = 1, si% o% n
            !--------------------------------------------------
            ! different options for nominal height of radiances
            !--------------------------------------------------
            select case (method)
            !-------------------------------------------------------------
            ! pl_method == 1 : take maximum value of temperature Jacobian
            !-------------------------------------------------------------
            case (1)
              !-------------------------------------------------------------
              ! set sensitivities to zero if not with respect to temperature
              ! OBS_TV: temperature; OBS_RH: relative humidity
              !-------------------------------------------------------------
              where (bi% t_int (si% i% i + 1       &
                                :si% i% i + si% i% n) /= OBS_TV)  H (j,:) = 0._wp
              !----------------------------------------------------------------
              ! get index of maximum absolute value of the sensitivity function
              !----------------------------------------------------------------
              i = sum (maxloc (abs (H (j,:))))
              !---------------------------------------------------
              ! set return value (nominal pressure of observation)
              !---------------------------------------------------
              if (i==1) then
                 !---------------------------------------------------
                 ! for uppermost levels set very small pressure value
                 !---------------------------------------------------
                 bi% body(i0 + j)% plev = 0.000001
              else
                 !---------------------
                 ! derive p from log(p)
                 !---------------------
                 bi% body(i0 + j)% plev = exp (bi% lev (si% i% i + i))
              end if

            !-------------------------------------------------------------
            ! pl_method == 2 : take mean value of temp.+rel.hum. Jacobian
            !-------------------------------------------------------------
            case (2,3,4)
              !-----------------
              ! loop over levels
              !-----------------
              sw    = 0._wp        ! sum of weights
              swz   = 0._wp        ! sum of weights * Jacobian
              swz2  = 0._wp        ! sum of weights * Jacobian^2
              wmax  = 0._wp        ! max weighted Jacobian
              wq    = 0._wp        ! temperature weight
              i     = 0            ! maximum Jacobian level index
              m     = 1000         ! highest level with temperature data
              do k = 1, si% i% n
                select case (bi% t_int (si% i% i + k))
                case (OBS_TV)
                  select case (method)
                  case (2)
                    w   = abs (H (j,k)) * e_bg_t
                  case (3)
                    w   =     (H (j,k)  * e_bg_t) ** 2
                  case (4)
                    w   = abs (H (j,k)  * e_bg_t * e_fi% s(ib)% x(si% i% i + k))
                  end select
                  m   = min (m,k)
                  if (w > wmax) then
                    wmax = w
                    i    = k
                  endif
!print *,'### OBS_TV :',j,k,e_fi% s(ib)% x(si% i% i + k),H (j,k),w,exp (bi% lev (si% i% i + k)),sw,swz, bi% lev (si% i% i + k)
                case (OBS_RH)
                  select case (method)
                  case (2)
                    w   = abs (H (j,k)) * e_bg_rh
                  case (3)
                    w   =     (H (j,k)  * e_bg_rh) ** 2
                  case (4)
                    w   = abs (H (j,k)  * e_bg_rh * e_fi% s(ib)% x(si% i% i + k))
                  end select
                  wq  = wq + w
                case default
                  w = 0
                end select
                wz  = w * z (k)
                sw  = sw   + w
                swz = swz  + wz
                swz2 =swz2 + wz * z (k)
              end do
              if (si% hd% obstype == OT_RAD .and. i == m) then
                !---------------------------------------------------
                ! for uppermost levels set very small pressure value
                !---------------------------------------------------
                bi% body(i0 + j)% plev = 0.000001
                swz2 = 0._wp
              else
                !---------------------
                ! derive mean Jacobian
                !---------------------
                swz  =       swz  / sw
                swz2 = sqrt (swz2 / sw - swz * swz)
                if (p_log) then
                  bi% body(i0 + j)% plev = exp (swz)
                else
                  bi% body(i0 + j)% plev =      swz
                  swz2 = swz2 / swz
                endif
              end if
              !----------------------------------------------
              ! sensitivity with respect to humidity
              !----------------------------------------------
              wq = wq / sw
              bi% body(i0 + j)% plev_width = swz2
            case(5,6,7)
              if (si% hd% obstype == OT_RAD) then
                ! calculate height/width based on transmission in process_tovs_mult
              else
                call finish('nominal_height','invalid value for pl_method')
              end if
            case default
              call finish('nominal_height','invalid value for pl_method')
            end select
          end do
          deallocate (H)
          deallocate (z)
        end select
      end do
    end do
  end subroutine nominal_height
!==============================================================================

  subroutine set_time_slice (obs, last_time, tinc, n_slot, interpol)
  type (t_obs)  ,intent(inout) :: obs(:)    ! observation data
  type (t_time) ,intent(in)    :: last_time ! time of last slot
  type (t_time) ,intent(in)    :: tinc      ! time increment in between slots
  integer       ,intent(in)    :: n_slot    ! number of time slots
  integer       ,intent(in)    :: interpol  ! -1:interpolate 0:NN >0:no action

    integer  :: ib  ! observation box index
    integer  :: is  ! observation spot index
    real(wp) :: o_a ! time: observation - analysis    (seconds)
    real(wp) :: f_a ! time: first slot  - analysis    (seconds)
    real(wp) :: o_f ! time: observation - first slot  (seconds)
    real(wp) :: ti  ! time increment in between slots (seconds)
    real(wp) :: w
    integer  :: t

    if (interpol <= 0) then
      ti  = seconds (tinc)
      f_a = - ti * (n_slot - 1)
#ifdef __NEC__
    else                ! Arrghhh!
      ti  = -999._wp    ! Work around excessive optimization in nfort-3.0.5
#endif
    endif

    do ib = 1, size(obs)
      obs(ib)% n_time = n_slot
      if (obs(ib)% pe /= dace% pe) cycle
      do is = 1, obs(ib)% n_spot
        select case (interpol)
        case (1:)
          obs(ib)% spot(is)% i_time = 1       ! time slice
          obs(ib)% spot(is)% w_time = 0._wp   ! time slice fraction
        case default
          o_a = seconds (obs(ib)% spot(is)% actual_time - last_time)
          o_f = o_a - f_a
          w   = o_f / ti
          t   = floor (w)
          w   = w - t
          t   = t + 1
          if (t < 1) then
            t = 1
            w = 0._wp
          else if (t >= n_slot) then
            t = n_slot
            w = 0._wp
          endif
          if (interpol == 0) then
            if (w > 0.5_wp) t = t + 1
            w = 0._wp
!         else
!           if (t == n_slot) then
!             t = n_slot - 1
!             w = 1._wp
!           endif
          endif
          obs(ib)% spot(is)% i_time = t     ! time slice
          obs(ib)% spot(is)% w_time = w     ! time slice fraction
        end select
      end do
    end do
  end subroutine set_time_slice

!==============================================================================

  subroutine set_rules_cosmo
  !---------------------------------------------
  ! set namelist groups /rules/ default settings
  ! for KENDA (COSMO LETKF)
  !---------------------------------------------

    type (t_set) :: ass       ! derived type for observation selection

    ass% use     = STAT_ACTIVE_1
    ass% m_rej   = 3
    ass% sgm_oc  = (/ huge(1._sp),huge(1._sp), huge(1._sp),0.0_sp /)
    ass% sgm_fg  = (/      3._sp ,     0._sp , huge(1.0_sp) /)
    ass% bnd_fg  = (/-huge(1._sp),huge(1._sp)/)
    ass% bnd_obs =    huge(1._sp)
    ass% sgm_vq  = svqc
    ass% frm_vq  = vqc_form
    ass% prc     = PRC_NONE

    call new_rule ('allow everything',                       &
                   use=STAT_ACTIVE_1, verb=9,                &
                   t=ass, gp=ass, uv=ass, q=ass, p=ass, o=ass)

  end subroutine set_rules_cosmo

!------------------------------------------------------------------------------

  subroutine set_rules_icon
  !---------------------------------------------
  ! set namelist groups /rules/ default settings
  ! for the global data assimilation system
  !---------------------------------------------

    type (t_set) :: ass
    type (t_set) :: no
    type (t_set) :: gns

    ass% use     = STAT_ACTIVE_1
    ass% m_rej   = 3
    ass% sgm_oc  = (/ huge(1._sp),huge(1._sp), 5._sp,      0.0_sp /)
    ass% sgm_fg  = (/      5._sp ,     0._sp , huge(1.0_sp) /)
    ass% bnd_fg  = (/-huge(1._sp),huge(1._sp)/)
    ass% bnd_obs =    huge(1._sp)
    ass% sgm_vq  = svqc
    ass% frm_vq  = vqc_form
    ass% prc     = PRC_NONE

    no           = ass
    no%  use     = STAT_DISMISS

    gns          = ass
    gns% bnd_fg  = (/0._sp, 0.02_sp/)
!   gns% bnd_obs =          0.02_sp

    call new_rule ('forbid everything',                &
                   use=STAT_DISMISS, verb=9,           &
                   t=no, gp=no, uv=no, q=no, p=no, o=no)

!   in this case: TEMP and PILOT (WINDPROFILER later RADAR) are processed together
    call new_rule ('define "Kennzahlen" for TEMP',                  &
                    type=TEMP,  use=STAT_ACTIVE_1,                  &
                    db_kz= [ 508,509,510,511, 512,513,514,515,      &
                             516,517,518,519, 520,521,522,523,      &
                             524,525,526,527, 536,548,              &
                             764,765,766,767, 768,769,770,771,      &
                             776,777,778,779, 780,781,782,783, 792, &
                             10520,10521,10526,10527,10553,10574,   &
                             10600,10776,10777,10780,10782,10783,   &
                             10785,10516,10517,10570               ])

!   same for COSMO observation operators

    call new_rule ('allow all "Kennzahlen" for COSMO obs.oper.', &
                    type=COSMO,  use=STAT_ACTIVE_1               )

!  508 PILOT,A geopt.height
!  509 PILOT,B geopt.height
!  510 PILOT,C geopt.height
!  511 PILOT,D geopt.height

!  512 PILOT,A
!  513 PILOT,B
!  514 PILOT,C
!  515 PILOT,D
!  516 TEMP, A , mobile (Pegasus)
!  517 TEMP, B , mobile (Pegasus)
!  518 TEMP, C , mobile (Pegasus)
!  519 TEMP, D , mobile (Pegasus)

!  520 TEMP, A
!  521 TEMP, B
!  522 TEMP, C
!  523 TEMP, D

!10520 TEMP BUFR
!10521 TEMP BUFR (reduced)
!10526 TEMP BUFR high resolution
!10527 TEMP BUFR high resolution (reduced)
!10574 TEMP DESCENDING BUFR high resolution

!  524 TEMP, A , mobile
!  525 TEMP, B , mobile
!  526 TEMP, C , mobile
!  527 TEMP, D , mobile
!  536 TEMP, ABCD , merged <1983
!10516 TEMP mobile BUFR high resolution
!10517 TEMP mobile BUFR high resolution (reduced)
!10570 TEMP mobile DESCENDING BUFR high resolution

!  548 WINDPROF, USA
!  549 WINDPROF, USA
!  550 WINDPROF, USA
!  551 WINDPROF, USA
!  552 WINDPROF, USA
!  553 WINDPROF, u,v;    Lindenberg,Europe
!  554 WINDPROF,    ;t,w Lindenberg,Europe
!  555 WINDPROF, Japan
!  556 WINDPROF, RASS Deutschl

!10553 WINDPROFILER
!10600 VAD WIND PROFILE

!  764 PILOT SHIP, A geopt.height
!  765 PILOT SHIP, B geopt.height
!  766 PILOT SHIP, C geopt.height
!  767 PILOT SHIP, D geopt.height

!  768 PILOT SHIP, A
!  769 PILOT SHIP, B
!  770 PILOT SHIP, C
!  771 PILOT SHIP, D

!  776 TEMP SHIP, A
!  777 TEMP SHIP, B
!  778 TEMP SHIP, C
!  779 TEMP SHIP, D
!  792 TEMP SHIP, ABCD, merged <1983

!10776 TEMP SHIP BUFR
!10777 TEMP SHIP BUFR (reduced)
!10782 TEMP SHIP BUFR high resolution
!10783 TEMP SHIP BUFR high resolution (reduced)
!10785 TEMP SHIP DESCENDING BUFR high resolution

!  780 TEMP DROP, A
!  781 TEMP DROP, B
!  782 TEMP DROP, C
!  783 TEMP DROP, D
!10780 TEMP DROP BUFR

    call new_rule ('define "Kennzahlen" for SYNOP',          &
                    type=SYNOP, use=STAT_ACTIVE_1,           &
                    db_kz= [ 0, 1, 5, 128, 131, 156, 166,    &
                             170, 182, 256, 384, 385, 386,   &
                             10000,10005,10015,10128,10143,  &
                             10150,10158,10170,              &
                             10256,10384,10385              ])

!    0 Bodenmeldungen SYNOP sect.1-4
!    1 Bodenmeldungen METAR
!    5 Bodenmeldungen SYNOP sect.5 (German national SYNOP BUFR)
!    9 Bodenmeldungen Pseudo-SYNOP (aus TEMP)
!  128 Bodenmeldungen SYNOP sect.1-3 automatisch
!  131 Bodenmeldungen SYNOP, recoded from METAR, USA
!  170 Bodenmeldungen SWIS,  road weather stations
!  182 Bodenmeldungen SNOW
!  256 Bodenmeldungen SHIP, manuell
!  265 Bodenmeldungen Pseudo-SYNOP (aus TEMP)
!  384 Bodenmeldungen SHIP, automatisch
!  385 Bodenmeldungen DRIBU, DRIFTER, BUOY
!  386                PAOB
!10000 Bodenmeldungen SYNOP sect.1-4 manuell     (new BUFR format)
!10005 Bodenmeldungen SYNOP sect.5 (German national SYNOP new BUFR)
!10015 Bodenmeldungen SYNOP          manuell     (BUFR, WIGOS ID)
!10128 Bodenmeldungen SYNOP sect.1-3 automatisch (new BUFR format)
!10143 Bodenmeldungen SYNOP          automatisch (BUFR, WIGOS ID)
!10150 Bodenmeldungen SYNOP          automatisch (BUFR, SEEMHEM )
!10158 Bodenmeldungen CMAN  coastal stations (USA;new BUFR format)
!10256 Bodenmeldungen SHIP           Manuell     (new BUFR format)
!10384 Bodenmeldungen SHIP           automatisch (new BUFR format)
!10385 Bodenmeldungen DRIBU, DRIFTER, BUOY       (new BUFR format)

!   call new_rule ('define "Kennzahlen" for SCATT',        &
!                   type=SCATT, use=STAT_ACTIVE_1,         &
!                   db_kz= [ 1697, 1698, 1699, 1700, 1701, 1702 ])

! 1697 Scatterometer  QUICKSCAT as DRIBU
! 1697 Scatterometer  QUICKSCAT from SKY in new format(needs some adaption)
! 1698 Scatterometer  ASCAT     from SKY in new format(from Eumetsat)
! 1699 Scatterometer  ASCAT     as DRIBU
! 1699 Scatterometer  ASCAT     from SKY in new format(needs some adaption)
! 1700 Scatterometer  OSCAT/RSCAT from SKY in new format(needs some adaption)
! 1701 Scatterometer  JASON-2/3
! 1702 Scatterometer  SARAL

    call new_rule ('define "Kennzahlen" for AMV',&
                     type=AMV, use=STAT_ACTIVE_1,&
                     db_kz= (/DK_SATOB_2,        &! SATOB section 2
                              DK_SATOB_3,        &! SATOB section 3
                              1674,1675,1677,    &! SATOB section 4,5,6 &
                              DK_AMV_EUMETSAT,   &! AMV, EUMETSAT
                              DK_AMV_GOES,       &! AMV, GOES
                              DK_AMV_MODIS,      &! AMV, MODIS
                              DK_AMV_MODIS_OLD,  &! AMV, MODIS
                              DK_AMV_MODIS_A,    &! AMV, MODIS
                              DK_AMV_MTV,        &! AMV, JAPAN
                              DK_AMV_NOAA,       &! AMV, NOAA polar orb.
                              DK_AMV_SENTINEL    &! AMV, SENTINEL (EUMETSAT)
                              /))
! 20080826: disable for routine equivalence:
!                             DK_AMV_FY_X,       &! AMV, CHINA
!                             DK_AMV_NOAA/))      ! AMV, NOAA polar orb.
! 20080515: add DK_AMV_MTV
! 20080515: for testing purposes: activate  DK_AMV_FY_X DK_AMV_NOAA
! MODIS : archive <200705: DK_AMV_MODIS_A   kz=1706 (Ascii format
! MODIS : globus  >200705: DK_AMV_MODIS_OLD kz=1705 (Bufr  format; identical with DK_AMV_GOES=1705)
! MODIS : SKY     >200804: DK_AMV_MODIS     kz=1710 (Bufr  format)
! don't use
!                             DK_AMV_FY_X/))      ! chinese satellite

! 1672 AMV Section 2
! 1673 AMV Section 3
! 1674 AMV Section 4
! 1675 AMV Section 5
! 1677 AMV Section 7
! 1704 AMV EUMETSAT geostationary sat.                  DK_AMV_EUMETSAT
! 1705 AMV NOAA/NESIDS              geostationary sat.  DK_AMV_GOES
! 1706 AMV NOAA/NESIDS MODIS(ASCII) polar sat.
! 1707 AMV CHINA                    geostationary sat.  DK_AMV_FY_X
! 1708 AMV JAPAN                    geostationary sat.  DK_AMV_MTV
! 1709 AMV NOAA/NESIDS NOAAxx       polar sat.          DK_AMV_NOAA
! 1710 AMV NOAA/NESIDS MODIS(BUFR)  polar sat.          DK_AMV_MODIS
! 1711,1712,1713,1714,              cloud,temp.,upper air humidity,ozone

!   routine (pegasus)
!   call new_rule ('define "Kennzahlen" for AIREP',&
!                    type=AIREP, use=STAT_ACTIVE_1,&
!                    db_kz= (/DK_AMDAR,          &! AMDAR
!                             DK_AIREP,          &! AIREP
!                             DK_ACARS_USA,      &! ACARS-Daten USA
!                             DK_ACARS_EU/))      ! ACARS-Daten EUROPA

!   routine (SKY    )
    call new_rule ('define "Kennzahlen" for AIREP',&
                     type=AIREP, use=STAT_ACTIVE_1,&
                     db_kz= (/DK_AMDAR,          &! AMDAR
                              DK_AIREP,          &! AIREP
                              DK_ACARS_USA,      &! ACARS-Daten USA
                              DK_ACARS_EU,       &! ACARS-Daten EUROPA
                              DK_ACARS_CH,       &! ACARS-Daten China (SKY)
                              DK_ACARS,          &! ACARS-Daten sonstige(SKY)
                              DK_ACARS_SINGLE,   &! ACARS-Daten single-level
                              DK_MODES_OLD,      &! MODES-Daten EUROPA
                              DK_MODES         /))! MODES-Daten EUROPA

! don't use
! redundant with DK_ACARS_EU: DK_ACARS_LH,       &! ACARS-Daten Lufthansa
! DK_CODAR          528  In-Situ Beob.,  CODAR
! DK_AMDAR          529  In-Situ Beob.,  AMDAR
! DK_AIREP          530  In-Situ Beob.,  AIREP
! DK_ACARS_LH       532  In-Situ Beob.,  ACARS-Daten Lufthansa
! DK_ACARS_USA      533  In-Situ Beob.,  ACARS-Daten USA
! DK_ACARS_EU       534  In-Situ Beob.,  ACARS-Daten EUROPA (Bracknell)
! DK_ACARS_CH       535, In-Situ Beob.,  ACARS-Daten China (SKY)
! DK_ACARS          535  In-Situ Beob.,  ACARS-Daten sonstige(PEGASUS)
! DK_ACARS          538  In-Situ Beob.,  ACARS-Daten sonstige(SKY)
! DK_ACARS_SINGLE 10532  In-Situ Beob.,  ACARS single-level (unified format)
! DK_MODES        10534  In-Situ Beob.,  MODES-Daten

    call new_rule ('define "Kennzahlen" for TOVS',&
                     type=TOVS, use=STAT_ACTIVE_1,&
                     db_kz= (/-1/),               &! undefined
!                    c=(/(ass,i=1,nc)/))
                     o=ass)

!    call new_rule ('define "Kennzahlen" for GPSGB',&
!                     type=GPSGB, use=STAT_ACTIVE_1,&
!                     db_kz= (/94,-1/))             ! undefined

    call new_rule ('define "Kennzahlen" for GPSGB',&
                     type=GPSGB, use=STAT_ACTIVE_1,&
                     db_kz= (/94,95,-1/))           ! ZTD + STD BUFR

    call new_rule ('define "Kennzahlen" for GPSRO',&
                     type=GPSRO, use=STAT_ACTIVE_1,&
!                    db_kz= (/-1/))                   ! undefined
                     db_kz= (/1694, 1695/))           ! DWD BUFR

    call new_rule ('define "Kennzahlen" for SATEM',&
                     type=SATEM, use=STAT_ACTIVE_1,&
!                    db_kz= (/-1/))                   ! undefined
                     db_kz= (/1794/))                 ! DWD BUFR

    call new_rule ('specific parameters for TEMP gp,rh',&
                    type=TEMP,                          &
                    bsurf=10000._wp,                    &
                    verb=0, t=ass, gp=ass, q=ass, p=no)

    call new_rule ('specific parameters for TEMP wind',&
                    type=TEMP,                         &
                    bsurf=2500._wp,                    &
                    verb=0, uv=ass, p=no)

    call new_rule ('specific parameters for SATEM',&
         type=SATEM,                               &
         verb=0, t=ass, gp=ass, q=ass, p=no)

!    call new_rule ('specific parameters for SYNOP gp,rh',&
!                   type=SYNOP,                           &
!                   bsurf=10000._wp,                       &
!                   asurf=10000._wp,                       &
!                   verb=0, msl=1, t=no, gp=ass, q=ass, p=no)
!
!    call new_rule ('specific parameters for SYNOP wind',&
!                   type=SYNOP,                          &
!                   bsurf=10000._wp,                      &
!                   asurf=10000._wp,                       &
!                   verb=0, msl=1, t=no, uv=ass, p=no)
!
!    call new_rule ('specific parameters for SYNOP-Land',&
!                    type=SYNOP, bf_type=WMOA_SURF_LAND, &
!                    bsurf=10000._wp,                     &
!                    asurf=10000._wp,                       &
!                    verb=0, msl=1, t=no, gp=ass, uv=no, q=ass, p=no)
!
!    call new_rule ('specific parameters for SYNOP',                 &
!                    type=SYNOP,                                     &
!                    bsurf=10000._wp,                                &
!                    asurf=10000._wp,                                &
!                    verb=0, msl=1, t=no, gp=ass, uv=ass, q=ass, p=no)

    call new_rule ('specific parameters for AMV', &
                    type=AMV,                     &
                    verb=0, uv=ass)

    call new_rule ('specific parameters for AIREP', &
                    type=AIREP,                     &
                    verb=0, uv=ass, t=ass, q=ass)

    call new_rule ('specific parameters for GPSRO', &
                    type=GPSRO,                     &
                    verb=0, o=gns)

    call new_rule ('specific parameters for GPSGB', &
                    type=GPSGB,                     &
                    verb=0, o=ass)

  end subroutine set_rules_icon

!==============================================================================
end module mo_obs
