!
!+ Routines specific to scatterometer and altimeter observations (SCATT)
!
MODULE mo_scatt
!
! Description:
!   Routines specific to scatterometer and altimeter observations (SCATT).
!   Split off from previously common module mo_synop.
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
! V1_4         2009/03/26 Alexander Cress
!  Read scatterometer data from NetCDF
! V1_5         2009/05/25 Alexander Cress
!  fix scatterometer bias correction
! V1_6         2009/06/10 Alexander Cress
!  subroutine mleqc_qscat: check for invalid index to proposed wind
! V1_7         2009/08/24 Andreas Rhodin
!  include dbkz=10385 in array 'kzahl' (BUOY, new BUFR format)
! V1_8         2009/12/09 Andreas Rhodin
!  optimizations, cleanup, workaround bug in xlf V12.1
! V1_9         2010/04/20 Harald Anlauf
!  qscat_bcor: improve numerical consistency between platforms/compilers
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  read_synop_netcdf: reduce diagnostic output for netcdf_verb==0
! V1_13        2011/11/01 Harald Anlauf
!  changes for BUOYs, METARs, OCEANSAT-2, TEMP merge; optimisations for NEC SX
! V1_14        2011/11/08 Harald Anlauf
!  Cleanup quality control for OSCAT
! V1_15        2011/12/06 Andreas Rhodin
!  remove unused component t_datum% u
! V1_16        2011/12/09 Alexander Cress
!  SCATTEROMETER: store fov in spti% phase
! V1_20        2012-06-18 Harald Anlauf
!  read_synop_netcdf: if pmsl (mpppp) is reported as 0, set to invalid value
!  Oscat: reject data if windspeed is <3m/s or >30m/s
! V1_22        2013-02-13 Harald Anlauf
!  Disambiguation of synops and ships via bufr_type
!  Enable ships reporting as "SHIP"
!  Bugfix for "new buoys" (KZ=10385) in the Antarctic
!  Separate treatment of ship and buoys
!  Additional diagnostics for double entries
! V1_24        2013/04/26 Andreas Rhodin
!  SCATTerometer: store satid in spot%hd%satid (required for correct thinning)
! V1_27        2013-11-08 Harald Anlauf
!  Dismiss SYNOPs with lat=0,lon=0
! V1_28        2014/02/26 Andreas Rhodin
!  new interface to new_int
! V1_29        2014/04/02 Harald Anlauf
!  Reduce default verbosity for scatterometer data
! V1_31        2014-08-21 Andreas Rhodin
!  changes for ECMWF SYNOP (BUFR2)NetCDF input
! V1_37        2014-12-23 Alexander Cress
!  preference BUFR SYNOP reports
! V1_40        2015-02-27 Harald Anlauf
!  changes for new SHIP BUFR reports (A.Cress)
!  bugfix for reading scatterometer reports in NetCDF (R.Faulwetter)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_45        2015-12-15 Harald Anlauf
!  Implement windspeed observations for Jason-2; Preparations for SARAL/Altika
! V1_46        2016-02-05 Harald Anlauf
!  Improve handling of SYNOP wind obs; enhance check of multiple reports;
!  check against proper missing value for station height;
!  Fix DRIBU station id; extend SYNOP/DRIBU code for monitoring of T2m
! V1_47        2016-06-06 Harald Anlauf
!  scatterometer bias correction;
!  extend altimetry data QC (Jason-2) to ECMWF recommendations;
!  handle Jason-3, Sentinel 3A
! V1_48        2016-10-06 Andreas Rhodin
!  passive monitoring for RR,GUST,T2M,FF,DD
!  RH2M assimilation: option to use model first guess
! V1_49        2016-10-25 Andreas Rhodin
!  extend sanity check for reported SYNOP pressure; bug fix for MEC
! V1_50        2017-01-09 Harald Anlauf
!  add checks for altimeter; changes to run 3dvar with COSMO/COMET data.
! V1_51        2017-02-24 Andreas Rhodin
!  changes to pass observations from COSMO fof-files through the 3dvar
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD 2000-2008  original source
! Oliver Schmid   DWD       2005       new obs data type
! Gerhard Paul    DWD       2008       input in NetCDF format
! Harald Anlauf   DWD       2008       modifications for SX8
!==============================================================================

!-------------------------------------------------
! uncomment to exit in case of unknown BUFR codes:
!
!#define CHECKCODES
!-------------------------------------------------

  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp, sp, i2, i8  ! real, integer kind parameters
  use mo_namelist,   only: position_nml,  &! position namelist
                           nnml,          &! namelist Fortran unit number
                           POSITIONED      ! ok    code from position_nml
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_bcast         ! broadcast routine
  use mo_dace_string,only: char3           ! conversion: int -> char(len=3)
  use mo_obs_set,    only: t_obs_block     ! obs. data type
  use mo_exception,  only: finish          ! abort on error condition
  use mo_atm_state,  only: t_atm           ! atmospheric state data type
  use mo_t_datum,    only: t_datum,       &! data type for one observed datum
                           inv_datum,     &! constant for invalid datum
                           set_datum,     &! set t_datum% o  (observed value)
                           set_qbits,     &! set t_datum% qc (quality bits)
                           print,         &! generic print routine
                           rvind,         &! invalid value
                           SRC_DER,       &! derived  value flag
                           QC_OK           ! quality flag: OK
!                          QC_MISS,       &!               missing
!                          QC_NOUSE,      &!               do not use
!                          QC_CLIM         !               climatological range
  use mo_t_use,      only: t_use,         &! status variable data type
                           use_0,         &! default values of type use
                           decr_use,      &! decrease the state of a datum
!                          STAT_ACTIVE,   &!
                           STAT_DISMISS,  &!
                           STAT_PASSIVE,  &!
!                          STAT_OBS_ONLY, &!
                           CHK_NOTUSED,   &!
                           CHK_INSDAT,    &!
                           CHK_DOMAIN,    &!
                           CHK_DATASET,   &!
!                          CHK_QI,        &! quality index
                           CHK_NOIMPL      !
  use mo_wmo_tables, only: WMO0_ECMWF      ! generating center
  use mo_obs_tables, only: decr_rpt_use,  &!
                           idb_dbk,       &! index in table rept_stat
!                          rept_char,     &! observation type characteristics
                           check_report_0,&! init. flags, standard checks
                           check_report_1  !
  use mo_t_obs,      only: t_obs,         &! observation data type
                           t_spot,        &! spot data type
                           t_head,        &! observation data type
                           derive_dbkz,   &! derive DBKZ if not present
                           new_par,       &! request space in parameter array
                           new_spot,      &! request space in metadata  array
                           new_obs,       &! request space in observation array
                           new_int,       &! request space in int.obs.    array
                           set_xuv,       &! set unit vectors, solar zenith
                           invalid,       &! invalid observation value (real)
!                          empty_spot,    &! empty spot
                           set_vqc_insitu,&! subroutine to set VQC bounds
                           shrink_report, &! remove passive observations
                           monitor_ff,    &! flag to monitor wind speed
                           monitor_dd,    &! flag to monitor wind direction
                           source,        &! list   of Report source files
                           SCATT,         &! module flag value
!                          CHR_ID,        &! H is the identity operator
                           CHR_NONL,      &! H is a non-linear operator
                           TSK_INIT,      &! task flag: initialize modules
                           TSK_READ,      &!   read observations
                           TSK_SET_CHR,   &!   set observation characteristics
                           TSK_SETUP_COLS,&!   determine model columns required
                           TSK_SETUP_FULL,&!   setup descriptionof PSAS-space
                           TSK_SETUP_FUL0,&!   not used for SCATT
                           TSK_SHRINK,    &!  release unused obs. in report
                           TSK_R,         &!   setup observational errors
                           TSK_K,         &!   evaluate linear operator
                           TSK_Y,         &!   evaluate nonlinear operator
                           ITY_ICOL,      &! interpolation type: column
                           netcdf_verb,   &! verbosity of NetCDF decoding
                           OBS_U,OBS_V,OBS_DUM
  use mo_scatt_bc   ,only: read_scatt_bc_nml! read scatterometer bc namelist
  use mo_fdbk_tables,only: VN_U, VN_U10M, &!          wind component code
                           VN_V, VN_V10M, &!                         code
                           VN_FF,         &!          wind speed     code
                           VN_DD,         &!          wind direction code
                           OT_DRIBU,      &! data type specification
                           OT_SCATT        ! data type specification
  use mo_t_col,      only: t_cols,        &! model columns data type
                           COL_UV          !
  use mo_time,       only: init_time,     &! initialise time data type
                           operator(<),   &! compare  times
                           operator(>=),  &! compare  times
                           operator(==),  &! compare  times
                           operator (-),  &! subtract times
!                          hours,         &! derive hours   from time variable
!                          minutes,       &! derive minutes from time variable
                           chhmm,         &!
                           print,         &! generic print routine
                           cyyyymmddhhmmss ! derive string from time
  use mo_dec_matrix, only: t_vector_segm, &! vector segment data type
                           mp              ! kind parameter for matrix elements
  use mo_physics,    only: fd_uv,         &! derive ff, dd from u,v
                           pi,            &! die Zahl pi
                           d2r             ! pi/180, factor degree  -> radians
  use mo_grid_intpol,only: idx_init
  use mo_bufr_dwd,   only: t_bufr,        &! BUFR record data type
                           inv_bufr,      &! indicator for invalid value
!                          nvind,         &! integer missing value indicator
#ifdef  CHECKCODES
                           bufr_get_entry_texts, &!
                           bufr_get_entry_units, &!
#endif
                           bufr_print_sections,  &! print sections 0-4
                           bufr_print_subset,    &! print body of BUFR record
                           bufr_get_character     ! decode character entry
  use mo_obstypes,only:    t_obsid,       &! observation id table entry
                           obstype_dbkz,  &! derive obsids from dbkz
                           obstype_bufr    ! derive obsids from bufr type
  use mo_obs_err,    only: scatt_obs_err   ! set SCATT obs.err.
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only: ncid,          &! NetCDF file id
                           dimids_max,    &! max number of NetCDF dimension ids
                           imissing,      &! NetCDF _FillValue for integer
                           rmissing,      &! NetCDF _FillValue for reals
                           s2ikz,         &! DWD-internal classifier
                           s1cat,         &! data category
                           s1catls,       &! local data sub category
                           s1cent,        &! data centre
                           stime,         &! header observation time (section1)
                           db_time,       &! data bank time
                           s1cents,       &! data sub centre
                           s1updat,       &! update sequence no.
                           mlah,          &! latitude
                           mloh,          &! longitude
                           obs_time,      &! body observation time
                           ystidn,        &! any type of station identifier as variable
                           istidn          ! WMO numeric station number combined
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,        only: nf90_Inquire,          &!
                           nf90_Inquire_Dimension,&!
                           nf90_Inquire_Variable, &!
                           nf90_inq_varid,        &!
                           nf90_get_var,          &!
                           NF90_MAX_NAME,         &!
                           NF90_NOERR,            &!
                           NF90_INT,              &!
                           NF90_BYTE,             &!
                           NF90_SHORT,            &!
                           NF90_FLOAT,            &!
                           NF90_DOUBLE             !
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_scatt
  public :: read_scatt_bufr
  public :: read_scatt_netcdf ! read observations from netCDF file
  public :: t_scatt           ! SCATT observation          data type
  public :: inva_sca          ! invalid SCATT obs          t_scatt constant
  public :: check_store_scatt ! subroutine passed to BUFR read routine
  public :: use_ps_model      ! invalid value: use model surface pressure
  public :: read_scatt_nml    ! read namelist /SCATT_OBS/
!------------------------------------------------------------------------------
  !======================
  ! data type definitions
  !======================
  type t_scatt
    !-------------------------------------------------------------
    type (t_datum) :: total  ! quality control flag
    !-------------------------------------------------------------
    type (t_datum) :: z      ! height of station used     [m]
    type (t_datum) :: z_msl  ! height mean sea level      [m]       ! always 0
    type (t_datum) :: p_ref  ! reference pressure         [Pa]
    !-------------------------------------------------------------
    type (t_datum) :: ff     ! wind speed                 [m/s]
    type (t_datum) :: dd     ! wind direction             [deg]
    type (t_datum) :: uu     ! wind component             [m/s]
    type (t_datum) :: vv     ! wind component             [m/s]
    !-------------------------------------------------------------
!   type (t_datum) :: ff2    ! wind speed      alternative[m/s]
!   type (t_datum) :: dd2    ! wind direction  choice     [deg]
!   type (t_datum) :: uu2    ! wind component  for        [m/s]
!   type (t_datum) :: vv2    ! wind component  scatt.     [m/s]
    !-------------------------------------------------------------
!   type (t_datum) :: ssm    ! surface soil moisture      [%]
!   type (t_datum) :: essm   ! error estimate of ssm      [%]
    !-------------------------------------------------------------
  end type t_scatt

  type (t_scatt) ,parameter :: inva_sca =                 &
        t_scatt (inv_datum,inv_datum,inv_datum,inv_datum, &
                 inv_datum,inv_datum,inv_datum,inv_datum  )! ff,dd,uu,vv
!                inv_datum,inv_datum,inv_datum,inv_datum,  ! ff2,..,vv2

! integer, parameter :: scatt_int_size = size (transfer (inva_sca,(/0/)))
  integer            :: scatt_int_size = 0
  real(sp),parameter :: use_ps_model = 989898._sp ! use model surface pressure
!------------------------------------------------------------------------------
  !===============================
  ! List of "Datenbank-Kennzahlen"
  !===============================
  integer, parameter :: kz_scat (8) = [ 1697, 1698, 1699, 1700, 1701, 1702, 1770, 1780 ]
!------------------------------------------------------------------------------
  !-------------------
  ! Namelist SCATT_OBS
  !-------------------
! integer  :: zls_to_z   = 3         ! set z from land surface   height (see
! integer  :: zmsl_to_z  = 4         ! set z from mean sea level height  below)
  integer  :: scatt_proc = 0         ! processing of scatterometer observations
                                     ! 0: use wind proposed by producer of data
                                     ! 1: use wind consistent with first guess
                                     ! 2: chose wind in analysis step
  logical  :: prt_data   = .false.   ! print data
  integer  :: verbose    = 0         ! Verbosity level of consistency checks
  !
  ! values for zls_to_z or zmsl_to_z :
  !   0: use station height
  !   1: use surface height if station height is not present
  !   2: use surface height if surface height < station height
  !   3: use surface height if surface height > station height
  !   4: use surface height if surface height is present
  !   5: use surface height
  !
  namelist  /SCATT_OBS/ scatt_proc, prt_data, verbose
  !-------------------------------------------
  ! Module parameters used in NetCDF interface
  !-------------------------------------------
  integer, parameter :: TY_I  = 1       ! type of field: integer
  integer, parameter :: TY_F  = 2       ! type of field: real (float)
  integer, parameter :: NDIM1 = 1       ! dimension of field: 1
  integer, parameter :: NDIM2 = 2       ! dimension of field: 2
!==============================================================================
contains ! subroutines
!=====================
  subroutine process_scatt (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
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

    !----------------
    ! local variables
    !----------------
    integer                :: ni         ! actual number of coef. for interpol.
    integer                :: i, k       ! loop index
    integer                :: tsk        ! task
    integer(i8)            :: iatm       ! model column flag
    integer       ,pointer :: ty (:)     ! unpacked types
    logical                :: change     ! argument from shrink_report
    integer                :: iii, ioi   !+++ work around NEC SX bug
    real(wp)               :: ff, dd     ! wind speed, direction
    real(wp)               :: ff_x, ff_y ! derivatives
    real(wp)               :: dd_x, dd_y ! derivatives
    integer                :: i_u,i_v,i_d! interpolation space indices
    integer                :: o_u,o_v    ! observation   space indices
    integer                :: o_d,o_f    ! observation   space indices
    real(wp)               :: sf         ! scaling factor, wind speed profile
    !---------------------------
    ! process optional arguments
    !---------------------------
    !  TSK_INIT       ! initialize modules
    !  TSK_READ       ! read observations
    !  TSK_SET_CHR    ! set observation characteristics
    !  TSK_SETUP_COLS ! setup columns
    !  TSK_SETUP_FUL0 ! setup description of PSAS-space
    !  TSK_SETUP_FULL ! setup description of PSAS-space
    !  TSK_Y          ! evaluate nonlinear operator
    !  TSK_K          ! evaluate linear operator
    !==================================================
    ! skip tasks not required for this observation type
    !==================================================
    tsk = task
    tsk = iand (task, not (&
          TSK_READ))        ! Input is read by separate BUFR reading routine
    if (tsk == 0) return
    !===============
    ! initialisation
    !===============
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_scatt_nml
      call read_scatt_bc_nml
      tsk = tsk - TSK_INIT
      if (tsk == 0) return
    endif

    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      tsk = tsk - TSK_SET_CHR
      spot% int_type  = ITY_ICOL
      spot% cost      = 1._wp
      spot% nr        = spot% o%n
      spot% char      = CHR_NONL
!     if (spot% hd% obstype == OT_SCATT .and. spot% nr == 1) then
!        spot% char   = CHR_NONL        ! Altimeter data: wind speed
!     else
!        spot% char   = CHR_ID
!     end if
      if (tsk == 0) return
    endif

    !==========================================
    ! tsk == TSK_SHRINK:
    ! release unused observations in the report
    !==========================================
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change)
      if (.not. change .or. spot% i%n==0) then
        tsk = tsk - TSK_SHRINK
        if (tsk == 0) return
      endif
    endif

    !===========
    ! PSAS-space
    !===========
    if (iand (TSK_SETUP_FUL0+TSK_SHRINK,tsk) /= 0) then
!     if (dace% pe==obs% o% pe) then
        ioi = spot% o% i
        i_u = 0; i_v = 0; i_d = 0
        do i = 1, spot% o% n
          select case (obs% o% varno (ioi+i))
          case (VN_U, VN_U10M, VN_V, VN_V10M, VN_FF, VN_DD)
            if (i_u == 0) then
              i_u = 1
            endif
            if (i_v == 0) then
              i_v = 1
            endif
            o_v   = i
          end select
        end do

        ni = 0
        if (i_u > 0) then
          ni = ni + 1; i_u = ni
        endif
        if (i_v > 0) then
          ni = ni + 1; i_v = ni
        endif
        if (i_d > 0) then
          ni = ni + 1; i_d = ni
        endif

!       if (ni == 0) then
!         ni  = ni + 1
!         i_d = ni
!       endif

        if (iand (TSK_SHRINK,tsk) /= 0) then
           k = 0
           if (ni == 0) then
             spot% i% n = ni
           else
             do i = 1, spot% i% n
               select case (obs% o% t_int (spot%i%i+i))
               case (OBS_DUM)
                 k = i_d
               case (OBS_U)
                 k = i_u
               case (OBS_V)
                 k = i_v
               end select
               if (k > 0) then
                 obs% o% t_int (spot%i%i+k) = obs% o% t_int (spot%i%i+i)
                 obs% o% lev   (spot%i%i+k) = obs% o% lev   (spot%i%i+i)
               endif
             end do
             spot% i% n = ni
           endif
        else
          call new_int (obs% o, spot, ni)

          if (i_d > 0) then
            obs% o% t_int (spot%i%i+i_d) = OBS_DUM
            obs% o% lev   (spot%i%i+i_d) = log (100000._wp)
          endif

          if (i_u > 0) then
            obs% o% t_int (spot%i%i+i_u) = OBS_U
            obs% o% lev   (spot%i%i+i_u) = log (obs% o% olev (spot%o%i+o_v))
          endif
          if (i_v > 0) then
            obs% o% t_int (spot%i%i+i_v) = OBS_V
            obs% o% lev   (spot%i%i+i_v) = log (obs% o% olev (spot%o%i+o_v))
          endif
        endif
!     endif
      if (iand (TSK_SETUP_FUL0,tsk) /= 0) tsk = tsk - TSK_SETUP_FUL0
      if (iand (TSK_SHRINK    ,tsk) /= 0) tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !========================
    ! determine model columns
    !========================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then

      iatm = COL_UV

      call idx_init (       &
             spot% col% c,  &! <-  column descriptor
             spot% col% h,  &!  -> interpolation coefficients
             obs% o% mc,    &! <-> model column descriptors
             iatm,          &! <-  fields required
             0,             &! <-  tracers required
             atm% grid,     &! <-  model grid
             spot% i_time,  &! <-  time slot
             spot% w_time   )! <-  time interpolation weight

      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)

      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif

    !=================================================================
    ! load temp data from observation data type for further processing
    !=================================================================
!   n = spot% o% n
!!! call load_scatt (obs% o, spot, o)

    !===========================================
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
      if (obs% o% pe == dace% pe) then
        !=================================
        ! setup SCATT observational errors
        !=================================
        ty => obs% o% varno (spot% o% i + 1 : spot% o% i + spot% o% n)
        call scatt_obs_err (obs, spot, ty)
        !=================
        ! setup vqc bounds
        !=================
        call set_vqc_insitu (spot, obs% o)
      endif
      tsk = tsk - TSK_R
      if (tsk == 0) return
    endif

    !==============================
    ! tsk == TSK_Y
    ! evaluate observation operator
    !==============================
    if (iand (TSK_Y,tsk) /= 0) then
      iii = spot% i% i
      ioi = spot% o% i
      !------------
      ! set indices
      !------------
      i_u = 0; i_v = 0
      do i = 1, spot% i%n
        select case (obs% o% t_int (iii+i))
        case (OBS_U)
          i_u = i
        case (OBS_V)
          i_v = i
        case default
          call finish ('process_scatt(TSK_Y)','invalid t_int')
        end select
      end do

      !-------------------------------------------------
      ! scaling factor for near-surface wind observation
      !-------------------------------------------------
      sf = 1._wp
      !------------------------------
      ! evaluate observation operator
      !------------------------------
      obs% o% body (ioi+1:ioi+spot% o% n)% op_na = 0
      ff = -1._wp
      do i = 1, spot% o% n
        select case (obs% o% varno (ioi+i))
        case (VN_FF)
          if (ff < 0._wp) call fd_uv (ff, dd, xi%x (iii+i_u), xi%x (iii+i_v))
          y%x (ioi+i) = ff * sf
        case (VN_DD)
          if (ff < 0._wp) call fd_uv (ff, dd, xi%x (iii+i_u), xi%x (iii+i_v))
          y%x (ioi+i) = dd
        case (VN_U,  VN_U10M)
          y%x (ioi+i) =  xi%x (iii+i_u) * sf
        case (VN_V,  VN_V10M)
          y%x (ioi+i) =  xi%x (iii+i_v) * sf
        case default
          call finish ('process_scatt(TSK_Y)',                        &
                       'invalid varno: '//char3(obs% o% varno (ioi+i)))
        end select
      end do

      tsk = tsk - TSK_Y
      if (tsk == 0) return
    endif

    !================
    ! tsk == TSK_K
    ! derive Jakobian
    !================
    if (iand (TSK_K,tsk) /= 0) then
      ioi  = spot% o% i
      iii  = spot% i% i
      !------------
      ! set indices
      !------------
      i_u = 0; i_v = 0
      do i = 1, spot% i%n
        select case (obs% o% t_int (iii +i))
        case (OBS_U)
          i_u = i
        case (OBS_V)
          i_v = i
        case default
          call finish ('process_scatt(TSK_K)','invalid t_int')
        end select
      end do

      o_u=0; o_v=0; o_f=0; o_d=0
      do i = 1, spot% o%n
        select case (obs% o% varno (ioi+i))
        case (VN_U,  VN_U10M)
          o_u  = i
        case (VN_V,  VN_V10M)
          o_v  = i
        case (VN_FF)
          o_f  = i
        case (VN_DD)
          o_d  = i
        end select
      end do

      !------------------------------------
      ! derivatives d ff / d u , d ff / d v
      !------------------------------------
      if (o_f > 0 .or. o_d > 0) then
        call fd_uv (ff, dd, obs% xi% x(iii+i_u), obs% xi% x(iii+i_v), &
                    ff_x, ff_y, dd_x, dd_y                            )
      endif
      !-------------------------------------------------
      ! scaling factor for near-surface wind observation
      !-------------------------------------------------
      sf = 1._wp
      !-------------
      ! set Jacobian
      !-------------
      k = obs% H% ia (iii + 1)
!NEC$ nomove
      do i = 1, spot% i%n
        select case (obs% o% t_int (iii+i))
        case (OBS_U)
          if (o_u > 0) then
            obs% H% ja     (k)     = ioi+o_u        ! row    index
            obs% H% packed (k)     = real (sf,mp)   ! coefficient
            k = k + 1
          endif
          if (o_f > 0) then
            obs% H% ja     (k)     = ioi+o_f        ! row    index
            obs% H% packed (k)     = ff_x * sf      ! coefficient
            k = k + 1
          endif
          if (o_d > 0) then
            obs% H% ja     (k)     = ioi+o_d        ! row    index
            obs% H% packed (k)     = dd_x           ! coefficient
            k = k + 1
          endif
        case (OBS_V)
          if (o_v > 0) then
            obs% H% ja     (k)     = ioi+o_v        ! row    index
            obs% H% packed (k)     = real (sf,mp)   ! coefficient
            k = k + 1
          endif
          if (o_f > 0) then
            obs% H% ja     (k)     = ioi+o_f        ! row    index
            obs% H% packed (k)     = ff_y * sf      ! coefficient
            k = k + 1
          endif
          if (o_d > 0) then
            obs% H% ja     (k)     = ioi+o_d        ! row    index
            obs% H% packed (k)     = dd_y           ! coefficient
            k = k + 1
          endif
        case default
          call finish ('process_scatt(TSK_K)','invalid t_int')
        end select
        obs% H% ia (iii + i + 1) = k
      end do

      !------------------------------
      ! evaluate observation operator
      !------------------------------
      do i = 1, spot% o% n
        select case (obs% o% varno (ioi+i))
        case (VN_FF)
          obs%yi%x (ioi+i) = ff * sf
        case (VN_DD)
          obs%yi%x (ioi+i) = dd
        case (VN_U,  VN_U10M)
          obs%yi%x (ioi+i) = obs%xi%x (iii+i_u) * sf
        case (VN_V,  VN_V10M)
          obs%yi%x (ioi+i) = obs%xi%x (iii+i_v) * sf
        case default
          call finish ('process_scatt(TSK_K)',                        &
                       'invalid varno: '//char3(obs% o% varno (ioi+i)))
        end select
      end do

      tsk = tsk - TSK_K
      if (tsk == 0) return
    endif

    !===========
    ! left tasks
    !===========
    if (tsk /= 0) then
      if (dace% lpio) write (6,*) 'process_scatt: unknown task',tsk
      call finish ('process_scatt','unknown task')
    endif
  end subroutine process_scatt

!==============================================================================

  subroutine read_scatt_bufr (bufr, spt, obs, lkeep, cc)

  type (t_bufr) ,intent(inout)        :: bufr  ! BUFR record to decode
  type (t_spot) ,intent(inout)        :: spt   ! meta information to set
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type to set
  logical       ,intent(out)          :: lkeep ! accept observation ?
  integer       ,intent(in) ,optional :: cc    ! part of year ccyy

    !----------------
    ! index variables
    !----------------
    integer        :: is               ! sub-set    index
    integer        :: ie               ! entry in sub-set
    integer        :: id               ! descriptor index

    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
    integer            :: ival  ! value decoded from BUFR (integer)
    real(sp)           :: rval  ! value decoded from BUFR (real)
    character(len=8)   :: ymnem ! value decoded from BUFR (string)

    type(t_scatt)      :: s                ! SCATT observation read
    integer            :: yyyy,mo,dd,hh,mi ! actual time read
    logical            :: unused = .false.

    !---------------------
    ! loop over data, copy
    !---------------------
    lkeep = .false.
    if (bufr% sec3% num_subsets /= 1) then
      call decr_rpt_use (spt, CHK_NOIMPL, &
                         comment='read_scatt_bufr: cannot handle subsets')
      return
    endif
    if (bufr% sec0% edition == 4) then
      spt% corme = bufr% sec1% update      ! Correction message?
    end if
    is = 1
    call construct_scatt (s)
    do ie=1,bufr% nbufdat(is)
      !-------------
      ! decode datum
      !-------------
      ival  = bufr% ibufdat (is,ie)
      !----------------------------------------
      ! if valid, copy datum from BUFR record
      !   to structures of type t_spot, t_datum
      !----------------------------------------
      if(ival /= inv_bufr) then
        !---------------------------------------
        ! decode real/character datum, mnemonics
        !---------------------------------------
        id    = bufr% idescidx(is,ie)
!       itype = bufr% itype(id)
        ymnem = bufr% ymnem(id)
        rval  = rvind
        if (bufr% is_char(id) == 0) then
          IF (ival /= inv_bufr) &
            rval = ival * bufr% scale(id)
        endif
        !----------------------
        ! copy individual datum
        !----------------------
        select case (ymnem)     ! BUFR data descriptors

        case ('MCORME') ! COR MESSAGE MARK
          spt% corme = ival
        !------------------------
        ! class 01 identification
        !------------------------
        case ('YDDDD') ! SHIP OR MOBILE LAND STATION IDENTIFIER
          call bufr_get_character (bufr, ival, spt% statid)
        case ('YCCCC') ! ICAO LOCATION INDICATOR
          call bufr_get_character (bufr, ival, spt% statid)
        case ('MABNN') ! BUOY/PLATFORM IDENTIFIER
          if(spt% hd% obstype == OT_DRIBU) then
            write(spt%  statid,'(i5.5)')       ival
          else
            write(spt%  statid,'("BP_",i5.5)') ival
          endif
        case ('YSSOSN') ! SHORT STATION OR SITE NAME
          if (spt% statid=='') call bufr_get_character (bufr,ival,spt% statid)
        case ('MDS') ! direction of motion
        case ('NVS') ! speed of motion [m/s]
        case ('MII') ! WMO   block number
          spt%  ident  = 1000      * int(rval)
        case ('NIII') ! WMO station number
          spt%  ident  = spt% ident + int(rval)
        !-------------------------
        ! class 02 instrumentation
        !-------------------------
        case ('NIX') ! type of station
        !-------------------------
        ! class 04 location (time)
        !-------------------------
        case ('MJJJ') ! year
          yyyy = ival
          if(yyyy<100 .and. present(cc)) yyyy = yyyy + cc * 100
        case ('MMM') ! month
          mo   = ival
        case ('MYYQ')   ! Q-bits for following value              [code table]
        case ('MYY') ! day
          dd   = ival
        case ('MGG') ! hour
          hh   = ival
        case ('NGG') ! minute
          mi   = ival
        !---------------------------------
        ! class 05 location (horizontal-1)
        !---------------------------------
        case ('MLAH')  ! latitude (high accuracy)                [deg]
          spt% col% c% dlat = rval
        case ('MLALAQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MLALA') ! latitude (coarse accuracy)              [deg]
          if (spt% col% c% dlat == rvind) spt% col% c% dlat = rval
        !---------------------------------
        ! class 06 location (horizontal-2)
        !---------------------------------
        case ('MLOH')  ! longitude (high accuracy)                [deg]
          spt% col% c% dlon = rval
        case ('MLOLOQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MLOLO') ! longitude (coarse accuracy)              [deg]
          if (spt% col% c% dlon == rvind) spt% col% c% dlon =rval
        !-----------------------------
        ! class 07 location (vertical)
        !-----------------------------
        case ('MHP') ! height of station                    [m]
!!!       call set_datum (s% zr ,rval ,spt% corme)
        case ('MPN') ! pressure (vert.location)             [Pa]
!!!       call set_datum (s% p ,rval ,spt% corme)
        case ('MPNQ')! Q-bits for following value           [code table]
!!!       call set_qbits (s% p ,ival)
        !---------------------------------
        ! class 08 significance qualifiers
        !---------------------------------
        !----------------------------------------
        ! class 10 vertical elements and pressure
        !----------------------------------------
        case ('MHHA')  ! height of land surface                  [m]
!!!       call set_datum (s% z_ls ,rval ,spt% corme)
        case ('MPPP')  ! pressure                                [Pa]
!!!       call set_datum (s% ps ,rval ,spt% corme)
        case ('MPPPQ') ! Q-bits for following value              [code table]
!!!       call set_qbits (s% ps ,ival)
        case ('MPPPP') ! pressure reduced to mean sea level      [Pa]
!!!       call set_datum (s% p_msl ,rval ,spt% corme)
        case ('MPPPPQ')! Q-bits for following value              [code table]
!!!       call set_qbits (s% p_msl ,ival)
        case ('MPHPH') ! altimeter setting (QNH) (METAR)         [Pa]
!!!       call set_datum (s% p_msl ,rval ,spt% corme)
        case ('NHNHN') ! geopotential                            [m**2/s**2]
!!!       call set_datum (s% gp    ,rval ,spt% corme)
!!!       s% gp% o = s% gp% o / gacc                        ! -> [gpm]
        case ('NHNHNQ')
!!!       call set_qbits (s% gp    ,ival)
        !-----------------------------
        ! class 11 wind and turbulence
        !-----------------------------
        case ('NDD') ! wind direction        at 10 m   [deg]!->[rad]
          call set_datum (s% dd ,rval ,spt% corme) ! * d2r
        case ('NDDQ')
          call set_qbits (s% dd ,ival)
        case ('NFF') ! wind speed            at 10 m        [m/s]
          call set_datum (s% ff ,rval ,spt% corme)
        case ('NFFQ')
          call set_qbits (s% ff ,ival)
        !---------------------
        ! class 12 temperature
        !---------------------
        case ('MTTT')  ! dry bulb  temperature at  2 m        [K]
!!!       call set_datum (s% t ,rval ,spt% corme)
        case ('MTTTQ')
!!!       call set_qbits (s% t ,ival)
        case ('MTFTF') ! wet bulb  temperature at  2 m        [K]
        case ('MTDTD') ! dew point temperature at  2 m        [K]
!!!       call set_datum (s% td ,rval ,spt% corme)
        case ('MTDTDQ')
!!!       call set_qbits (s% td ,ival)
        !--------------------------------------------
        ! class 13 hydrographic/hydrological elements
        !--------------------------------------------
        case ('MUUU')  ! relative humidity     at  2 m   [%]->[ ]
!!!       call set_datum (s% rh ,rval * 0.01_sp ,spt% corme)
        case ('MUUUQ')
!!!       call set_qbits (s% rh ,ival)
        !----------------
        ! BUFR4 templates
        !----------------
        case ('MHOSNN') ! height of station ground above mean sea  [M]
!!!       call set_datum (s% z_ls ,rval ,spt% corme)
        case ('MHOBNN') ! height of barometer above mean sea level [M]
!!!       call set_datum (s% zr ,rval ,spt% corme)
        case ('MTDBTQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% t ,ival)
        case ('MTDBT')  ! temperature/dry bulb temperature         [K]
!!!       call set_datum (s% t ,rval ,spt% corme)
        case ('NDNDNQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% dd ,ival)
        case ('NDNDN')  ! wind direction                           [DEGREE_TRUE]
          call set_datum (s% dd ,rval ,spt% corme)
        case ('NFNFNQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% ff ,ival)
        case ('NFNFN')  ! wind speed                               [M/S]
          call set_datum (s% ff ,rval ,spt% corme)
        case ('MTDNHQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% td ,ival)
        case ('MTDNH')  ! dew-point temperature                    [K]
!!!       call set_datum (s% td ,rval ,spt% corme)
        case ('NHHHN')  ! geopotential height                      [GPM]
!!!       call set_datum (s% gp    ,rval ,spt% corme)
#ifdef  CHECKCODES
        !-----------------------
        ! check for unused codes
        !-----------------------
        case ('MOBITQ') ! overall quality bits                    [code table]
        case ('MADDF')  ! associated field significance           [code table]
        case ('MA1')    ! WMO region number                       [code table]
        case ('NIW')    ! type of instrument.for wind measurement [flag table]
        case ('MADDF0') ! associated field significance           [code table]
        case ('NPPPQ')  ! Q-bits for following value              [code table]
        case ('NPPP')   ! 3 hour pressure change                  [Pa]
        case ('NA')     ! characteristic of pressure tendency     [code table]
        case ('MVVQ')   ! Q-bits for following value              [code table]
        case ('MVV')    ! horizontal visibility                   [m]
        case ('NWWQ')   ! Q-bits for following value              [code table]
        case ('NWW')    ! present weather                         [code table]
        case ('MW1Q')   ! Q-bits for following value              [code table]
        case ('MW1')    ! past weather (1)                        [code table]
        case ('MW2Q')   ! Q-bits for following value              [code table]
        case ('MW2')    ! past weather (2)                        [code table]
        case ('MNQ')    ! Q-bits for following value              [code table]
        case ('MN')     ! cloud cover (total)                     [%]
        case ('NHQ')    ! Q-bits for following value              [code table]
        case ('NH')     ! height of base of cloud                 [m]
        case ('MCCQ')   ! Q-bits for following value              [code table]
        case ('MCC')    ! cloud type                              [code table]
        case ('MCC0Q')  ! Q-bits for following value              [code table]
        case ('MCC0')   ! cloud type                              [code table]
        case ('MCC1Q')  ! Q-bits for following value              [code table]
        case ('MCC1')   ! cloud type                              [code table]
        case ('NFXGU')  ! maximum wind speed(gusts)               [m/s]
        case ('MGGTR1') ! durat.of time relat.to following value  [hour]
        case ('MRRRQ')  ! Q-bits for following value              [code table]
        case ('MRRR')   ! total precipitation/total water equiv.  [kg/m**2]
        case ('MNHQ')   ! Q-bits for following value              [code table]
        case ('MNH')    ! cloud amount                            [code table]
        case ('MVTSU')  ! vertical significance (surface observ.) [code table]
        case ('MCC2')   ! cloud type                              [code table]
        case ('NH0Q')   ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('NH0')    ! height of base of cloud                 [m]
        case ('MNH0Q')  ! Q-bits for following value              [code table]
        case ('MNH0')   ! cloud amount                            [code table]
        case ('NSSS')   ! total snow depth                        [m]
        case ('MGGACT') ! actual hour of observation              [hour]
        case ('NGGACT') ! actual minute of observation            [minute]
        case ('MPW')    ! period of waves                         [s]
        case ('MHW')    ! height of waves                         [m]
        case ('MCOTR')  ! duration of time relat.to follow.value  [code table]
        case ('MGGTR3') ! durat.of time relat.to following value  [hour]
        case ('MSSSS')  ! global radiation integr.o.period specif.[J/m**2]
        case ('YSUPL')  ! 008 characters                          [CCITT IA5]
        case ('NFXME')  ! maximum wind speed(10 min mean wind)    [m/s]
        case ('MMOSTM') ! method of sea-surface temperature measu.[code table]
        case ('MTSQ')   ! Q-bits for following value              [code table]
        case ('MTS')    ! sea/water temperature                   [K]
        case ('MGGTR2') ! durat.of time relat.to following value  [hour]
        case ('MSSSM0') ! total sunshine                          [minute]
        case ('ME')     ! state of ground (with or without snow)  [code table]
        case ('MGGTR4') ! durat.of time relat.to following value  [hour]
        case ('MDFSR')  ! diff.solar radiation integr.o.per.spec. [J/m**2]
        case ('MPWPW')  ! period of wind waves                    [s]
        case ('MHWHW')  ! height of wind waves                    [m]
        case ('MGGTR5') ! durat.of time relat.to following value  [hour]
        case ('MLWR')   ! long wave radiat.,integr./period specif.[J/m**2]
        case ('ME0')    ! state of ground (with or without snow)  [code table]
        case ('MGGTR0Q')! Q-bits for following value              [code table]
        case ('MGGTR0') ! durat.of time relat.to following value  [hour]
        case ('MTNTNQ') ! Q-bits for following value              [code table]
        case ('MTNTN')  ! minimum temp., height/period specified  [K]
        case ('MSSSMQ') ! Q-bits for following value              [code table]
        case ('MSSSM')  ! total sunshine                          [minute]
        case ('MTXTXQ') ! Q-bits for following value              [code table]
        case ('MTXTX')  ! maximum temp., height/period specified  [K]
        case ('MGGTRQ') ! Q-bits for following value              [code table]
        case ('MGGTR')  ! durat.of time relat.to following value  [hour]
        case ('MTGTG')  ! ground minimum temp., past 12 hours     [K]
        case ('MVTSU0') ! vertical significance (surface observ.) [code table]
        case ('MNH1')   ! cloud amount                            [code table]
        case ('MCC3')   ! cloud type                              [code table]
        case ('MCT')    ! cloud top description                   [code table]
        case ('NDW12')  ! direction of swell waves                [degree true]
        case ('NHT')    ! height of top of cloud                  [m]
        case ('NP24Q')  ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('NP24')   ! 24 hour pressure change                 [Pa]
        case ('MEEE')   ! evaporation/evapotranspiration          [kg/m**2]
        case ('NIE')    ! instrument/crop for evapo(transpi)ration[code table]
        case ('MPW12')  ! period of swell waves                   [s]
        case ('MHW12')  ! height of swell waves                   [m]
        case ('MMOWTM') ! method of wet-bulb temperature measurem.[code table]
        case ('Loop000':'Loop999')  ! Start of Loop
        case ('Lcnt000':'Lcnt999')  ! Loop Counter
        case ('MDREP')  ! delayed descriptor replication factor      numeric
        !------
        ! METAR
        !------
        case ('NWWIP')  ! intensity or proximity of weather       [CODE_TABLE]
        case ('YWWP')   ! weather phenomena (METAR)               [CCITT_IA5]
        case ('NWWIP0') ! intensity or proximity of weather       [CODE_TABLE]
        case ('YWWP0')  ! weather phenomena (METAR)               [CCITT_IA5]
        case ('YWWD')   ! descriptor of weather                   [CCITT_IA5]
        case ('YWWD0')  ! descriptor of weather                   [CCITT_IA5]
        !-----------------------
        ! BUFR4 template (DRIBU)
        !-----------------------
        case ('MWRSA')  ! WMO REGION SUB-AREA                      [NUMERIC]
        case ('MTODB')  ! TYPE OF DATA BUOY                        [CODE_TABLE]
        case ('MTISI')  ! TIME SIGNIFICANCE                        [CODE_TABLE]
        case ('MJJJ0')  ! YEAR                                     [YEAR]
        case ('MMM0')   ! MONTH                                    [MONTH]
        case ('MYY0')   ! DAY                                      [DAY]
        case ('MGG0')   ! HOUR                                     [HOUR]
        case ('NGG0')   ! MINUTE                                   [MINUTE]
        case ('MQBST')  ! QUALITY OF BUOY SATELLITE TRANSMISSION   [CODE_TABLE]
        case ('MQOBL')  ! QUALITY OF BUOY LOCATION                 [CODE_TABLE]
        case ('MLQC')   ! LOCATION QUALITY CLASS (..66 CONFIDENCE) [CODE_TABLE]
        case ('NDRTY')  ! DROGUE TYPE                              [CODE_TABLE]
        case ('MDCI')   ! DEPTH CORRECTION INDICATOR               [CODE_TABLE]
        case ('NZNZN')  ! DEPTH BELOW SEA/WATER SURFACE            [M]
        case ('MTN00Q') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MTN00')  ! SEA/WATER TEMPERATURE                    [K]
        case ('MDREP0') ! DELAYED DESCRIPTOR REPLICATION FACTOR    [NUMERIC]
        case ('MTISI1') ! TIME SIGNIFICANCE                        [CODE_TABLE]
        case ('NGGTP')  ! TIME PERIOD OR DISPLACEMENT              [MINUTE]
        case ('MTISI3') ! TIME SIGNIFICANCE                        [CODE_TABLE]
        case ('NOMDP')  ! OPERATOR OR MANUFACTURER DEF. PARAMETER  [NUMERIC]
        case ('MDREP1') ! DELAYED DESCRIPTOR REPLICATION FACTOR    [NUMERIC]
        case ('NDRDE')  ! DROGUE DEPTH                             [M]
        case ('NK2')    ! METHOD OF SALINITY/DEPTH MEASUREMENT     [CODE_TABLE]
        case ('MSNSN')  ! SALINITY                                 [CODE_TABLE]
        case ('MHAWAS1')! HEIGHT OF SENSOR ABOVE WATER SURFACE     [M]
        case ('MANTYP') ! ANEMOMETER TYPE                          [CODE_TABLE]
        !------------------------
        ! BUFR4 templates (SYNOP)
        !------------------------
        case ('MGGTP')  ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP1') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP2') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP3') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP4') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('NCI')    ! SEA ICE CONCENTRATION                    [CODE_TABLE]
        case ('MAMTI')  ! AMOUNT AND TYPE OF ICE                   [CODE_TABLE]
        case ('MICSI')  ! ICE SITUATION                            [CODE_TABLE]
        case ('MICDE')  ! ICE DEVELOPMENT                          [CODE_TABLE]
        case ('MSI')    ! ICE DEPOSIT (THICKNESS)                  [M]
        case ('MRS')    ! RATE OF ICE ACCREATION                   [CODE_TABLE]
        case ('NGGTP0') ! TIME PERIOD OR DISPLACEMENT              [MINUTE]
        case ('MTNTNHQ')! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MTNTNH') ! MINIMUM TEMP., HEIGHT/PERIOD SPECIFIED   [K]
        case ('YSOSN')  ! STATION OR SITE NAME                     [CCITT_IA5]
        case ('MHOSEN') ! HEIGHT OF SENSOR ABOVE LOCAL GROUND      [M]
        case ('MVTSU2') ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)  [CODE_TABLE]
        case ('MGGTP5') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('NGGTP1') ! TIME PERIOD OR DISPLACEMENT              [MINUTE]
        case ('MGGTP6') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP0') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP7') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGLSR')  ! GLOBAL SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]
        case ('MTXTXHQ')! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MTXTXH') ! MAXIMUM TEMP., HEIGHT/PERIOD SPECIFIED   [K]
        case ('MRR24Q') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MRR24')  ! TOTAL PRECIPITATION, PAST 24 HOURS       [KG/M**2]
        case ('MDSRH')  ! DIFFUSE SOLAR RADIAT. (HIGH ACC.) INTEGR [J/M**2]
        case ('MEEEV')  ! EVAPORATION/EVAPOTRANSPIRATION           [KG/M**2]
        case ('MTGTGH') ! GROUND MINIMUM TEMP., PAST 12 HOURS      [K]
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
      end if
    end do
    !--------------------------------------
    ! in case of unknown codes:
    ! print suspicious BUFR record and exit
    !--------------------------------------
    if (unused) then
      call bufr_print_sections (bufr)
      call bufr_print_subset   (bufr, is)
      call finish ('read_scatt_bufr','code(s) not implemented')
    endif
    !-----------------------------------------------
    ! if station name is not set, use station number
    !-----------------------------------------------
    if (spt% statid   =='')   write(spt% statid,'(i5.5)') spt% ident
    !----------------------------------------
    ! 'QSCAT', 'ASCAT', ... for Scatterometer
    !----------------------------------------
    s% z = s% z_msl             ! Preset used height from reported height
    select case (spt% hd% dbkz)
    case (1697)
      spt% statid = 'QSCAT'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1698,1699)
      spt% statid = 'ASCAT'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1700)
      spt% statid = 'OSCAT'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1701)
      spt% statid = 'JASON-2'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1702)
      spt% statid = 'SARAL'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1780)
      spt% statid = 'SENTINEL-3'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    end select
    !----------------
    ! standard checks
    !----------------
    lkeep = .true.
    call init_time (spt% actual_time, yyyy, mo, dd, hh, mi)
    call check_report_1 (spt)
    if (spt% use% state <= STAT_DISMISS) lkeep = .false.

    if (lkeep) call check_store_scatt (s, spt, obs, lkeep)

    if (prt_data) call print_scatt (spt, s)

  end subroutine read_scatt_bufr
!------------------------------------------------------------------------------
  subroutine check_store_scatt (scatt, spot, obs, lkeep, repl)
  type(t_scatt),intent(inout)        :: scatt ! SCATT level information
  type(t_spot) ,intent(inout)        :: spot  ! meta data of this observation
  type(t_obs)  ,intent(inout)        :: obs   ! data of all observations
  logical      ,intent(out)          :: lkeep ! flag: keep or reject
  integer      ,intent(in) ,optional :: repl  ! observation to replace

    !----------------
    ! local variables
    !----------------
    target                  :: scatt
    integer, parameter      :: mv = 4     ! max. no. variables used
    integer                 :: typ (mv)   ! observation type
    real(wp)                :: lev (mv)   ! level of observation
    type(t_datum)           :: bod (mv)   ! body entry
    integer                 :: id
    type (t_spot)  ,pointer :: s
    integer                 :: n          ! number of valid observations
    integer                 :: i
    type (t_datum) ,pointer :: zs         ! link to scatt% z_msl
!   integer                 :: zs_to_z    ! copy of zls_to_z    or zmsl_to_z
    logical                 :: scat

    n     = 0
    lkeep = .false.
    scat  = .false.

    !----------------
    ! set preferences
    !----------------
    select case (spot% hd% buf_type)
    case (12)
!     zs_to_z =  zmsl_to_z    ! station height to use
      zs      => scatt% z_msl
    case default
      call finish('check_store_scatt','BUFR-type not 12')
    end select
    scat = any (spot% hd% dbkz == kz_scat)

    !--------------------------
    ! always use surface height
    !--------------------------
    scatt% z = zs

    !--------------------------------
    ! store station or surface height
    !--------------------------------
    if (scatt% z% qc == 0) spot% z = scatt% z% o

    if(prt_data) print *, repeat('=',64)

    !-------------------------------------------
    ! select reference pressure for ff, dd, u, v
    !-------------------------------------------
    scatt% p_ref % o = use_ps_model

    !--------------------------------------------
    ! insert in list if information is sufficient
    !--------------------------------------------
    if (spot% col% c% dlon     /= invalid .and. &
        spot% col% c% dlat     /= invalid       ) then
      !------------------
      ! use wind over sea
      !------------------
      if (scatt% ff% qc == 0 .and. scatt% dd% qc == 0) then
        scatt% uu% src = SRC_DER
        scatt% uu% qc  = QC_OK
        scatt% uu% o   = scatt% ff% o * (-sin (d2r * scatt% dd% o))
        scatt% vv% src = SRC_DER
        scatt% vv% qc  = QC_OK
        scatt% vv% o   = scatt% ff% o * (-cos (d2r * scatt% dd% o))

        n = n + 1
        if(prt_data) print *,'u : n =',n
        typ (n)             = VN_U
        bod (n)             = scatt% uu
        bod (n)% plev       = scatt% p_ref% o
        lev (n)             = scatt% p_ref% o
        n = n + 1
        if(prt_data) print *,'v : n =',n
        typ (n)             = VN_V
        bod (n)             = scatt% vv
        bod (n)% plev       = scatt% p_ref% o
        lev (n)             = scatt% p_ref% o
      end if

      if (scatt% ff% qc == 0 .and.                                         &
          ((scat   .and. scatt% dd% qc /= 0) .or.                          &
           (monitor_ff .and. scatt% ff% qc == 0 .and. scatt% dd% qc == 0)))&
        then
        !----------------------------------------------------------
        ! Scatterometer-like wind speed data without wind direction
        ! or FF monitoring
        !----------------------------------------------------------
        n = n + 1
        if (prt_data) print *,'ff: n =',n
        typ (n)             = VN_FF
        bod (n)             = scatt% ff
        bod (n)% plev       = scatt% p_ref% o
        lev (n)             = scatt% p_ref% o
      end if

      if (scatt% ff% qc == 0 .and. scatt% ff% o > 0. .and. &
          scatt% dd% qc == 0 .and. monitor_dd              ) then
        n = n + 1
        if (prt_data) print *,'dd: n =',n
        typ (n)             = VN_DD
        bod (n)             = scatt% dd
        bod (n)% plev       = scatt% p_ref% o
        lev (n)             = scatt% p_ref% o
      end if
    end if

    !-----------------------------------------------------------------
    ! Catch forgotten adjustment of mv if compiled w/o bounds checking
    !-----------------------------------------------------------------
    if (n > mv) call finish('check_store_scatt','increase mv')

!   if (prt_data) call print_scatt (spot, scatt)

    !-----------------------------------------
    ! Drop bogus reports at exactly (0°N,0°E).
    !-----------------------------------------
    if (spot% col% c% dlat == 0._wp .and. &
        spot% col% c% dlon == 0._wp) then
       n = 0
    end if

    !--------------------
    ! store data into OBS
    !--------------------
    if (n>0) then
      lkeep = .true.
      if (present (repl)) then
        s => obs% spot (repl)
      else
        call new_spot (obs, 1, set_id=.true.)
        s => obs% spot (obs% n_spot)
      endif
      id = s% id
      s = spot
!     s% int_type  = ITY_ICOL
      s% id        = id
      s% col% nlev = 1
!     s% cost      = 1._wp
!     s% char      = CHR_ID
!     s% nr        = n
      call new_obs (obs, n, s)
!     call new_int (obs, s, n)
      spot% o      = s% o      ! return observation space info
      obs% varno (s%o%i+1:s%o%i+s%o%n) =      typ (1:n)
      obs% body  (s%o%i+1:s%o%i+s%o%n) =      bod (1:n)
!     obs% t_int (s%i%i+1:s%i%i+s%i%n) =      typ (1:n)
      obs%  olev (s%o%i+1:s%o%i+s%o%n) =      lev (1:n)
!     obs%   lev (s%i%i+1:s%i%i+s%i%n) = log (lev (1:n))
      call store_scatt (obs, s, scatt)
      call set_xuv (s)
      !--------------------------------------
      ! certain quantities are always passive
      !--------------------------------------
        do i = 1, s% o% n
          select case   (obs% varno(s% o% i+i))
          case (VN_FF)
            !----------------------------------------
            ! FF is active only for SCATT without u,v
            !----------------------------------------
            if (s% o% n > 1)                                           &
              call decr_use (obs% body (s% o% i+i)% use, STAT_PASSIVE, &
                             check=CHK_NOTUSED                         )
          case (VN_DD)
            call decr_use (obs% body (s% o% i+i)% use, STAT_PASSIVE, &
                           check=CHK_NOTUSED                         )
          end select
        end do
    else
      call decr_rpt_use (spot ,CHK_INSDAT, STAT_DISMISS)
    endif
  end subroutine check_store_scatt
!==============================================================================
  pure subroutine construct_scatt (s)
  type (t_scatt) ,intent(out) :: s

    s% total %mn =''
    s% z     %mn ='z'
    s% z_msl %mn ='zmsl'
    s% ff    %mn ='ff'
    s% dd    %mn ='dd'
    s% uu    %mn ='u'
    s% vv    %mn ='v'
    s% z_msl %o  = 0.
    s% z_msl %qc = 0_i2

  end subroutine construct_scatt
!==============================================================================
  subroutine print_scatt (s, d)
  type (t_spot ) ,intent(in) :: s
  type (t_scatt) ,intent(in) :: d

    write (6,*)
    write (6,*)' statid               = ',s% statid
    write (6,*)' lat lon              = ',s% col% c% dlat, s% col% c% dlon
    write (6,*)' buftype subtype dbkz = ',s% hd% buf_type,               &
                                          s% hd% buf_subtype, s% hd% dbkz
    write (6,*)' time actual db       = ',chhmm(s% hd% time) ,' ',       &
                                          chhmm(s% actual_time),' ',     &
                                          chhmm(s% hd% db_time)
    write (6,*)' file record          = ',s% hd% source, s% hd% record
    write (6,*)
    write (6,'(a4,4x,a15,3a4)')'name','value','qc','use','src'

    call print (d% total ,'quality control flags')
    call print (d% z_msl ,'height   mean sea level    [m]')
    call print (d% p_ref ,'reference pressure         [Pa]')
    call print (d% ff    ,'wind speed                 [m/s]')
    call print (d% dd    ,'wind direction             [deg]')
    call print (d% uu    ,'wind component             [m/s]')
    call print (d% vv    ,'wind component             [m/s]')
  end subroutine print_scatt
!==============================================================================
  subroutine load_scatt (obs, spot, scatt)
  type (t_obs)   ,intent(in)  :: obs   ! data of all observations
  type (t_spot)  ,intent(in)  :: spot  ! meta data of this observation
  type (t_scatt) ,intent(out) :: scatt ! SCATT information
  !------------------------------------------------------------------------
  ! Load the data from components PAR, OBS of OBS from position provided by
  ! SPOT. Store into SCATT.
  !------------------------------------------------------------------------
  scatt = transfer (obs% par (spot% p% i+1 : spot% p% i + spot% p% n), scatt)
  end subroutine load_scatt
!------------------------------------------------------------------------------
  subroutine store_scatt  (obs, spot, scatt)
  type (t_obs)   ,intent(inout) :: obs   ! data of all observations
  type (t_spot)  ,intent(inout) :: spot  ! meta data of this observation
  type (t_scatt) ,intent(in)    :: scatt ! SCATT level information
  !-----------------------------------------------------------------------
  ! Store the data from variable SCATT in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------
    integer :: n, par(1)  !+++ work around bug in NEC sxf90/2.0rev360-ia64
    if (scatt_int_size==0) scatt_int_size = size (transfer (inva_sca,(/0/)))
    n = spot% col% nlev * scatt_int_size
    call new_par (obs, n, spot=spot)
    obs % par (spot% p% i+1 : spot% p% i + spot% p% n) = &
      transfer(scatt, par)
  end subroutine store_scatt
!==============================================================================
  subroutine read_scatt_nml
  !-------------------------
  ! read namelist /SCATT_OBS/
  !-------------------------
    integer :: ierr
    !-------------
    ! set defaults
    !-------------
    prt_data   = .false.   ! print data
!   zls_to_z   = 3         ! set z from land surface height
!   zmsl_to_z  = 4         ! set z from mean sea level height
    scatt_proc = 0         ! processing of scatterometer observations
    verbose    = 0         ! Verbosity level of consistency checks

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,*) 'Namelist /SCATT_OBS/:'
      write(6,'()')
      call position_nml ('SCATT_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=SCATT_OBS, iostat=ierr)
        if (ierr/=0) call finish ('read_scatt_nml',              &
                                  'ERROR in namelist /SCATT_OBS/')
#else
        read (nnml ,nml=SCATT_OBS)
#endif
      case default
        write(6,*) 'Namelist not present, defaults used'
        write(6,'()')
      end select
      !---------
      ! printout
      !---------
!     write(6,'(a,i9,a)')    'zls_to_z     =',zls_to_z    ,'    set z from land surface height'
!     write(6,'(a,i9,a)')    'zmsl_to_z    =',zmsl_to_z   ,'    set z from mean sea level height'
      write(6,'(a,i9,a)')    'scatt_proc   =',scatt_proc  ,'    processing of scatterometer observations'
      write(6,'(a,l9,a)')    'prt_data     =',prt_data    ,'    print data'
      write(6,'(a,i9,a)')    'verbose      =',verbose     ,'    verbosity'
    endif
!   call p_bcast (zls_to_z,    dace% pio)
!   call p_bcast (zmsl_to_z,   dace% pio)
    call p_bcast (scatt_proc,  dace% pio)
    call p_bcast (prt_data,    dace% pio)
    call p_bcast (verbose,     dace% pio)
!   if (dace% pe /= dace% pio) prt_data = .false.
  end subroutine read_scatt_nml

!==============================================================================
  subroutine read_scatt_netcdf (ifile, i_source, obs, rec1, recl, &
                                lkeep, nkeep)
  !====================================================================
  ! Read SCATT / surface observations from netCDF (converted BUFR data)
  !====================================================================
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in source-file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data type to set
  integer       ,intent(in)           :: rec1      ! first record to read
  integer       ,intent(in)           :: recl      ! last  record to read
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.

  type (t_spot)        :: spt0, spti   ! observation meta data
  type (t_spot) ,save  :: empty        !
  type (t_use)         :: use          ! status variable
  type (t_head)        :: head         ! meta information

  integer              :: bufr_type    ! BUFR message type    read
  integer              :: bufr_subtype ! BUFR message subtype read
  integer              :: centre       ! generating centre
  integer              :: ierm         ! error code

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
  integer              :: dimid_report   !                     for  Reports
  integer              :: dimid_l_fxg    !                     for  gusts
! integer              :: dimid_l_nff    !                     for  scat windspeeds

  integer              :: len_layer      ! number of vertical layers in NetCDF-File 1.type of field
  integer              :: len_report     ! number of reports in NetCDF file         1.type of field
  integer              :: len_l_fxg      ! number of gusts           in NetCDF-File 3.type of field
! integer              :: len_l_nff      ! number of scat windspeeds in NetCDF-File 3.type of field
  integer              :: len_report1    ! number of reports in NetCDF file         3.type of field
  integer              :: j              ! loop index
  integer              :: nc1            ! first  dimension for netcdf 1.type of field getvar in start / count
  integer              :: nc2            ! second dimension for netcdf 1.type of field getvar in start / count
  integer              :: nc3            ! first  dimension for netcdf 2.type of field getvar in start / count
  integer              :: nc4            ! second dimension for netcdf 2.type of field getvar in start / count
  integer              :: nc5            ! first  dimension for netcdf 3.type of field getvar in start / count
  integer              :: nc6            ! second dimension for netcdf 3.type of field getvar in start / count
  integer              :: entry1,entry   ! position in source file (subset)

  integer ,allocatable :: ifd_1d  (:)    !
  real    ,allocatable :: rfd_1d  (:)    !
  integer ,allocatable :: ifd_2d1 (:,:)  !
  real    ,allocatable :: rfd_2d1 (:,:)  !
  integer ,allocatable :: ifd_2d2 (:,:)  !
  real    ,allocatable :: rfd_2d2 (:,:)  !

! dummy field
  integer ,allocatable :: ifd_1dm (:)    !
  real    ,allocatable :: rfd_1dm (:)    !
  integer ,allocatable :: ifd_2dm (:,:)  !
  real    ,allocatable :: rfd_2dm (:,:)  !

  !----------------
  ! index variables
  !----------------
  integer       :: nreport          ! number of observations (default)

  integer       :: is               ! sub-set    index
  integer       :: i, jj            !, index180

  logical       :: lpr_scatt        ! extended  printing of scatts
  logical       :: lpr_extd         ! extended  printing of scatts
  integer       :: npr_extd         ! number of extended printing of aireps
  logical       :: l_scat           ! report is scatterometer
  logical       :: l_altim          ! report is scatterometer-like altimeter

! logical for meteorological variables(over all reports)
! logical       :: l_mpppp
  logical       :: l_nff, l_ndd

! logical for buoys (over all reports)
! logical       :: l_mdrep

! logical for special scatterometer fields (over all reports)
  logical       :: l_norbnu, l_nctcn , l_nwvcq , l_mswcq , l_nnova , l_niswv , &
                   l_mpror , l_ntnsm,  l_mlcfs,  l_nrsst , l_nadqf , l_nacqf , &
                   l_narf  , l_nipi ,  l_nfffa,  l_nkswh,  l_nrasst, l_nrms20, &
                   l_nr20sw, l_nsbswh, l_nrmsor, l_nrmswh, l_nodle

  character(NF90_MAX_NAME)   :: yname_v  ! NetCDF variable name

  integer       :: i_dft           ! integer default value
  real          :: r_dft           ! real    default value
  real          :: zdir1, zdir2, zdel, zlargest

  !-----------------------------------------
  ! quantities derived from the BUFR message
  !-----------------------------------------
! real    ,allocatable      :: mpppp   (:)    ! PRESSURE REDUCED TO MSL [Pa]
  !---------
  ! loop 000
  !---------
  !---------
  ! loop 001
  !---------
! integer ,allocatable      :: mdrep   (:)    ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !--------------
  ! scatterometer
  !--------------
  integer ,allocatable      :: norbnu  (:)    ! ORBIT NUMBER
  integer ,allocatable      :: nctcn   (:)    ! CROSS-TRACK CELL NUMBER
! real,    allocatable      :: nssm    (:)    ! SURFACE SOIL MOISTURE (MS) [%]
! real,    allocatable      :: neesm   (:)    ! ESTIMATED ERROR IN SURFACE SOIL MOISTURE [%]
! real,    allocatable      :: nmssm   (:)    ! MEAN SURFACE SOIL MOISTURE
! real,    allocatable      :: nsmq    (:)    ! SOIL MOISTURE QUALITY [%]
! integer, allocatable      :: nsnco   (:)    ! SNOW COVER [%]
! real,    allocatable      :: nflsf   (:)    ! FROZEN LAND SURFACE FRACTION [%]
  integer ,allocatable      :: nwvcq   (:)    ! ASCAT WIND VECTOR CELL QUALITY
  integer ,allocatable      :: mswcq   (:)    ! SEAWIND WIND VECTOR CELL QUALITY
  real    ,allocatable      :: mpror   (:)    ! PROBABILITY OF RAIN
  integer ,allocatable      :: nnova   (:)    ! NUMBER OF VECTOR AMBIGUITIES
  integer ,allocatable      :: niswv   (:)    ! INDEX OF SELECTED WIND VECTOR
  integer ,allocatable      :: ntnsm   (:)    ! TOTAL NUMBER OF SIGMA-0 MEASUREMENTS
! integer ,allocatable      :: mdreps  (:)    ! DELAYED REPLICATION FACTOR
  real    ,allocatable      :: nff     (:,:)  ! WIND SPEED AT 10 M
  real    ,allocatable      :: ndd     (:,:)  ! WIND DIRECTION AT 10 M
  real    ,allocatable      :: mlcfs   (:,:)  ! LIKELIHOOD OF COMP. FOR SOLUTION
  integer ,allocatable      :: mle_flag(:)    ! QUALITY CONTROL FLAG DER SCATTEROMETERDATEN
  !----------
  ! altimeter
  !----------
  integer ,allocatable      :: nrsst   (:)    ! REMOTELY SENSED SURFACE TYPE
  integer ,allocatable      :: nrasst  (:)    ! RADIOMETER SENSED SURFACE TYPE
  integer ,allocatable      :: nadqf   (:)    ! ALTIMETER DATA QUALITY FLAG
  integer ,allocatable      :: nacqf   (:)    ! ALTIMETER CORRECTION QUALITY FLAG
  integer ,allocatable      :: narf    (:)    ! ALTIMETER RAIN FLAG
  integer ,allocatable      :: nipi    (:)    ! ICE PRESENCE INDICATOR
  real    ,allocatable      :: nfffa   (:)    ! WIND SPEED FROM ALTIMETER
  real    ,allocatable      :: nodle   (:)    ! OCEAN DEPTH/LAND ELEVATION
  !---------------
  ! Jason template
  !---------------
  real    ,allocatable      :: nrms20  (:)    ! RMS OF 20 HZ KU BAND OCEAN RANGE
  real    ,allocatable      :: nkswh   (:)    ! KU BAND SIGNIFICANT WAVE HEIGHT
  real    ,allocatable      :: nr20sw  (:)    ! RMS 20 HZ KU BAND SIGNIFICANT WAVE HEIGHT
  !---------------
  ! SARAL template
  !---------------
  real    ,allocatable      :: nrmsor  (:)    ! RMS OF SPECIFIC BAND OCEAN RANGE
  real    ,allocatable      :: nsbswh  (:)    ! SPECIFIC BAND SIGNIFICANT WAVE HEIGHT
  real    ,allocatable      :: nrmswh  (:)    ! RMS SPECIFIC BAND SIGNIFICANT WAVE HEIGHT
!
! not expanded variables:
! d  defined values
!  * in field list above integrated
!                NBVLR:long_name = "BATTERY VOLTAGE (LARGE RANGE)" ;
!                MTYOE0:long_name = "TYPE OF EQUIPMENT" ;
!                NBVLR0:long_name = "BATTERY VOLTAGE (LARGE RANGE)" ;
!                MTYOE1:long_name = "TYPE OF EQUIPMENT" ;
!                NBVLR1:long_name = "BATTERY VOLTAGE (LARGE RANGE)" ;
!                MTYOE2:long_name = "TYPE OF EQUIPMENT" ;
! d              NDRTY:long_name = "DROGUE TYPE" ;
!                MLDDS:long_name = "LAGRANGIAN DRIFTER DROGUE STATUS" ;
! d              NDRDE:long_name = "DROGUE DEPTH" ;
!                NLDS:long_name = "LAGRANGIAN DRIFTER SUBMERGENCE" ;
! d              MDCI:long_name = "DEPTH CORRECTION INDICATOR" ;
!                NCABL:long_name = "CABLE LENGTH" ;
!                MHPEC:long_name = "HYDROSTATIC PRESS. OF LOWER END OF CABLE" ;
!                MSI:long_name = "ICE DEPOSIT (THICKNESS)" ;
!                MMOSTM:long_name = "METHOD OF SEA-SURFACE TEMPERATURE MEASU." ;
!                NK1:long_name = "INDICATOR FOR DIGITIZATION" ;
! d              NK2:long_name = "METHOD OF SALINITY/DEPTH MEASUREMENT" ;
!                MDREP:long_name = "DELAYED DESCRIPTOR REPLICATION FACTOR" ;
! d              NZNZN:long_name = "DEPTH BELOW SEA/WATER SURFACE" ;
! d*             MTN00:long_name = "SEA/WATER TEMPERATURE" ;
!                MSNSN:long_name = "SALINITY" ;
!                MMOCM:long_name = "METHOD OF CURRENT MEASUREMENT" ;
!                NK3K4:long_name = "DURATION AND TIME OF CURRENT MEASUREMENT" ;
!
! d              MDREP0:long_name = "DELAYED DESCRIPTOR REPLICATION FACTOR" ;
!                MHOBNN:long_name = "HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL" ;
!                MTYOE3:long_name = "TYPE OF EQUIPMENT" ;
!                MTINST:long_name = "INSTRUMENT TEMPERATURE" ;
! d*             MPPP:long_name = "PRESSURE" ;
! d*             MPPPP:long_name = "PRESSURE REDUCED TO MSL" ;
! d*             NPPP:long_name = "3 HOUR PRESSURE CHANGE" ;
! d*             NA:long_name = "CHARACTERISTIC OF PRESSURE TENDENCY" ;
!                MTYOE4:long_name = "TYPE OF EQUIPMENT" ;
!                MHOSEN:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
!                MHAWAS:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
! d*             MTDBT:long_name = "TEMPERATURE/DRY BULB TEMPERATURE" ;
!  *             MTDNH:long_name = "DEW-POINT TEMPERATURE" ;
! d*             MUUU:long_name = "RELATIVE HUMIDITY" ;
!                MHOSEN0:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
!                MHAWAS0:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
!                NACH2V:long_name = "ART. CORREC. OF SENSOR H 2 ANOTHER VALUE" ;
! d              MHAWAS1:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
! d              MANTYP:long_name = "ANEMOMETER TYPE" ;
! d              NIW:long_name = "TYPE OF INSTRUMENT.FOR WIND MEASUREMENT" ;
! d*             MTISI1:long_name = "TIME SIGNIFICANCE" ;
! d*             NGGTP:long_name = "TIME PERIOD OR DISPLACEMENT" ;
! d*             NDNDN:long_name = "WIND DIRECTION" ;
! d*             NFNFN:long_name = "WIND SPEED" ;
!
!                MTISI2:long_name = "TIME SIGNIFICANCE" ;
!                NGGTP0:long_name = "TIME PERIOD OR DISPLACEMENT" ;
!  *             NMWGD:long_name = "MAXIMUM WIND GUST DIRECTION" ;
!  *             NFXGU:long_name = "MAXIMUM WIND SPEED(GUSTS)" ;
!                NACH2V0:long_name = "ART. CORREC. OF SENSOR H 2 ANOTHER VALUE" ;
!                MHAWAS2:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
!                MHOSEN1:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
!                MGGTP:long_name = "TIME PERIOD OR DISPLACEMENT" ;
!  *             MRRR:long_name = "TOTAL PRECIPITATION/TOTAL WATER EQUIV." ;
!                MHOSEN2:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
! d              MTISI3:long_name = "TIME SIGNIFICANCE" ;
!                MGGTP0:long_name = "TIME PERIOD OR DISPLACEMENT" ;
!  *             MSSSS:long_name = "GLOBAL RADIATION INTEGR.O.PERIOD SPECIF." ;
!                MTISI4:long_name = "TIME SIGNIFICANCE" ;
! d              NOMDP:long_name = "OPERATOR OR MANUFACTURER DEF. PARAMETER" ;
!                MDREP1:long_name = "DELAYED DESCRIPTOR REPLICATION FACTOR" ;
!                YSUPL:long_name = "008 CHARACTERS" ;
!                DATE:long_name = "Date as YYYYMMDD" ;
!                TIME:long_name = "Time as HHMMSS" ;
!
!                NORBNU:long_name = "ORBIT NUMBER" ;
!                NCTCN:long_name = "CROSS-TRACK CELL NUMBER" ;
!                NWVCQ:long_name = "ASCAT WIND VECTOR CELL QUALITY" ;
!                MSWCQ:long_name = "SEAWIND WIND VECTOR CELL QUALITY" ;
!                MPROR:long_name = "PROBABILITY OF RAIN" ;
!                NNOVA:long_name = "NUMBER OF VECTOR AMBIGUITIES" ;
!                NISWV:long_name = "INDEX OF SELECTED WIND VECTOR" ;
!                NTNSM:long_name = "TOTAL NUMBER OF SIGMA-0 MEASUREMENTS" ;
!                NFF:long_name = "WIND SPEED AT 10 M" ;
!                NDD:long_name = "WIND DIRECTION AT 10 M" ;
!                MLCFS:long_name = "LIKELIHOOD COMPUTED FOR SOLUTION" ;

  !---------------------------------------
  ! NetCDF variable IDs for body   section
  !---------------------------------------
  ! variable ID's in NetCDF file for
  !    expansion  in NetCDF file BUFR- data section4
  !
! integer    :: varid_mpppp    ! PRESSURE REDUCED TO MSL [Pa]
  !---------
  ! loop 000
  !---------
  !---------
  ! loop 001
  !---------
! integer    :: varid_mdrep    ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !--------------
  ! Scatterometer
  !--------------
  integer    :: varid_norbnu   ! Orbit Number
  integer    :: varid_nctcn    ! Cross-Track Cell Number
  integer    :: varid_nwvcq    ! Ascat Wind Vector Cell Quality
  integer    :: varid_mswcq    ! Seawind Wind Vector Cell Quality
  integer    :: varid_mpror    ! Probability of Rain
  integer    :: varid_nnova    ! Number of Vector Ambiguities
  integer    :: varid_niswv    ! Index of selected Wind Vector
  integer    :: varid_ntnsm    ! TOTAL NUMBER OF SIGMA-0 MEASUREMENTSor
  integer    :: varid_nff      ! Scatterometer Wind Speed [M/S]
  integer    :: varid_ndd      ! Scatterometer Wind Direction [Degree]
  integer    :: varid_mlcfs    ! Likelihood of computed solution
  !----------
  ! Altimeter
  !----------
  integer    :: varid_nrsst    ! Remotely Sensed Surface Type
  integer    :: varid_nrasst   ! Radiometer Sensed Surface Type
  integer    :: varid_nadqf    ! Altimeter Data Quality Flag
  integer    :: varid_nacqf    ! Altimeter Correction Quality Flag
  integer    :: varid_narf     ! Altimeter Rain Flag
  integer    :: varid_nipi     ! Ice Presence Indicator
  integer    :: varid_nfffa    ! Wind Speed from Altimeter [m/s]
  integer    :: varid_nodle    ! Ocean depth/land elevation [m]
  integer    :: varid_nrms20   ! RMS of 20 Hz Ku band ocean range [m]
  integer    :: varid_nkswh    ! Ku Band significant wave height  [m]
  integer    :: varid_nr20sw   ! RMS 20 HZ Ku band significant wave height [m]
  integer    :: varid_nrmsor   ! RMS of specific band ocean range
  integer    :: varid_nsbswh   ! Specific band significant wave height
  integer    :: varid_nrmswh   ! RMS specific band significant wave height [m]

  !-----------------------------------------
  ! quantities derived from the BUFR message
  !-----------------------------------------
  integer            :: icount ! Zaehler
  type(t_scatt)      :: s      ! SCATT observation read
  !------------
  ! temporaries
  !------------
  integer  :: dimids (dimids_max)
!------------------------------------------------------------------------------

  lpr_scatt = .false.; if (netcdf_verb > 1) lpr_scatt = .true.
  lpr_extd  = .true.
  npr_extd  =   2
  icount    =   0
  if (lpr_extd) npr_extd  =  100

  !------------------------------
  ! get default number of reports
  !------------------------------
  nreport     =  recl - rec1 + 1

  lkeep       = .false.
  nkeep       = 0

  report_subt = s2ikz(1)
  l_scat      = any (report_subt == kz_scat)
  l_altim     = .false.
  !---------------------------
  ! get variables ids and name
  !---------------------------
  status = nf90_Inquire (ncid, nVariables=ncvars)
  if ( lpr_scatt .and. lpr_extd ) then
    do j = 1 , ncvars
    status = nf90_Inquire_Variable(ncid, j,name=yname_v )
    write (6,'(a,i3,a,a)') 'nf90_Inquire_Variable(', j,          &
                           ')  Variable name(o): ', trim (yname_v)
    end do
  endif
  !========================
  ! get dimension of fields
  !========================
! set some defaults
  len_layer  = 0
  len_report = nreport
  nc1 = len_layer
  nc2 = len_report
  nc3 = 1
  nc4 = 1
  nc5 = 1
  nc6 = 1
  len_l_fxg   = 1
  len_report1 = nreport

  if ( l_scat) then

    !--------------------------
    ! default for scatterometer
    !--------------------------
    nc5       = 1
    nc6       = 1
    len_l_fxg = 1
    ! get dimension of gust loop
    status = nf90_inq_varid (ncid, 'NFF' ,  varid_nff)
    if (status /= nf90_noerr) then
      nc3         = len_l_fxg
      nc4         = len_report1
    else
      status = nf90_Inquire_Variable(ncid, varid_NFF, ndims=numDims, dimids=dimids, natts=numAtts)

      nc3 = 0
      nc4 = 0

      if(numDims >= 2) then
        dimid_l_fxg  = dimids(1)
        dimid_report = dimids(2)
        status = nf90_Inquire_Dimension(ncid, dimid_l_fxg,  len=len_l_fxg )
!       status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report1)
        nc3       = len_l_fxg
        nc4       = len_report1
      endif
    endif
    if (netcdf_verb > 0) WRITE(6,*) 'NC3/NC4', nc3, nc4
  endif

  !===================================
  ! allocate fields for reading netcdf
  !===================================

  allocate ( ifd_1d  (len_report) )
  allocate ( rfd_1d  (len_report) )
  allocate ( ifd_2d1 (nc1,nc2   ) )
  allocate ( rfd_2d1 (nc1,nc2   ) )
  allocate ( ifd_2d2 (nc3,nc4   ) )
  allocate ( rfd_2d2 (nc3,nc4   ) )

  allocate ( ifd_1dm (  1       ) )
  allocate ( rfd_1dm (  1       ) )
  allocate ( ifd_2dm (  1,  1   ) )
  allocate ( rfd_2dm (  1,  1   ) )

  if ( (.not. l_scat .and. lpr_scatt) ) then
    write (6,'()')
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' BEGINN  MO_SCATT        !!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    write (6,'(a,i4,/,6x,a,i4,/,6x,a,i4,/,6x,a,10i4)')   'pe=',dace% pe,            &
                            ' varid_MRRR  number of dimensions: ',numDims,          &
                            ' varid_MRRR  number of attributes: ',numAtts,          &
                            ' varid_MRRR  ids    of dimensions: ',dimids(1:numDims)

    write (6,'(a,i3,a,i8,/ , &
             & 1x  ,a,i8,/ , &
             & 1x  ,a,i10,a,i10) ') &
                        'pe=',dace% pe,  ' Reports in      BUFR File:',nreport,      &
                                         ' Periods fxgu in BUFR File:',len_l_fxg  ,  &
                                         ' nc5=',nc5,' nc6=',nc6
  endif
  !----------------------------------------
  ! define number of reports in source-file
  !----------------------------------------
  i_source = len_report

  !----------------
  ! allocate fields
  !----------------
! allocate ( mpppp   (len_report))     ! PRESSURE REDUCED TO MSL [Pa]
!
! scatterometer
!
  allocate (norbnu   (len_report))              ! ORBIT NUMBER
  allocate (nctcn    (len_report))              ! CROSS-TRACK CELL NUMBER
  allocate (nwvcq    (len_report))              ! ASCAT WIND VECTOR CELL QUALITY
  allocate (mswcq    (len_report))              ! SEAWIND WIND VECTOR CELL QUALITY
  allocate (mpror    (len_report))              ! PROBABILITY OF RAIN
  allocate (nnova    (len_report))              ! NUMBER OF VECTOR AMBIGUITIES
  allocate (niswv    (len_report))              ! INDEX OF SELECTED WIND VECTOR
  allocate (ntnsm    (len_report))              ! TOTAL NUMBER OF SIGMA-0 MEASUREMENTS
! allocate (mdreps   (len_report))              ! DELAYED DESCRIPTOR REPLICATION FACTOR
  allocate (nff      (len_l_fxg,len_report))    ! SCATTEROMETER WIND VELOCITY
  allocate (ndd      (len_l_fxg,len_report))    ! SCATTEROMETER WIND DIRECTION
  allocate (mlcfs    (len_l_fxg,len_report))    ! LIKELIHOOD COMPUTED FOR SOLUTION
  allocate (mle_flag (len_report))              ! QUALITY FLAG DER SCATTEROMETER DATEN
  !----------
  ! Altimeter
  !----------
  allocate (nrsst    (len_report))              ! REMOTELY SENSED SURFACE TYPE
  allocate (nrasst   (len_report))              ! RADIOMETER SENSED SURFACE TYPE
  allocate (nadqf    (len_report))              ! ALTIMETER DATA QUALITY FLAG
  allocate (nacqf    (len_report))              ! ALTIMETER CORRECTION QUALITY FLAG
  allocate (narf     (len_report))              ! ALTIMETER RAIN FLAG
  allocate (nipi     (len_report))              ! ICE PRESENCE INDICATOR
  allocate (nfffa    (len_report))              ! WIND SPEED FROM ALTIMETER
  allocate (nodle    (len_report))              ! OCEAN DEPTH/LAND ELEVATION
  allocate (nrms20   (len_report))              ! RMS OF 20 HZ KU BAND OCEAN RANGE
  allocate (nkswh    (len_report))              ! KU BAND SIGNIFICANT WAVE HEIGHT
  allocate (nr20sw   (len_report))              ! RMS 20 HZ KU BAND SIGNIFICANT WAVE HEIGHT
  allocate (nrmsor   (len_report))              ! RMS OF SPECIFIC BAND OCEAN RANGE
  allocate (nsbswh   (len_report))              ! SPECIFIC BAND SIGNIFICANT WAVE HEIGHT
  allocate (nrmswh   (len_report))              ! RMS SPECIFIC BAND SIGNIFICANT WAVE HEIGHT

  !=====================
  ! get data from netcdf
  !=====================
  i_dft =    -1
  r_dft = rvind
  !-------------------------------
  ! pressure reduced to msl  [Pa]
  !-------------------------------
! mpppp   = r_dft
! l_mpppp = .false.
! l_mdrep = .false.

  !--------------
  ! scatterometer
  !--------------
  l_nrsst  = .false.; l_nrasst = .false.; l_nadqf  = .false.; l_nacqf = .false.
  l_narf   = .false.; l_nipi   = .false.; l_nfffa  = .false.; l_nodle = .false.
  l_nrms20 = .false.; l_nkswh  = .false.; l_nr20sw = .false.  ! Jason
  l_nrmsor = .false.; l_nsbswh = .false.; l_nrmswh = .false.  ! SARAL
  if (l_scat) then
    !--------------
    ! Wind velocity
    !--------------
    status = nf90_inq_varid (ncid, 'NFF' ,  varid_nff)
    if (status /= nf90_noerr) then
       status = nf90_inq_varid (ncid, 'NFFFA', varid_nff)
       if (status == nf90_noerr) l_altim = .true.
    end if
    status = nf90_Inquire_Variable(ncid, varid_nff, ndims=numDims, dimids=dimids, natts=numAtts)
    if (netcdf_verb > 0) then
      WRITE(6,*) 'NORBNU:', numDims
      WRITE(6,*) 'NORBNU:', dimids(1),dimids(2)
      WRITE(6,*) 'NORBNU:', numAtts,'  ',len_report
    end if
    if(numDims >= 2) then
      dimid_l_fxg  = dimids(1)
      dimid_report = dimids(2)
      status = nf90_Inquire_Dimension(ncid, dimid_l_fxg,  len=len_l_fxg )
      status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report1)
    endif
    if (netcdf_verb > 0) WRITE(6,*) 'NFF:   ', len_l_fxg, len_report1

    call get_dat(ncid, norbnu , rfd_1dm, ifd_2dm, rfd_2dm, 'NORBNU',  varid_norbnu, ty_i, ndim1, &
                       ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_norbnu ,rec1)
    if (netcdf_verb > 0) &
      write(6,*) 'MAX/MIN FIELD NORBNU:', maxval(norbnu), minval(norbnu)

    !----------------
    ! Altimeter data:
    !----------------
    if (l_altim) then
     !-----------------------------
     ! Remotely Sensed Surface Type
     !-----------------------------
     call get_dat(ncid, nrsst  , rfd_1dm, ifd_2dm, rfd_2dm, 'NRSST',  varid_nrsst, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nrsst,rec1)
     !-------------------------------
     ! Radiometer Sensed Surface Type
     !-------------------------------
     call get_dat(ncid, nrasst , rfd_1dm, ifd_2dm, rfd_2dm, 'NRASST', varid_nrasst,ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nrasst,rec1)
     !----------------------------
     ! Altimeter Data Quality Flag
     !----------------------------
     call get_dat(ncid, nadqf  , rfd_1dm, ifd_2dm, rfd_2dm, 'NADQF',  varid_nadqf, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nadqf,rec1)
     if (.not. l_nadqf) then
      ! Try "Band specific altimeter data quality flag" instead:
      call get_dat(ncid,nadqf  , rfd_1dm, ifd_2dm, rfd_2dm, 'MBSADF', varid_nadqf, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nadqf,rec1)
     end if
     !----------------------------------
     ! Altimeter Correction Quality Flag
     !----------------------------------
     call get_dat(ncid, nacqf  , rfd_1dm, ifd_2dm, rfd_2dm, 'NACQF',  varid_nacqf, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nacqf,rec1)
     if (.not. l_nacqf) then
      ! Try "Band specific altimeter correction quality flag" instead:
      call get_dat(ncid,nacqf  , rfd_1dm, ifd_2dm, rfd_2dm, 'MBSACF', varid_nacqf, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nacqf,rec1)
     end if
     !--------------------
     ! Altimeter Rain Flag
     !--------------------
     call get_dat(ncid, narf   , rfd_1dm, ifd_2dm, rfd_2dm, 'NARF' ,  varid_narf , ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_narf ,rec1)
     !-----------------------
     ! Ice Presence Indicator
     !-----------------------
     call get_dat(ncid, nipi   , rfd_1dm, ifd_2dm, rfd_2dm, 'NIPI' ,  varid_nipi , ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nipi ,rec1)
     !--------------
     ! Wind velocity
     !--------------
     call get_dat(ncid, ifd_1dm, nfffa  , ifd_2dm, rfd_2dm, 'NFFFA',  varid_nfffa, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nfffa,rec1)
     !---------------------------
     ! Ocean depth/land elevation
     !---------------------------
     call get_dat(ncid, ifd_1dm, nodle  , ifd_2dm, rfd_2dm, 'NODLE',  varid_nodle, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nodle,rec1)
     !---------------------------------
     ! RMS of 20 Hz Ku band ocean range
     !---------------------------------
     call get_dat(ncid, ifd_1dm, nrms20 , ifd_2dm, rfd_2dm, 'NRMS20', varid_nrms20, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nrms20, rec1)
     if (.not. l_nrms20) &
     call get_dat(ncid, ifd_1dm, nrmsor , ifd_2dm, rfd_2dm, 'NRMSOR', varid_nrmsor, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nrmsor, rec1)
     !--------------------------------
     ! Ku band significant wave height
     !--------------------------------
     call get_dat(ncid, ifd_1dm, nkswh  , ifd_2dm, rfd_2dm, 'NKSWH',  varid_nkswh, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nkswh, rec1)
     if (.not. l_nkswh) &
     call get_dat(ncid, ifd_1dm, nsbswh , ifd_2dm, rfd_2dm, 'NSBSWH', varid_nsbswh,ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nsbswh, rec1)
     !------------------------------------------
     ! RMS 20 Hz Ku band significant wave height
     !------------------------------------------
     call get_dat(ncid, ifd_1dm, nr20sw , ifd_2dm, rfd_2dm, 'NR20SW', varid_nr20sw, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nr20sw, rec1)
     if (.not. l_nr20sw) &
     call get_dat(ncid, ifd_1dm, nrmswh , ifd_2dm, rfd_2dm, 'NRMSWH', varid_NRMSWH, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_NRMSWH, rec1)

    !--------------------
    ! Scatterometer data:
    !--------------------
    else
     !---------------------
     ! Index of Cell Number
     !---------------------
     call get_dat(ncid, nctcn  , rfd_1dm, ifd_2dm, rfd_2dm, 'NCTCN',  varid_nctcn, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nctcn,rec1)
     if (netcdf_verb > 0) &
       write(6,*) 'MAX/MIN FIELD NCTCN:', maxval(nctcn), minval(nctcn)
     !-----------------------------
     ! Number of vector Ambiguities
     !-----------------------------
     call get_dat(ncid, nnova  , rfd_1dm, ifd_2dm, rfd_2dm, 'NNOVA',  varid_nnova, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nnova,rec1)
     if (netcdf_verb > 0) &
       write(6,*) 'MAX/MIN FIELD NNOVA:', maxval(nnova), minval(nnova)
     !------------------------------
     ! Index of selected Wind Vector
     !------------------------------
     call get_dat(ncid, niswv  , rfd_1dm, ifd_2dm, rfd_2dm, 'NISWV',  varid_niswv, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_niswv,rec1)
     if (netcdf_verb > 0) &
       write(6,*) 'MAX/MIN FIELD INDEX:', maxval(niswv), minval(niswv)
     !-------------------------------------
     ! Total number of sigma-0 measurements
     !-------------------------------------
     call get_dat(ncid, ntnsm , rfd_1dm, ifd_2dm, rfd_2dm, 'NTNSM',  varid_ntnsm, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_ntnsm,rec1)
     if (netcdf_verb > 0) &
       write(6,*) 'MAX/MIN FIELD SIGMA:', maxval(ntnsm), minval(ntnsm)
     !--------------------
     ! probability of rain
     !--------------------
     call get_dat(ncid, ifd_1dm, mpror  , ifd_2dm, rfd_2dm, 'MPROR',  varid_mpror, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mpror,rec1)
     if (netcdf_verb > 0) &
       write(6,*) 'MAX/MIN FIELD RAIN:', maxval(mpror), minval(mpror)

     call get_dat(ncid, nwvcq  , rfd_1dm, ifd_2dm, rfd_2dm, 'NWVCQ',  varid_nwvcq, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nwvcq,rec1)
     call get_dat(ncid, mswcq  , rfd_1dm, ifd_2dm, rfd_2dm, 'MSWCQ',  varid_mswcq, ty_i, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mswcq,rec1)
     !--------------
     ! Wind velocity
     !--------------
     call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, nff    , 'NFF'  ,  varid_nff  , ty_f, ndim2, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nff  ,rec1)
     !---------------
     ! Wind direction
     !---------------
     call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, ndd    , 'NDD' ,  varid_ndd   , ty_f, ndim2, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_ndd   ,rec1)
     !---------------------------------------
     ! likelihood computed for solution MLCFS
     !---------------------------------------
     call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, mlcfs  , 'MLCFS' , varid_mlcfs, ty_f, ndim2, &
                        ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mlcfs ,rec1)
    end if

    if (netcdf_verb > 0) then
     if (l_altim) then
      write(6,*) 'MAX/MIN FIELD NFFFA:', maxval (nfffa,mask=nfffa/=rvind), &
                                         minval (nfffa,mask=nfffa/=rvind)
     else
      write(6,*) 'MAX/MIN FIELD NFF1:', maxval (nff,mask=nff/=rvind), &
                                        minval (nff,mask=nff/=rvind)
      PRINT*, 'MAX/MIN NWVCQ = ', maxval(nwvcq),minval(nwvcq)   ! ASCAT only
     end if
    end if

    select case (report_subt)
    case (kz_scat(2),kz_scat(3))
      call mleqc_ascat(nwvcq,len_report,mle_flag,ierm)
      if (verbose > 1) PRINT*, 'ANZAHL DISMISS ASCAT P. = ', sum(mle_flag)
    case (kz_scat(1))
      call mleqc_qscat(nctcn,ntnsm,mswcq,niswv,nff,ndd,mlcfs,len_report,mle_flag,ierm)
      if (verbose > 1) PRINT*, 'ANZAHL DISMISS QSCAT P. = ', sum(mle_flag)
      call qscat_bcor(len_report,nnova,mle_flag,mpror,nff,ierm)
    case (kz_scat(4),kz_scat(7))
      call mleqc_oscat(nctcn,ntnsm,mswcq,niswv,nff,ndd,mlcfs,len_report,mle_flag,ierm)
      if (verbose > 1) PRINT*, 'ANZAHL DISMISS OSCAT/HSCAT P. = ', sum(mle_flag)
    case (kz_scat(5))
      call rawqc_altim (nrsst, nrasst, nadqf, nacqf, narf, nipi,     &
                        nfffa, nrms20, nkswh, nr20sw, nodle, mle_flag)
      if (verbose > 1) PRINT*, 'Total dismissed JASON   = ', sum (mle_flag)
    case (kz_scat(6))
      if (.not. l_narf)   narf   = 0        ! Rain flag not present in SARAL data???
      if (.not. l_nrms20) nrms20 = 0        ! Missing RMS of 20 Hz Ku band ocean range
      if (.not. l_nkswh)  nkswh  = 0        ! Missing Ku band SWH
      if (.not. l_nr20sw) nr20sw = 0        ! Missing RMS 20 Hz Ku band SWH
      if (      l_nrmsor) nrms20 = nrmsor   ! Use SARAL specific RMS of OR
      if (      l_nsbswh) nkswh  = nsbswh   ! Use SARAL specific SWH
      if (      l_nrmswh) nr20sw = nrmswh   ! Use SARAL specific RMS of SWH
      call rawqc_altim (nrsst, nrasst, nadqf, nacqf, narf, nipi,     &
                        nfffa, nrms20, nkswh, nr20sw, nodle, mle_flag)
      if (verbose > 1) PRINT*, 'Total dismissed SARAL   = ', sum (mle_flag)
    case (kz_scat(8))
      if (.not. l_nacqf)  nacqf  = 0        ! Missing altimeter correction quality
      if (.not. l_narf)   narf   = 0        ! Rain flag not present in SENTINEL data???
      if (.not. l_nrms20) nrms20 = 0        ! Missing RMS of 20 Hz Ku band ocean range
      if (.not. l_nkswh)  nkswh  = 0        ! Missing Ku band SWH
      if (.not. l_nr20sw) nr20sw = 0        ! Missing RMS 20 Hz Ku band SWH
      if (      l_nrmsor) nrms20 = nrmsor   ! Use SENTINEL specific RMS of OR
      if (      l_nsbswh) nkswh  = nsbswh   ! Use SENTINEL specific SWH
      if (      l_nrmswh) nr20sw = nrmswh   ! Use SENTINEL specific RMS of SWH
      call rawqc_altim (nrsst, nrasst, nadqf, nacqf, narf, nipi,     &
                        nfffa, nrms20, nkswh, nr20sw, nodle, mle_flag)
      if (verbose > 1) PRINT*, 'Total dismissed SENTINEL = ', sum (mle_flag)
    end select

    if (netcdf_verb > 0) then
     if (l_altim) then
      write(6,*) 'MAX/MIN FIELD NFFFA2:', maxval (nfffa,mask=(mle_flag==0)), &
                                          minval (nfffa,mask=(mle_flag==0))
     else
      write(6,*) 'MAX/MIN FIELD NFF2:', maxval (nff,mask=nff/=rvind), &
                                        minval (nff,mask=nff/=rvind)
     end if
    end if
  endif

  !----------------------------------
  ! list of defined variables in file
  !----------------------------------
  if (netcdf_verb > 0) then
! meteorological part
    write (6,'(a,i3,1x,a,/,2x, 8l9)')  'pe=',dace% pe , &
     'l_nff  l_nfffa',&
      l_nff, l_nfffa

! technological part
    write (6,'(a,i3,1x,a,/,3x, 9l9)') 'pe=',dace% pe , &
     'l_nrsst  l_nrasst l_nadqf  l_nacqf  l_narf   l_nipi   l_nkswh',&
      l_nrsst ,l_nrasst,l_nadqf ,l_nacqf, l_narf , l_nipi,  l_nkswh
    write (6,'(a,i3,1x,a,/,3x, 9l9)') 'pe=',dace% pe , &
     'l_nrms20 l_nr20sw l_nsbswh l_nrmsor l_nrmswh l_nodle',&
      l_nrms20,l_nr20sw,l_nsbswh,l_nrmsor,l_nrmswh,l_nodle
  end if

  !-------------------------------
  ! preset total number of reports
  !-------------------------------
  entry   = sum (source(1:ifile-1)% entries) + rec1 -1

  !------------------
  ! loop over reports
  !------------------
! netcdf  : is number of report
! do is = 1,min( 2,len_report)
! do is = 1,min(10,len_report)
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
  if (bufr_type   <0) bufr_type    = obsid% bufrtype
  if (bufr_subtype<1 .and. obsid% centre == WMO0_ECMWF) &
                      bufr_subtype = obsid% subtype
  obstype                          = obsid% obstype
  if (obstype < 0) cycle
  report_subti      = idb_dbk (report_subt, obstype)

  head% obstype     = obstype
  head% dbkz        = report_subt
! head% modtype     = rept_char(obstype)% mod
  head% modtype     = SCATT
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

  if ( lpr_scatt .and. is < npr_extd ) then
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

  !--------------------------------------------
  ! perform simple generic check on report type
  !--------------------------------------------
  call check_report_0 (use, head, 1)
  if (use% state <= STAT_DISMISS) then
    if (lpr_scatt .and. is < npr_extd) write (6,*) '  check_report_0: state =', use% state
    lkeep = .false.
    cycle
  endif
  !------------------
  ! create new report
  !------------------
  spt0             = empty
  spt0% use        = use
  spt0% hd         = head
  spti             = spt0
! spti% hd% subset = is

  lkeep = .false.
  call construct_scatt (s)

  !------------------------------
  ! process the following entries
  !------------------------------
  spti% corme        = max ( s1updat(is), 0)
  spti% col% c% dlat = mlah(is)
  spti% col% c% dlon = mloh(is)
  spti% z            = -999.
  spti% actual_time  = obs_time(is)
  spti% statid       = ystidn(is)
  spti% ident        = istidn(is)

  if ( lpr_scatt  .and. is < npr_extd ) then
    write (6,'()')
    write (6,'(   a21, i6  ,a, /, &
              &   a21, i6  ,   /, &
              & 2(a21,f8.3 ,   /),&
              &   a21,2i8  ,   / ,&
              &   a21, a   ,   / ,&
              &   a21, i5      / )' )                    &
          'pe= ',dace% pe,'  spti',                      &
          'spti% corme        = ', spti% corme        ,  &
          'spti% col% c% dlat = ', spti% col% c% dlat ,  &
          'spti% col% c% dlon = ', spti% col% c% dlon ,  &
          'spti% actual_time  = ', spti% actual_time  ,  &
          'spti% statid       = ', spti% statid       ,  &
          'spti% ident        = ', spti% ident
  endif

  !-----------------------------------------------------------
  ! Determine the Scatterometer Wind Ambiguity which opposes
  ! the selected Wind by 180 Degrees (or closest to 180)
  !----------------------------------------------------------
  if (l_scat .and. .not. l_altim) then
    if (nnova(is) >= 2) then
       zlargest = -999.
!      index180 = -999

!+++++++++++++++++++++++++++++
! fix for niswv = -1 CHECK !!!
!+++++++++++++++++++++++++++++
       if ( niswv(is) /= -1 ) then
         i = niswv(is)
       else
         i = 1
       endif
!      zdir1    = ndd(niswv(is),is) * d2r
       zdir1    = ndd(        i,is) * d2r

       do jj = 1 , nnova(is)
       zdir2 = ndd(jj,is) * d2r
       zdel  = pi - abs(pi - abs(zdir1 - zdir2))

       if (zdel > zlargest) then
         zlargest = zdel
!        index180 = jj
       endif
!      write(6,*) 'HIER:',jj, niswv(is),ndd(jj,is),zdir1,zdir2,zlargest,index180
       enddo
    else
!      index180 = 1
    endif
  endif
  !--------------------------
  ! Scatterometer cell number
  !---------------------------
  if (l_scat .and. .not. l_altim) then
    spti% phase  = nctcn(is)
  endif

  !-----------
  ! wind, gust
  !-----------
  if (l_scat) then
   if (l_altim) then
    !-----------------
    ! wind speed [m/s]
    !-----------------
    call set_datum (s% ff ,nfffa(is) ,spti% corme)
   else
    icount = icount + 1
    if ( niswv(is) /= -1 ) then
      i = niswv(is)
    else
      i = 1
    endif
!   write(6,*) icount, is, i, niswv(is), ndd(i,is), nff(i,is), report_subt
    !------------------------------------------------------
    ! wind direction [degree_true]                ! ->[rad]
    !------------------------------------------------------
    call set_datum (s% dd ,ndd(i,is) ,spti% corme) ! * d2r
    !-------------------------------------------------------------
    ! second wind direction [degree_true]                ! ->[rad]
    !-------------------------------------------------------------
!   call set_datum (s% dd2 ,ndd(index180,is) ,spti% corme) ! * d2r
    !-----------------
    ! wind speed [m/s]
    !-----------------
    call set_datum (s% ff ,nff(i,is) ,spti% corme)
    !------------------------
    ! second wind speed [m/s]
    !------------------------
!   call set_datum (s% ff2 ,nff(index180,is) ,spti% corme)
!   print*,'BIN NACH WINDSPEED '
   end if
  endif

  !-----------------------------------------------
  ! if station name is not set, use station number
  !-----------------------------------------------
  if (spti% statid == '' .and. spti% ident > 0) &
    write(spti% statid,'(i5.5)') spti% ident

  !--------------------------------------------------
  ! Scatterometer: 'QSCAT', 'ASCAT', 'OSCAT', 'RSCAT'
  !--------------------------------------------------
  spti% hd% satid = spti% ident ! store WMO satellite id
  s% z = s% z_msl               ! Preset used height from reported height
  select case (spti% hd% dbkz)
  case (1697)
    spti% statid = 'QSCAT'
    s% z% o  = 0._sp
    s% z% qc = 0_i2
  case (1698,1699)
    spti% statid = 'ASCAT'
    s% z% o  = 0._sp
    s% z% qc = 0_i2
  case (1700)
    select case (spti% hd% satid)
    case (421)
       spti% statid = 'OSCAT'
    case (422)
       spti% statid = 'SCATSAT1'
    case (502,990)              ! HY-2A, old dummy ID
       spti% statid = 'HSCAT'
    case (801,991)              ! ISS, old dummy ID
       spti% statid = 'RSCAT'
    end select
    s% z% o  = 0._sp
    s% z% qc = 0_i2
  case (1701,1702,1780)
    select case (spti% hd% satid)
    case ( 61)
       spti% statid = 'SENTI_3A'
    case ( 65)
       spti% statid = 'SENTI_3B'
    case (261)
       spti% statid = 'JASON-2'
    case (262)
       spti% statid = 'JASON-3'
    case (441)
       spti% statid = 'SARAL'
    end select
  case (1770)
    select case (spti% hd% satid)
    case (503)              ! HY-2-B/C/D
       spti% statid = 'HSCAT-B'
    case (504)              ! HY-2-B/C/D
       spti% statid = 'HSCAT-C'
    case (505)              ! HY-2-B/C/D
       spti% statid = 'HSCAT-D'
    end select
    s% z% o  = 0._sp
    s% z% qc = 0_i2
  end select
  !----------------
  ! standard checks
  !----------------
  lkeep = .true.
! call init_time (spti% actual_time, yyyy, mo, dd, hh, mi)
  call check_report_1 (spti)

  if (l_scat) then
    if ( mle_flag(is) == 1 ) then
       call decr_rpt_use (spti ,CHK_DATASET)
    endif
  endif


  if (spti% use% state <= STAT_DISMISS) lkeep = .false.

  if (lkeep) then
    call check_store_scatt (s, spti, obs, lkeep)
    if (lkeep) then
       nkeep = nkeep + 1
    end if
  endif

  if (prt_data) then
    print *, 'pe=',dace% pe,' mo_scatt: station-id lat lon typ subtyp lkeep'
    write(*,'(1x,a,2x,2f15.5,2i8,l6)') spti% statid, &
         spti% col% c% dlat, spti% col% c% dlon,     &
         spti% hd% buf_type, spti% hd% buf_subtype, lkeep
    print *, 'pe=',dace% pe,' mo_scatt: station-time  actual time  db-time'
    call print (spti% hd% time)
    call print (spti%     actual_time)
    call print (spti% hd% db_time)
!   print *, spti% hd% dbkz

    call print_scatt (spti, s)
  endif
  !-------------------------
  ! end of loop over reports
  !-------------------------
  end do

  !-----------
  ! deallocate
  !-----------
  if (allocated ( ifd_1d  )) deallocate ( ifd_1d  )
  if (allocated ( rfd_1d  )) deallocate ( rfd_1d  )
  if (allocated ( ifd_2d1 )) deallocate ( ifd_2d1 )
  if (allocated ( rfd_2d1 )) deallocate ( rfd_2d1 )
  if (allocated ( ifd_2d2 )) deallocate ( ifd_2d2 )
  if (allocated ( rfd_2d2 )) deallocate ( rfd_2d2 )

  if (allocated ( ifd_1dm )) deallocate ( ifd_1dm )
  if (allocated ( rfd_1dm )) deallocate ( rfd_1dm )
  if (allocated ( ifd_2dm )) deallocate ( ifd_2dm )
  if (allocated ( rfd_2dm )) deallocate ( rfd_2dm )

! if (allocated ( mpppp   )) deallocate ( mpppp   )
!
! scatterometer
!
  if (allocated (norbnu   )) deallocate ( norbnu  )
  if (allocated (nctcn    )) deallocate ( nctcn   )
  if (allocated (nwvcq    )) deallocate ( nwvcq   )
  if (allocated (mswcq    )) deallocate ( mswcq   )
  if (allocated (nnova    )) deallocate ( nnova   )
  if (allocated (niswv    )) deallocate ( niswv   )
  if (allocated (ntnsm    )) deallocate ( ntnsm   )
  if (allocated (mpror    )) deallocate ( mpror   )
! if (allocated (mdreps   )) deallocate ( mdreps  )
  if (allocated (nff      )) deallocate ( nff     )
  if (allocated (ndd      )) deallocate ( ndd     )
  if (allocated (mlcfs    )) deallocate ( mlcfs   )
  if (allocated (nfffa    )) deallocate ( nfffa   )
  if (allocated (nodle    )) deallocate ( nodle   )
  if (allocated (nrms20   )) deallocate ( nrms20  )
  if (allocated (nkswh    )) deallocate ( nkswh   )
  if (allocated (nr20sw   )) deallocate ( nr20sw  )
  if (allocated (nsbswh   )) deallocate ( nsbswh  )
  if (allocated (nrmsor   )) deallocate ( nrmsor  )
  if (allocated (nrmswh   )) deallocate ( nrmswh  )

  end subroutine read_scatt_netcdf
!------------------------------------------------------------------------------
!   subroutine get_dat(ncid, ito1  , rto1  , ito2  , rto2  , charid, varid, ty_fto, ty_ffr, &
!                      ndim, ifd_1d, rfd_1d, ifd_2d, rfd_2d, i_dft, r_dft, l_var,           &
!                      start)
  subroutine get_dat(ncid, ito1  , rto1  , ito2  , rto2  , charid, varid, ty_fto, &
                     ndim, ifd_1d, rfd_1d, ifd_2d, rfd_2d, i_dft, r_dft, l_var,   &
                     start)

  integer     ,intent(in)                    :: ncid    ! NetCDF file handle

  integer     ,intent(inout), dimension(:)   :: ito1    ! 1dim integer field (to)
  real        ,intent(inout), dimension(:)   :: rto1    ! 1dim real    field (to)
  integer     ,intent(inout), dimension(:,:) :: ito2    ! 2dim integer field (to)
  real        ,intent(inout), dimension(:,:) :: rto2    ! 2dim real    field (to)
  character(*),intent(in   )                 :: charid  ! netCDF variable name
  integer     ,intent(  out)                 :: varid   ! netCDF variable id
  integer     ,intent(in   )                 :: ty_fto  ! type of field (to) ; 1 integer ; 2 real
! integer     ,intent(in   )                 :: ty_ffr  ! type of field (from)
  integer     ,intent(in   )                 :: ndim    ! dimension of field (to)
  integer     ,intent(inout), dimension(:)   :: ifd_1d  ! 1dim integer field (from)
  real        ,intent(inout), dimension(:)   :: rfd_1d  ! 1dim integer field (from)
  integer     ,intent(inout), dimension(:,:) :: ifd_2d  ! 2dim integer field (from)
  real        ,intent(inout), dimension(:,:) :: rfd_2d  ! 2dim integer field (from)
  integer     ,intent(in   )                 :: i_dft   ! default value for integer field (to)
  real        ,intent(in   )                 :: r_dft   ! default value for real    field (to)
  logical     ,intent(  out)                 :: l_var   ! logical for variable exists and any values defined
  integer     ,intent(in   )                 :: start   ! parameter for nf90_get_var


  integer              :: status         ! NetCDF status variable
  integer              :: itype          ! NetCDF type of variable
  integer              :: ndim_aux       ! NetCDF number of dimensions
  integer              :: ty_ffr         ! type of field
  character(len=120)   :: msg = ''

  l_var = .FALSE.

  status = nf90_inq_varid (ncid,trim( charid ),  varid )
  if (status == nf90_noerr) then
    status = nf90_inquire_variable (ncid, varid, xtype=itype, ndims=ndim_aux)
  end if

  if (status == nf90_noerr) then
    l_var = .true.
    if (ndim /= ndim_aux) then
      write(msg, '(A,I1,A,I1,A)') 'unexpected number of dimensions for variable '//&
           trim(charid)//': ',ndim,' expected, but ',ndim_aux,' in NetCDF file.'
      call finish('get_dat', msg)
    end if
    if (itype==nf90_int .or. itype==nf90_byte .or. itype==nf90_short) then
      ty_ffr = TY_I
    else if (itype==nf90_float .or. itype==nf90_double) then
      ty_ffr = TY_F
    end if

    !--------------
    ! 1 dimensional
    !--------------
    if ( ndim == 1 ) then
      if ( ty_ffr == TY_I ) then
        status = nf90_get_var (ncid, varid , ifd_1d ,start=(/start/))
        l_var  = any( ifd_1d /= imissing )
        if ( ty_fto == TY_I ) then
          where ( ifd_1d == imissing )
            ito1 = i_dft
          elsewhere
            ito1 = ifd_1d
          end where
        elseif ( ty_fto == TY_F ) then
          where ( ifd_1d == imissing )
            rto1 = r_dft
          elsewhere
            rto1 = ifd_1d
          end where
        end if
      else if ( ty_fto == TY_F ) then
        status = nf90_get_var (ncid, varid , rfd_1d ,start=(/start/))
        l_var  = any( rfd_1d /= rmissing )
        where ( rfd_1d == rmissing )
          rto1 = r_dft
        elsewhere
          rto1 = rfd_1d
        end where
      end if
    !--------------
    ! 2 dimensional
    !--------------
    elseif ( ndim == 2 ) then
      if ( ty_ffr == TY_I ) then
        status = nf90_get_var (ncid, varid , ifd_2d ,start=(/1,start/))
        l_var  = any( ifd_2d /= imissing )
        if ( ty_fto == TY_I ) then
          where ( ifd_2d == imissing )
            ito2 = i_dft
          elsewhere
            ito2 = ifd_2d
          end where
        elseif ( ty_fto == TY_F ) then
          where ( ifd_2d == imissing )
            rto2 = r_dft
          elsewhere
            rto2 = ifd_2d
          end where
        end if
      else if ( ty_fto == TY_F ) then
        status = nf90_get_var (ncid, varid , rfd_2d ,start=(/1,start/))
        l_var  = any( rfd_2d /= rmissing )
        where ( rfd_2d == rmissing )
          rto2 = r_dft
        elsewhere
          rto2 = rfd_2d
        end where
      end if
    end if
  else
  ! if not defined , preset field with default
    if ( ndim == 1 ) then
      if ( ty_fto == TY_I ) then
        ito1 = i_dft
      else if ( ty_fto == TY_F ) then
        rto1 = r_dft
      endif
    else if ( ndim == 2 ) then
      if ( ty_fto == TY_I ) then
        ito2 = i_dft
      else if ( ty_fto == TY_F ) then
        rto2 = r_dft
      endif
    endif
  endif
  end subroutine get_dat
!------------------------------------------------------------------------------
  subroutine mleqc_ascat (qc_field, nn, mle_flag, ierm)
  integer     ,intent(in),    dimension(:)     ::   qc_field
  integer     ,intent(in)                      ::   nn
  integer     ,intent(out),   dimension(nn)    ::   mle_flag
  integer     ,intent(out)                     ::   ierm

  integer                                      ::   is
  integer                                      ::   NICE
  integer                                      ::   NDW
  integer                                      ::   NKNMIQC
  integer                                      ::   IMONITOR
  integer                                      ::   IVARQC
  integer                                      ::   IWSP30
  integer                                      ::   IWSP3


      ierm  = 0

!     +----------------------------------------------------------------+
!     | Initialize NMLE quality control flag                           |
!     +----------------------------------------------------------------+

       NICE     = 0
       NDW      = 24
       NKNMIQC  = 0
       IMONITOR = 0
       IVARQC   = 0
       IWSP30   = 0
       IWSP3    = 0
       mle_flag = 0

       DO IS = 1 , NN

!      +-------------------------------------------------------------+
!      | Check the Monitor not used flag                             |
!      +-------------------------------------------------------------+
       if(.not.btest(qc_field(is),NDW-4-1) .and. &
          .not.btest(qc_field(is),NDW-5-1) ) then
          IMONITOR = IMONITOR + 1
       endif
!      +-------------------------------------------------------------+
!      | Check Ice flag and Land Flag in Bufr Data                   |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),NDW-8-1) .or. &
          btest(qc_field(is),NDW-9-1) ) then
!         WRITE(*,*) 'ICE/LAND Flag: ', btest(qc_field(is),NDW-8-1), &
!         btest(qc_field(is),NDW-9-1)
          NICE = NICE + 1
          mle_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if KNMI QC FLAG is set                                |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),NDW-6-1) ) then
          NKNMIQC =  NKNMIQC + 1
          mle_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if VAR QC FLAG is set                                 |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),NDW-7-1) ) then
          IVARQC =  IVARQC + 1
          mle_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if windspeed is greater than 30 m/s                   |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),NDW-11-1) ) then
          IWSP30 =  IWSP30 + 1
          mle_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if windspeed is <= 3 m/s                              |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),NDW-12-1) ) then
          IWSP3 =  IWSP3 + 1
          mle_flag(is) = 1
          cycle
       endif

!

       ENDDO
!
    if (verbose > 1) then
       PRINT*, 'Monitoring Punkte = ', IMONITOR
       PRINT*, 'Land/Eispunkte    = ', NICE
       PRINT*, 'KNMI QC PUNKTE    = ', NKNMIQC
       PRINT*, 'VAR  QC PUNKTE    = ', IVARQC
       PRINT*, 'WSP30   PUNKTE    = ', IWSP30
       PRINT*, 'WSP3    PUNKTE    = ', IWSP3
       PRINT*, 'ANZAHL DISMISS P. = ', sum(mle_flag)
    end if

  end subroutine mleqc_ascat
!------------------------------------------------------------------------------
  subroutine mleqc_qscat (cell_number, c_incell, qc_field, w_select, &
                          w_speed, w_dir, l_sol, nn, mle_flag, ierm)
  integer     ,intent(in),    dimension(:)     ::   cell_number
  integer     ,intent(in),    dimension(:)     ::   c_incell
  integer     ,intent(in),    dimension(:)     ::   qc_field
  integer     ,intent(in),    dimension(:)     ::   w_select
  real        ,intent(in),    dimension(:,:)   ::   w_speed
  real        ,intent(in),    dimension(:,:)   ::   w_dir
  real        ,intent(in),    dimension(:,:)   ::   l_sol
  integer     ,intent(in)                      ::   nn
  integer     ,intent(out),   dimension(nn)    ::   mle_flag
  integer     ,intent(out)                     ::   ierm

  integer                                      ::   NICE
  integer                                      ::   NWDIR
  integer                                      ::   JPL_RAIN

  real :: A,v,x0,x1,y0,thresh,thresh2
  real :: coef_2ord(6,3),coef_fin(73,6)
!
  real :: exp_mle ! expected MLE
  real :: mle_n   ! normalised MLE
  real :: dtwindir

  integer :: jw, jj, jk, isel_jpl, icoef, cell_reject, is !, ispeed
#ifdef __SX__
  integer :: ibb17
#endif
  real :: maxdtwindir
!
!
  parameter(A=2./100., y0=4., x0=5., x1=15., thresh2=2.)

  DATA coef_2ord/ 0.55,1.5,2.75,0.2121,-0.0074,0.00012,     &
                  0.0,0.0,0.0,-0.00248,0.00031,-4.75e-06,   &
                  0.0,0.0,0.0,3.018e-05,-4.084e-06,6.235e-08/
!
  ierm  = 0

!     +----------------------------------------------------------------+
!     | Initialize NMLE quality control flag                           |
!     +----------------------------------------------------------------+

       NICE     = 0
       JPL_RAIN = 0
       NWDIR    = 0
       mle_flag = 0
       cell_reject = 0

       DO IS = 1 , NN

!      +-------------------------------------------------------------+
!      | Check Ice flag and Land Flag in Bufr Data                                 |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),7) .or. btest(qc_field(is),8)) then
          mle_flag(is) = 1
          NICE = NICE + 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check JPL rain flag on the nadir swath                      |
!      +-------------------------------------------------------------+
!
       if((btest(qc_field(is),2) .and. (.not. btest(qc_field(is),3))) .and. &
          (cell_number(is) > 28 .and. cell_number(is) < 49 )) then
         mle_flag(is)  = 1
         JPL_RAIN = JPL_RAIN + 1
         cycle
       endif
!      +-------------------------------------------------------------+
!      | Reject outer swath data                                     |
!      +-------------------------------------------------------------+
       if ( cell_number(is) < 11 .or. cell_number(is) > 66 ) then
          mle_flag(is)  = 1
          cell_reject = cell_reject + 1
          cycle
       endif

!      +-------------------------------------------------------------+
!      | Check if wind diversity is smaller than 135 Grad
!      +-------------------------------------------------------------+
!
       if (cell_number(is) > 28 .and. cell_number(is) < 49 &
           .and. mle_flag(is) == 0 .and. c_incell(is) == 4) then

       maxdtwindir = -999.9
       DO jw = 1 , 6
          jj = 1
          jk = jw
          IF (jw .ge. 4 .and. jw .le. 5) then
             jj = 2
             jk = jk - 3
          else IF (jw .eq. 6) then
             jj = 3
             jk = jk - 5
          endif
          dtwindir = abs(w_dir(jj,is) - w_dir(jk+jj,is))
          if(dtwindir .GT. 180.) dtwindir = dtwindir - 360.
          if(ABS(dtwindir) .GT. maxdtwindir) maxdtwindir = ABS(dtwindir)
       ENDDO

       IF(maxdtwindir .LT. 135.) THEN
          NWDIR = NWDIR + 1
          mle_flag(is)  = 1
       ENDIF

       endif
!      +--------------------------------------------------------------+
!      | Compute normalized residual Rn                               |
!      +--------------------------------------------------------------+
#ifdef __SX__
       ibb17 = bufr_bits17(qc_field(is),8,3)
       if(          ibb17               == 0 .and. c_incell(is) == 4) then
#else
       if(bufr_bits17(qc_field(is),8,3) == 0 .and. c_incell(is) == 4) then
#endif

       isel_jpl = w_select(is)

       IF (isel_jpl .LT. 0) THEN
          mle_flag(is)  = 1
          cycle
       ENDIF


       if ( w_speed(isel_jpl,is) /= rvind) then
        v = w_speed(isel_jpl,is)
!       ispeed = int(v)+1
        do icoef=1,6
           coef_fin(cell_number(is),icoef) = coef_2ord(icoef,1)   + &
                                             coef_2ord(icoef,2)*cell_number(is) + &
                                             coef_2ord(icoef,3)*cell_number(is)**2
        enddo
        exp_mle = coef_fin(cell_number(is),1) *               &
                  exp(-((v-coef_fin(cell_number(is),2))       &
                        /  coef_fin(cell_number(is),3))**2/2) &
                +     (coef_fin(cell_number(is),4) &
                + v * (coef_fin(cell_number(is),5) &
                + v *  coef_fin(cell_number(is),6)))
        mle_n   = -l_sol(isel_jpl,is)/exp_mle
!
!           +----------------------------------------------------------+
!           | Apply Rn threshold                                       |
!           +----------------------------------------------------------+
!
            if(v > x1) then
              thresh = thresh2
            else
              thresh = y0-A*(v-x0)**2
            endif
!
            if(mle_n > thresh) then
              mle_flag(is)  = 1
            endif
!
         endif !v
       endif
!
       enddo ! is

    if (verbose > 1) then
       PRINT*, 'Messpunkte_qscat     =', NN
       PRINT*, 'Land/Eispunkte_qscat =', NICE
       PRINT*, 'Wind Diversity Points=', NWDIR
       PRINT*, 'JPL RAIN FLAG        =', JPL_RAIN
       PRINT*, 'OUTER SWATH          =', cell_reject
    end if

  CONTAINS

!     function bufr_bits17
!
!CDIR INLINE
      integer function bufr_bits17 (var, bit_bufr, num)
!
!     +----------------------------------------------------------------+
!     | Declaration of variables and data statements                   |
!     +----------------------------------------------------------------+
!
        implicit   none
        integer :: i, var, bit_bufr, num, bit_fort_i, bit_bufr_i
        logical :: lexceed
!
        bufr_bits17 = 0
        lexceed     = .false.
!
        do i=1,num
           bit_bufr_i = bit_bufr+(i-1)
           bit_fort_i = 16-bit_bufr_i
           if(bit_bufr_i == 17) bit_fort_i = 16
           if(bit_bufr_i >  17) lexceed    = .true.
           bufr_bits17 = bufr_bits17 + ibits (var,bit_fort_i,1) * 2**(i-1)
        enddo
        if (lexceed) call finish ('bufr_bits17','Bit number exceeding 17')
!
        return
      end function bufr_bits17
  end subroutine mleqc_qscat
!------------------------------------------------------------------------------
  subroutine mleqc_oscat (cell_number, c_incell, qc_field, w_select, &
                          w_speed, w_dir, l_sol, nn, mle_flag, ierm)
  integer     ,intent(in),    dimension(:)     ::   cell_number
  integer     ,intent(in),    dimension(:)     ::   c_incell
  integer     ,intent(in),    dimension(:)     ::   qc_field
  integer     ,intent(in),    dimension(:)     ::   w_select
  real        ,intent(in),    dimension(:,:)   ::   w_speed
  real        ,intent(in),    dimension(:,:)   ::   w_dir
  real        ,intent(in),    dimension(:,:)   ::   l_sol
  integer     ,intent(in)                      ::   nn
  integer     ,intent(out),   dimension(nn)    ::   mle_flag
  integer     ,intent(out)                     ::   ierm

    integer                                    ::   nice
    integer                                    ::   nwsp30, nwsp3
    integer                                    ::   nkmniflag
    integer                                    ::   nkmnimon
    integer                                    ::   nkmnivar

    integer :: is, isel_jpl
!
    ierm  = 0

!   +----------------------------------------------------------------+
!   | Initialize NMLE quality control flag                           |
!   +----------------------------------------------------------------+

    nice      = 0
    nkmniflag = 0
    nkmnimon  = 0
    nkmnivar  = 0
    mle_flag  = 0
    nwsp30    = 0
    nwsp3     = 0

    DO IS = 1 , NN

!      +-------------------------------------------------------------+
!      | Check Ice flag and Land Flag in Bufr Data                   |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),7) .or. &
          btest(qc_field(is),8)) then
          mle_flag(is) = 1
          nice = nice + 1
          cycle
       endif

!      +-------------------------------------------------------------+
!      | Check Rain/sea-state KNMI Monitoring flag                   |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),10)) then
         mle_flag(is) = 1
         nkmniflag = nkmniflag + 1
         cycle
       endif

!      +-------------------------------------------------------------+
!      | Check KNMI Monitoring Flag and value                        |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),11) .or. &
          btest(qc_field(is),12)) then
          mle_flag(is) = 1
          nkmnimon = nkmnimon + 1
          cycle
       endif

!      +-------------------------------------------------------------+
!      | Check Variational QC flag                                   |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),9)) then
         mle_flag(is) = 1
         nkmnivar = nkmnivar + 1
         cycle
       endif

!      +-------------------------------------------------------------+
!      | Check if windspeed is greater than 30 m/s                   |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),5)) then
          nwsp30 =  nwsp30 + 1
          mle_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if windspeed is <= 3 m/s                              |
!      +-------------------------------------------------------------+
!
       if(btest(qc_field(is),4)) then
          nwsp3 =  nwsp3 + 1
          mle_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check validity of index of selected wind                    |
!      +-------------------------------------------------------------+

       isel_jpl = w_select(is)
       IF (isel_jpl < 1) THEN
          mle_flag(is)  = 1
          cycle
       ENDIF
!      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      | Reject outer swath data +++ DISABLED; use /rules/ phase=... |
!      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      if ( cell_number(is) < 5 ) then
!         mle_flag(is)  = 1
!         cycle
!      endif

    enddo ! is

    if (verbose > 1) then
       PRINT*, 'Messpunkte Oscat     =', NN
       PRINT*, 'Land/Eispunkte Oscat =', nice
       PRINT*, 'KNMI QC FLAG         = ',nkmniflag
       PRINT*, 'KNMI MON Flag        = ',nkmnimon
       PRINT*, 'KNMI VAR Flag        = ',nkmnivar
       PRINT*, 'WSP > 30 m/sec       = ',nwsp30
       PRINT*, 'WSP <= 3 m/sec       = ',nwsp3
    end if

  end subroutine mleqc_oscat
!---------------------------------------------------------------------------------
 subroutine qscat_bcor (nn, c_incell, mle_flag, rain_index, w_speed, ierm)
  integer     ,intent(in)                      ::   nn
  integer     ,intent(in),    dimension(:)     ::   c_incell
  integer     ,intent(in),    dimension(:)     ::   mle_flag
  real        ,intent(in),    dimension(:)     ::   rain_index
  real        ,intent(inout), dimension(:,:)   ::   w_speed
  integer                                      ::   ierm

  real     :: a1,ra1
  real(wp) :: c1, b1                !+++ required for xlf V12.1
  real(wp) :: v0,v1,v2
! real     :: delv
  integer  :: ii, kk

  ierm  = 0
!Vorbelegen der Felder---------------------------------------------------------
 a1    = 0.112;                  ra1 = 1./a1
 c1    = 8.87

 PRINT*,' BIN IN BCOR NN: ', NN, rvind

 DO ii = 1 , NN
  IF (c_incell(ii) > 0 .and. mle_flag(ii) == 0 ) then
!   b1 = rain_index(ii)
!!  if(b1 < 0.) b1 = 0.
!   if(b1 == rvind) b1 = 0.
!   b1 = b1 * c1
   if (rain_index(ii) == rvind) then
     b1 = 0.
   else
     b1 = rain_index(ii) * c1   !+++ evaluated by xlf V12.1 in any case
   endif
   DO kk = 1, c_incell(ii)
     v0 = w_speed(kk,ii)
     v1 = v0 * 0.96
     v2 = v1-v1*b1/(ra1+v1)
     w_speed(kk,ii) = v2
!    delv = v2-v0
!    print*,'HIER: ', ii,kk,b1,v0,v1,v2,delv
   ENDDO
! ELSE
! delv    = 0.0
  ENDIF
 ENDDO

 end subroutine qscat_bcor
!==============================================================================
  subroutine rawqc_altim (nrsst, nrasst, nadqf, nacqf, narf, nipi,    &
                          nfffa, nrms20, nkswh, nr20sw, nodle, qc_flag)
    !----------------------------------------------------
    ! Quality control of raw and auxiliary altimeter data
    !----------------------------------------------------
    integer ,intent(in)  :: nrsst  (:)  ! Remotely Sensed Surface Type
    integer ,intent(in)  :: nrasst (:)  ! Radiometer Sensed Surface Type
    integer ,intent(in)  :: nadqf  (:)  ! Altimeter data quality flag
    integer ,intent(in)  :: nacqf  (:)  ! Altimeter correction quality flag
    integer ,intent(in)  :: narf   (:)  ! Altimeter rain flag
    integer ,intent(in)  :: nipi   (:)  ! Ice presence indicator
    real    ,intent(in)  :: nfffa  (:)  ! Wind speed from altimeter
    real    ,intent(in)  :: nrms20 (:)  ! RMS of 20 Hz Ku band ocean range
    real    ,intent(in)  :: nkswh  (:)  ! Ku band significant wave height
    real    ,intent(in)  :: nr20sw (:)  ! RMS 20 Hz Ku band sign. wave height
    real    ,intent(in)  :: nodle  (:)  ! Ocean depth/land elevation
    integer ,intent(out) :: qc_flag(:)  ! Resulting Overall QC flag

    integer :: nsurf
    integer :: nice
    integer :: nrain
    integer :: ndq
    integer :: ncq
    integer :: IWSP3
    integer :: IWSP30
    integer :: irms20
    integer :: ikswh
    integer :: ir20sw
    integer :: iodle
    integer :: is

    nsurf   = 0
    ndq     = 0
    ncq     = 0
    nice    = 0
    nrain   = 0
    IWSP3   = 0
    IWSP30  = 0
    irms20  = 0
    ikswh   = 0
    ir20sw  = 0
    iodle   = 0
    qc_flag = 0

    do is = 1, size (nfffa)
!      +-------------------------------------------------------------+
!      | Check Remotely Sensed Surface Type   (code table 0 08 029)  |
!      | 0 = open ocean or semi-enclosed sea                         |
!      +-------------------------------------------------------------+
       if (nrsst(is) /= 0) then
          nsurf       = nsurf + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check Radiometer Sensed Surface Type (code table 0 08 077)  |
!      | 3 = open ocean or semi-enclosed sea                         |
!      +-------------------------------------------------------------+
       if (nrasst(is) /= 3) then
          nsurf       = nsurf + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check Ice presence flag              (code table 0 21 169)  |
!      +-------------------------------------------------------------+
       if (nipi(is) /= 0) then
          nice        = nice + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check Rain flag                                             |
!      +-------------------------------------------------------------+
       if (narf(is) /= 0) then
          nrain       = nrain + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check Altimeter Data Quality flag                           |
!      +-------------------------------------------------------------+
       if(nadqf(is) /= 0) then
          ndq         = ndq + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check Altimeter Correction Quality flag                     |
!      +-------------------------------------------------------------+
       if(nacqf(is) /= 0) then
          ncq         = ncq + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if RMS of 20 Hz Ku band ocean range > 20 cm           |
!      +-------------------------------------------------------------+
       if(nrms20(is) > 0.2_sp) then
          irms20      = irms20 + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if Ku band significant wave height > 20 m             |
!      +-------------------------------------------------------------+
       if(nkswh(is) > 20._sp) then
          ikswh       = ikswh + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if RMS 20 Hz Ku band significant wave height > 1.0 m  |
!      +-------------------------------------------------------------+
       if(nr20sw(is) > 1.0_sp) then
          ir20sw      = ir20sw + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if ocean depth/land elevation above -20.0 m           |
!      +-------------------------------------------------------------+
       if(nodle(is) > -20._sp) then
          iodle       = iodle + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if windspeed is <= 3 m/s                              |
!      +-------------------------------------------------------------+
       if(nfffa(is) <= 3._sp) then
          IWSP3       = IWSP3 + 1
          qc_flag(is) = 1
          cycle
       endif
!      +-------------------------------------------------------------+
!      | Check if windspeed is greater than 30 m/s                   |
!      +-------------------------------------------------------------+
       if(nfffa(is) > 30._sp) then
          IWSP30      = IWSP30 + 1
          qc_flag(is) = 1
          cycle
       endif

    end do

    if (verbose > 1) then
       PRINT*, 'Sensed surface type  = ', nsurf
       PRINT*, 'Altimeter data qual. = ', ndq
       PRINT*, 'Altimeter corr.qual. = ', ncq
       PRINT*, 'Ice present          = ', nice
       PRINT*, 'Rain flagged         = ', nrain
       PRINT*, 'RMS Ku OR   > 0.2 m  = ', irms20
       PRINT*, 'Ku band SWH >  20 m  = ', ikswh
       PRINT*, 'RMS Ku SWH  > 1.0 m  = ', ir20sw
       PRINT*, 'Ocean depth > -20 m  = ', iodle
       PRINT*, 'Wind speed <= 3 m/s  = ', IWSP3
       PRINT*, 'Wind speed > 30 m/s  = ', IWSP30
       PRINT*, 'Total dismissed data = ', sum (qc_flag)
    end if
  end subroutine rawqc_altim
!==============================================================================
end module mo_scatt
