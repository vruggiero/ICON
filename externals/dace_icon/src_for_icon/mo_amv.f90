!
!+ Routines specific to the AMV (SATOB) observation operator
!
MODULE mo_amv
!
! Description:
!   Routines specific to the AMV (Atmospheric Motion Vector)
!   observation operator
!
! Current Maintainer: DWD, Alexander Cress
!    phone: +49 69 8062 2716
!    fax:   +49 69 8062 3721
!    email: alexander.cress@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  Changes for verification mode
! V1_5         2009/05/25 Andreas Rhodin
!  Read observations from feedback-file
! V1_7         2009/08/24 Andreas Rhodin
!  new namelist parameter "cnt_redu"
! V1_8         2009/12/09 Andreas Rhodin
!  add WMO0_CIMSS to generating centers (in addition to WMO0_NASA)
! V1_9         2010/04/20 Andreas Rhodin
!  preparations for GOES-13
!  TSK_SHRINK in subroutines process: pass parameter 'state' to shrink_report
!  remove obsolete namelist variable cnt_redu (use namelist THINNING instead)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  read_amv_netcdf: reduce verbosity for netcdf_verb==0
! V1_13        2011/11/01 Harald Anlauf
!  New namelist /OBSERR_AMV/
! V1_15        2011/12/06 Andreas Rhodin
!  option to specify codetypes for observation error tables, account for GOES-15
! V1_16        2011/12/09 Alexander Cress
!  account for GOES-15 6.5um water vapour channel
! V1_21        2013/01/18 Andreas Rhodin
!  extension for Meteosat 10
! V1_22        2013-02-13 Harald Anlauf
!  implement forward operator for layer wind; misc. code cleanups
! V1_28        2014/02/26 Andreas Rhodin
!  changed interface to new_int
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_43        2015-08-19 Harald Anlauf
!  Implement AMV height correction method of Kathrin Folger/LMU
!  prepare for AMV MISR data
! V1_45        2015-12-15 Alexander Cress
!  changes for Himawari-8/9
! V1_47        2016-06-06 Harald Anlauf
!  qi-index for NPP-AMV processed by NOAA; changes for AMV height correction
! V1_49        2016-10-25 Harald Anlauf
!  make sure nfnfn is always initialized
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2003-2008  original code
! Oliver Schmid   DWD  2005       new obs data type
! Harald Anlauf   DWD  2008       optimizations for SX8
! Gerhard Paul    DWD  2008       NetCDF input
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
  use mo_exception, only: finish           ! abort routine
  use mo_kind,      only: wp, sp, i8       ! working, single precision, 64Bit integer
  use mo_namelist,  only: position_nml,   &! position namelist
                          nnml,           &! namelist Fortran unit number
                          POSITIONED       ! ok    code from position_nml
  use mo_mpi_dace,  only: dace,           &! MPI group info
                          p_sum,          &! generic sum routine
                          p_bcast          ! broadcast routine
  use mo_time,      only: init_time,      &! initialise time data type
                          operator (==)    ! compare times
  use mo_dec_matrix,only: t_vector_segm,  &! vector segment data type
                          mp               ! (sparse) matrix precision
  use mo_physics,   only: d2r,            &! pi/180.
                          gacc             ! gravity acceleration
  use mo_usstd,     only: p_h_usstd        ! pressure from geopotential height
  use mo_run_params,only: input,          &! path name of input  files
                          path_file        ! concatenate: path/file.sufx
  use mo_dace_string,only:toupper          ! Convert to upper case
  use mo_fortran_units, &
                    only: get_unit_number,  &! get a free unit number
                          return_unit_number ! release the unit number
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_state, only: t_atm            ! atm. state data type
  use mo_t_col,     only: t_cols,         &! model columns data type
                          COL_UV,         &! specification of fields in column
                          COL_GEO          !
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_t_obs,     only: t_obs,          &!
                          t_spot,         &!
                          t_head,         &! observation data type
                          derive_dbkz,    &! derive DBKZ if not present
                          new_spot,       &! reserve memory
                          new_obs,        &! reserve memory
                          new_int,        &! reserve memory
                          set_xuv,        &! set unit vectors, solar zenith
                          invalid,        &! invalid observation value
                          set_vqc_insitu, &! subroutine to set VQC bounds
                          source,         &! list   of Report source files
                          pcc_conv,       &! per cent confidence usage flag
                          TSK_INIT,       &! FLAGS: initialize module
                          TSK_READ,       &!  read observations
                          TSK_SET_CHR,    &!  set observation characteristics
                          TSK_SETUP_COLS, &!  determine model columns required
                          TSK_SETUP_FUL0, &!  setup description of PSAS-space
                          TSK_SETUP_FULL, &!  setup description of PSAS-space
                          TSK_SHRINK,     &!  release unused obs. in report
                          TSK_K,          &!  evaluate    linear operator
                          TSK_Y,          &!  evaluate nonlinear operator
                          TSK_R,          &!   setup observational errors
                          tsk_name,       &! derive mnemonic for TSK_.. flags
                          CHR_ID,         &! H is the identity operator
                          CHR_LIN,        &! H is a linear operator
                          ITY_ICOL,       &! interpolation type: column
                          netcdf_verb,    &! verbosity of NetCDF decoding
                          z2p_amv,        &! derive pressure from height ?
                          varno_oq         ! varno to observed quantity
  use mo_fdbk_tables,only:VN_U,           &!          wind component code
                          VN_V,           &!                         code
                          VN_P,           &!                pressure code
                          VN_HEIGHT,      &!                  height code
                          OT_SATOB         ! obstype id
  use mo_obs_set,   only: t_obs_block      ! obs data type
  use mo_wmo_tables,only: WMO0_ECMWF,     &! generating center
                          WMO0_JMA,       &! generating center
                          WMO0_RSMC,      &! generating center
                          WMO0_NOAA,      &! generating center
                          WMO0_NASA,      &! generating center
                          WMO0_CIMSS,     &! generating center
                          WMO0_EUMET       ! generating center
  use mo_t_use,     only: t_use,          &! status variable data type
                          use_0,          &! default values of type use
                          STAT_DISMISS,   &!
                          CHK_INSDAT,     &!
                          CHK_NOTUSED,    &!
                          CHK_SURF,       &!
                          CHK_DOMAIN,     &!
                          CHK_HEIGHT,     &!
                          CHK_QI           !
  use mo_obs_tables,only: check_report_0, &! init. flags, standard checks
                          check_report_1, &! standard checks
                          idb_dbk,        &! index in table rept_stat
                          rept_char,      &! observation type characteristics
                          decr_rpt_use     ! change use-flags of report
  use mo_t_datum,   only: t_datum,        &! data type for one observed datum
                          rvind,          &! missing value indicator (real)
!#if defined(__SX__)
!                         inv_datum,      &! invalid datum
!#endif
                          set_datum,      &! set t_datum% o  (observed value)
                          QC_OK,QC_MISS,  &! quality control flag value
                          SRC_DOK          ! data source     flag value
  use mo_bufr_dwd,  only: t_bufr,         &! BUFR record data type
!                         bufr3_get_character_value, &!
#ifdef CHECKCODES
                          bufr_get_entry_texts,      &!
                          bufr_get_entry_units,      &!
#endif
                          bufr_print_sections,       &!
                          bufr_print_subset,         &!
                          inv_bufr         ! indicator for invalid value
  use mo_satid,     only: satname,        &! derive satellite name from satid
                          satid,          &! derive satellite id from name
                          satid_longname   ! derive satellite id from longname
  use mo_obs_err,   only: obs_err          ! get observation error from table
  use mo_obs_rules, only: iud,            &! undefined integer value
                          rud              ! undefined real    value
  use mo_obstypes,only:   t_obsid,        &! observation id table entry
                          obstype_dbkz,   &! derive obsids from dbkz
                          obstype_bufr     ! derive obsids from bufr type
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only:ncid,           &! NetCDF file id
                          dimids_max,     &! max number of NetCDF dimension ids
                          imissing,       &! NetCDF _FillValue for integer
                          rmissing,       &! NetCDF _FillValue for reals
                          s1date,         &! nominal (synoptic) date
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
                          ystidn           ! any type of station identifier as variable
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,       only: nf90_Inquire,          &!
                          nf90_Inquire_Dimension,&!
                          nf90_Inquire_Variable, &!
                          nf90_inq_varid,        &!
                          nf90_get_var,          &!
                          nf90_strerror,         &!
                          NF90_MAX_NAME,         &!
                          NF90_NOERR              !
  !------------------------------------
  ! Interface to interpolation routines
  !------------------------------------
  use mo_grid_intpol,only: idx_init
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_amv      ! general purpose AMV processing routine
  public :: read_nml_amv_obs ! read AMV specific namelist /AMV_OBS/
  public :: read_amv_bufr    ! read AMV observation from BUFR record
  public :: read_amv_netcdf  ! read AMV observation from netCDF file
  public :: check_amv_surf   ! check for sea/land surface
  public :: check_amv_qual   ! check AMV quality index
  public :: clean_amv_nml    ! deallocate linked list
  public :: amv_heightcorr   ! apply AMV height correction
!------------------------------------------------------------------------------
  !======================
  ! data type definitions
  !======================
  !-----------------------------------------------
  ! type t_amv: temporarily store AMV observations
  !-----------------------------------------------
  type t_amv
    type (t_datum) :: total ! quality control flag
    type (t_datum) :: p     ! pressure level             [Pa]
    type (t_datum) :: z     ! height of top of cloud     [m]
    type (t_datum) :: t     ! temperature                [K]
    type (t_datum) :: ff    ! wind speed                 [m/s]
    type (t_datum) :: dd    ! wind direction             [degree]
    type (t_datum) :: f     ! frequency
  end type t_amv
!------------------------------------------------------------------------------
  !-----------------
  ! Namelist AMV_OBS
  !-----------------
  integer ,parameter :: nsat = 8 ! max number of satellites
  integer ,parameter :: nm   = 8 ! max number of computational methods
  integer ,parameter :: nc   = 2 ! max number of centers with reduced pcc

  type t_amv_select ! parameters for AMV data selection
    integer                     :: satids(nsat)=  -1        ! satellite ids
    integer                     :: meths (nm)  =  -1        ! retrieval methods
    integer                     :: qimin (3,3) =  999       ! min QI value
    real(wp)                    :: p_low       =  70000._wp ! lower pres.bound
    real(wp)                    :: p_high      =  40000._wp ! high  pres.bound
    real(wp)                    :: lat_N       =   20._wp   ! latitude TR <> NH
    real(wp)                    :: lat_S       =  -20._wp   ! latitude TR <> SH
    logical                     :: seaonly(3)  =  .false.   ! use only over sea
    integer                     :: fwmodel     =    0       ! forward model 0-3
    real(wp)                    :: params (6)  =    0._wp   ! parameter list
    type(t_amv_select) ,pointer :: next        => NULL()    ! link(ed list)
  end type t_amv_select

  type(t_amv_select) ,pointer :: amv_select    => NULL()

  integer  :: chk_amv_pass  = 1       ! check AMV quality index
  logical  :: chk_amv_cdfin = .true.  ! .. when reading cdfin
  logical  :: chk_amv_post  = .false. ! .. later
!______________________________________________________________________________
! Remarks on the parameters of the AMV forward model/layer wind operator:
!
!   fwmodel: operator selection (0 = standard, n=1..3 "n-point formula")
!   params:  sets of pairs of values for layer offset and layer depth for
!            low, medium and high clouds (or one common set for all levels).
!
! Example:
!   fwmodel = 3
!   params  = -25 100   0 100   20 100
!
!   Forward model: Simpson's rule (3-point) over a symmetric interval
!   (width: 100 hPa) centered about:
!
!     p_obs - 25 hPa           (p_obs > p_low),
!     p_obs +  0 hPa   (p_low > p_obs > p_high),
!     p_obs + 20 hPa           (p_obs < p_high)
!______________________________________________________________________________
!------------------------------------------------------------------------------
  !-----------------------------------------------
  ! Satellite specific observation errors for AMVs
  !-----------------------------------------------
  integer  ,parameter :: nlev      = 15
  real(wp) ,parameter :: pt (nlev) = &                  ! Pressure levels [Pa]
       (/100000.,85000.,70000.,50000.,40000., &
          30000.,25000.,20000.,15000.,10000., &
           7000., 5000., 3000., 2000., 1000. /)

  type t_obserr_amv
    integer                 :: satids(nsat) = -1        ! Satellite ids
    integer                 :: meths (nm)   = -1        ! Retrieval methods
    real(wp)                :: lat_N        =   20._wp  ! Latitude TR <> NH
    real(wp)                :: lat_S        =  -20._wp  ! Latitude TR <> SH
    real(wp)                :: err  (nlev)  = -999._wp  ! Obs.error on levels
    real(wp)                :: scale(3)     =   1.0_wp  ! Scale factor NH/TR/SH
  end type t_obserr_amv

  integer                     :: nerr = 0               ! Number of entries
  type(t_obserr_amv), pointer :: amv_obserr(:) => NULL()
!------------------------------------------------------------------------------
  !----------------------------------------------
  ! Satellite specific height correction for AMVs
  !----------------------------------------------
  integer            :: verbose      = 2        ! Verbosity level
  real(wp)           :: p_const_hcor = 0._wp    ! Constant correction above

  integer, parameter :: nhlev = 30
  type t_amv_hgtcor
    integer                     :: satid       = -1      ! Satellite ID
    integer                     :: method      = -1      ! Retrieval method
    real(wp)                    :: lat_nb      =  90._wp ! Northern bound
    real(wp)                    :: lat_sb      = -90._wp ! Southern bound
    type(t_amv_hgtcor) ,pointer :: next        => NULL() ! link(ed list)
    real(wp)                    :: lev (nhlev) =  0._wp  ! Levels [Pa]
    real(wp)                    :: bias(nhlev) =  0._wp  ! Bias   [Pa]
  end type t_amv_hgtcor

  type(t_amv_hgtcor), pointer :: amv_hgtcor => NULL ()
!==============================================================================
contains
!==============================================================================
  subroutine process_amv (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
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
    integer                :: i, j, k, l, n ! Indices
    integer                :: tsk           ! task (local copy)
    real(sp)      ,pointer :: lv (:)        ! unpacked levels
    integer       ,pointer :: ty (:)        ! unpacked types
    integer                :: iii, ioi      !++ work around NEC SX bug
    real(wp)               :: scale_amv_obserr = 1._wp
    real(wp)               :: err
    real(wp)               :: d, w          ! offset, width parameters
    integer(i8)            :: icol          ! columns to retrieve
    real(wp)               :: z             ! height(gpm) for p.-interp.
    real(wp)               :: wzp           ! weight for p.-interpolation
    real(sp)               :: plev          ! pressure from height
    !--------------------------
    ! Weights of Simpson's rule
    !--------------------------
    real(mp), parameter    :: w3_2  = 2 / 3._mp
    real(mp), parameter    :: w3_1  = (1._mp - w3_2) / 2        ! sum(w3)=1 @ mp
    real(mp), parameter    :: w3(3) = (/ w3_1, w3_2, w3_1 /)
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
        + TSK_SHRINK       &! release unused obs. in report
    ))
    if (tsk == 0) return
    !===================================
    ! tsk == TSK_INIT:
    ! initialisation of module variables
    ! read AMV specific namelists
    !===================================
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_nml_amv_obs
      call read_nml_obserr_amv ()
      call read_nml_amv_heightcorr ()
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
      spot% char      = CHR_ID
      if (spot% fwmodel /= 0) &
          spot% char  = CHR_LIN     ! Linear operator for layer winds
      if (tsk == 0) return
    endif

    !===========
    ! PSAS-space
    !===========
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
!print *, "###TSK_SETUP_FUL0:",spot%ident,spot%stret,spot%fwmodel,spot%params
      select case (spot% fwmodel)
      case (0)
        call new_int (obs% o, spot, spot% o%n)
        obs%   o% t_int (spot%i%i+1:spot%i%i+spot%i%n) = varno_oq (&
          obs% o% varno (spot%o%i+1:spot%o%i+spot%o%n) )
!NEC$ novector
        obs%   o% lev   (spot%i%i+1:spot%i%i+spot%i%n) = log ( &
          real(obs% o% body  (spot%o%i+1:spot%o%i+spot%o%n)% plev, kind=wp))
      case (1,2,3)
        !----------------------------------------------------------
        ! n equidistant support points (levels) for linear operator
        !----------------------------------------------------------
        w = spot% params(2)
        d = spot% params(1) - 0.5_wp * w
        !------------------------
        ! Safeguard log(p+delta):
        ! Do we also need to adjust check_domain???
        !------------------------
        if (minval (obs% o% body (spot%o%i+1:spot%o%i+2)% plev) + d <= 0) then
          call decr_rpt_use (spot, CHK_HEIGHT, use=STAT_DISMISS, comment="layer")
          return
        end if
        n = spot% fwmodel
        call new_int (obs% o, spot, 2 * n)      ! Note: spot% o%n = 2
        w = w / max (n-1,1)
        k = 0
        do i = 1, n
          obs%   o% t_int (spot%i%i+k+1:spot%i%i+k+2) = varno_oq (&
            obs% o% varno (spot%o%i  +1:spot%o%i  +2) )
          obs%   o% lev   (spot%i%i+k+1:spot%i%i+k+2) = log ( &
            obs% o% body  (spot%o%i  +1:spot%o%i  +2)% plev + d + (i-1)*w)
          k = k + 2
        end do
      case default
        call finish('process_amv','invalid forward operator id')
      end select
      tsk = tsk - TSK_SETUP_FUL0
      if (tsk == 0) return
    endif

    !================================================================
    ! tsk == TSK_SETUP_COLS:
    ! determine the model column indices required by AMV observations
    ! specify the required fields (wind components)
    !================================================================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      icol = COL_UV
      if (obs% o% body (spot%o%i+1)% lev_typ == VN_HEIGHT) icol = icol + COL_GEO
      call idx_init (      &
            spot% col% c,  &! <-  column descriptor
            spot% col% h,  &!  -> interpolation coefficients
            obs% o% mc,    &! <-> model column descriptors
            icol,          &! <-  fields required
            0,             &! <-  tracers required
            atm% grid,     &! <-  model grid
            spot% i_time,  &! <-  time slot
            spot% w_time   )! <-  time interpolation weight

      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)

      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif

    !======================================================
    ! tsk == TSK_SETUP_FULL:
    ! setup description of PSAS-space
    ! interpolate pressure from height using the bg-profile
    !======================================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      if (dace% pe == obs% o% pe .and. &
          z2p_amv  == 2                ) then
        j = sum (maxloc (spot% col% h% w))
        k = spot% col% h% imc (j,1)
        do i = spot%o%i+1, spot%o%i+spot%o%n
          !-----------------------------------------------------------
          ! interpolate pressure if height is the independent quantity
          !-----------------------------------------------------------
          if (obs% o% body (i)% lev_typ == VN_HEIGHT) then
            n = size(cols% col(k)% geo)
            plev    = -1._wp
            z       = obs% o% olev (i) * gacc
            if (z <= cols% col(k)% geo(n)) then
              plev = exp (cols% col(k)% p(n))
            elseif (z >= cols% col(k)% geo(1)) then
              plev = exp (cols% col(k)% p(1))
            else
              do l = 2, size(cols% col(k)% geo)
                if (z > cols% col(k)% geo(l)) then
                  wzp     = (z                      - cols% col(k)% geo(l)) / &
                            (cols% col(k)% geo(l-1) - cols% col(k)% geo(l))
                  plev = exp (      wzp  * cols% col(k)% p(l-1) + &
                             (1._wp-wzp) * cols% col(k)% p(l  )   )
                  exit
                endif
              end do
            endif
            obs% o% body (i)% plev    = plev
            !---------------------------------------------------------
            ! currently we set 'olev' to 'plev' as some parts
            ! of the code still assume 'olev' is a pressure coordinate
            !---------------------------------------------------------
            obs% o% body (i)% lev_typ = VN_P
            obs% o% olev (i)          = plev
            spot% ps                  = plev
          endif
        end do
      endif

      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) return
    endif

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
      ty => obs% o% varno (spot% o% i + 1 : spot% o% i + spot% o% n)
      lv => obs% o% body  (spot% o% i + 1 : spot% o% i + spot% o% n)% plev
      n = spot% o% n
      !===============================
      ! setup AMV observational errors
      !===============================
      if (obs% o% pe == dace% pe) then
        k = obs% R% ia (spot% o% i + 1)
        do i=1,n
          iii = spot% o% i + i
          obs% R% ia (iii) = k
          select case (ty(i))
          case (VN_U,VN_V)
            err = obserr (spot% ident, spot% stret, spot% col% c% dlat, real(lv(i),wp))
            if (err > 0._wp) then
               obs% R% packed(k) = err**2
            else
               ! Fallback to old treatment of AMV observation error
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! temporary modifications for obs.error scaling
!
 if ((spot% statid .eq. 'TERRA') .or. (spot% statid .eq. 'AQUA')) then
    scale_amv_obserr = 0.8_wp
 else
    if (abs(spot% col% c% dlat) .ge. 20._wp)  then
       scale_amv_obserr = 0.7_wp
    else
       scale_amv_obserr = 0.5_wp
    endif
 endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            obs% R% packed(k) = (scale_amv_obserr) ** 2               &
              * obs_err (OT_SATOB,                                    &
                         spot%hd% codetype, VN_U, real(lv(i),wp),0._sp) ** 2
            end if

          case default
            call finish('process_amv','invalid observation type')
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

    !=========================================
    ! tsk == TSK_Y:
    ! run nonlinear operator (identity so far)
    !=========================================
    if (iand (TSK_Y,tsk) /= 0) then
if (spot%o%n /= 2) call finish ("process_amv","tsk_y:n/=2")
      k  = spot% o% n
      iii= spot% i% i
      ioi= spot% o% i
      select case (spot% fwmodel)
      case (0,1)
        y%x (ioi+1:ioi+k) = xi%x (iii+1:iii+k)
      case (2) ! Trapezoidal rule
        i = 1
!       do i = 1, k, 2
           j = iii + 2*i - 1
           y%x(ioi+i  ) = (xi%x(j  ) + xi%x (j+2)) * 0.5_wp   ! u
           y%x(ioi+i+1) = (xi%x(j+1) + xi%x (j+3)) * 0.5_wp   ! v
!       end do
      case (3) ! Simpson's rule
        i = 1
!       do i = 1, k, 2
           j = iii + 3*i - 2
           y%x(ioi+i  ) = w3(1) * (xi%x(j  ) + xi%x(j+4)) + w3(2) * xi%x(j+2) ! u
           y%x(ioi+i+1) = w3(1) * (xi%x(j+1) + xi%x(j+5)) + w3(2) * xi%x(j+3) ! v
!       end do
      end select
      tsk = tsk - TSK_Y
      if (tsk == 0) return
    endif

    !==================================
    ! tsk == TSK_K:
    ! set up H matrix (identity so far)
    !==================================
    if (iand (TSK_K,tsk) /= 0) then
      select case (spot% fwmodel)
      case (0,1)
        k = obs% H% ia (spot% i% i + 1)
        do i=1,spot% i% n                               ! rows = columns
          iii = spot% i% i + i                          !+++ workaround
          ioi = spot% o% i + i                          !+++ for NEC SX
          obs% H% ia (iii) = k                          ! column index
          obs% H% ja (k) = ioi                          ! row index
          obs% H% packed (k) = 1._mp                    ! coefficient
          k = k + 1
        end do
        obs%   H% ia (spot% i% i + spot% i% n + 1) = k  ! column index
        obs%   yi% x (spot% o% i+1:spot% o% i+spot% o% n) = &
          obs% xi% x (spot% i% i+1:spot% i% i+spot% i% n)
      case (2) ! Trapezoidal rule
if (spot%i%n /= 4) call finish ("process_amv","tsk_k:n/=4")
        k = obs% H% ia (spot% i% i + 1)
        do j=1,spot% i% n                             ! columns
          obs% H% ia (spot% i% i +j) = k              ! column index
          i = mod (j-1,2) + 1                         ! row
          obs% H% ja     (k) = spot% o% i + i         ! row index
          obs% H% packed (k) = 0.5_mp                 ! coefficient
          k = k + 1
        end do
        obs% H% ia (spot% i% i + spot% i% n + 1) = k  ! column index
        k  = 2 ! ==spot% o% n
        iii= spot% i% i
        ioi= spot% o% i
        i  = 1
!       do i = 1, k, 2
           j = iii + 2*i - 1
           obs%yi%x(ioi+i  ) = (obs%xi%x(j  ) + obs%xi%x (j+2)) * 0.5_wp   ! u
           obs%yi%x(ioi+i+1) = (obs%xi%x(j+1) + obs%xi%x (j+3)) * 0.5_wp   ! v
!       end do
      case (3) ! Simpson's rule
if (spot%i%n /= 6) call finish ("process_amv","tsk_k:n/=6")
        k = obs% H% ia (spot% i% i + 1)
        do j=1,spot% i% n                             ! columns
          obs% H% ia (spot% i% i +j) = k              ! column index
          i = mod (j-1,2) + 1                         ! row
          l =     (j-1)/2 + 1
          obs% H% ja     (k) = spot% o% i + i         ! row index
          obs% H% packed (k) = w3(l)                  ! coefficient
          k = k + 1
        end do
        obs% H% ia (spot% i% i + spot% i% n + 1) = k  ! column index
        k  = 2 ! ==spot% o% n
        iii= spot% i% i
        ioi= spot% o% i
        i  = 1
!       do i = 1, k, 2
           j = iii + 3*i - 2
           obs%yi%x(ioi+i  ) = w3(1)*(obs%xi%x(j  ) + obs%xi%x(j+4)) &
                             + w3(2)* obs%xi%x(j+2)                  ! u
           obs%yi%x(ioi+i+1) = w3(1)*(obs%xi%x(j+1) + obs%xi%x(j+5)) &
                             + w3(2)* obs%xi%x(j+3)                  ! v
!       end do
      end select
      tsk = tsk - TSK_K
      if (tsk == 0) return
    endif

    !===============================
    ! abort if invalid task is given
    !===============================
    call finish('process_amv','unknown task: '//tsk_name(tsk))

  end subroutine process_amv
!==============================================================================
  subroutine read_amv_bufr (bufr, spt, obs, lkeep, nkeep, cc)
  type (t_bufr) ,intent(inOUT)        :: bufr  ! BUFR record to decode
  type (t_spot) ,intent(inout)        :: spt   ! meta information to set
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type to set
  logical       ,intent(out)          :: lkeep ! accept observation ?
  integer       ,intent(out)          :: nkeep ! number of accepted obsvs.
  integer       ,intent(in) ,optional :: cc    ! part of year ccyy
  !=====================================
  ! Read AMV observations from BUFR data
  !=====================================
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
    integer         :: ir         ! record nr. (subset or loopcount)
    integer         :: irend      ! last record
    logical         :: fy_x       ! .true. for chinese satellite,different BUFR
    integer         :: is         ! sub-set    index
    integer         :: ie         ! entry in sub-set
    integer         :: id         ! descriptor index

    type (t_amv)    :: s          ! AMV observation read
    type (t_spot)   :: spt0, spti ! observation meta data
    logical         :: unused

    lkeep = .false.
    nkeep = 0
    spt0         = spt
    spt0% hd% id = spt% hd% id - 1
    !---------------------
    ! loop over data, copy
    !---------------------
    irend = max(1,bufr% sec3% num_subsets)
    is   = 1
    ir   = 0
    fy_x = .false.
    do
      ir = ir + 1
      if (bufr% sec3% num_subsets > 1) is = ir
      spti         = spt0
      spti% hd% id = spt0% hd% id + ir
      if (bufr% sec3% num_subsets>0) spti% hd% subset = ir
      !----------------------------------------------------------
      ! construct 'empty' Report data structure
      ! create additional 'reports' for subsequent subset members
      !----------------------------------------------------------
      call construct_amv (s)
      if (ir > 1) call check_report_0 (spti% use, spti% hd, 1)
      !--------------------------------
      ! loop over subsets or loopcounts
      !--------------------------------
      unused = .false.
      if (.not.fy_x) ie = 0
      do
        ie = ie + 1
        if (ie > bufr% nbufdat(is)) exit
        !---------------------------------------
        ! decode real/character datum, mnemonics
        !---------------------------------------
        ival  = bufr% ibufdat (is,ie)
        id    = bufr% idescidx(is,ie)
        if (bufr% is_char(id) == 0 .and. ival == inv_bufr) cycle
        itype = bufr% itype(id)
        ymnem = bufr% ymnem(id)
        rval  = rvind
        if (bufr% is_char(id) == 0) then
          IF (ival /= inv_bufr) &
            rval = ival * bufr% scale(id)
        else
!         call bufr3_get_character_value (bufr% ihandle, text, is, ie)
        endif

        select case (ymnem)
        !------------------------------
        ! process the following entries
        !------------------------------
        case ('MCORME')
          spti% corme = ival
          if (ival/=0) call finish ('read_amv_bufr',&
                                    ymnem//' cannot handle corrections')
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
        case ('MGG')    ! hour
          hh   = ival
        case ('NGG')    ! minute
          mi   = ival
        case ('MI1I2')  ! satellite identifier
          spti% ident     = ival
          spti% hd% satid = ival
!         fy_x = ival==513
        case ('MSIDP')  ! SATELLITE INSTRUMENT USED IN DATA PROCE.
          if (ival<0.or.ival>15) ival = -1 ! enforce missing value
          spti% sttyp    = ival
        case ('MSDWCM') ! satellite derived wind computation meth.
          spti% stret    = ival
        case ('MSZAN')  ! SATELLITE ZENITH ANGLE
          spti% stzen    = ival
        case ('MLASE')  ! LAND/SEA QUALIFIER
          spti% stlsf    = ival
        case ('MHEAM')  ! HEIGHT ASSIGNMENT METHOD
          spti% stclf    = ival
        case ('MOGC')   ! ORIGINATING/GENERATING CENTRE
          spti% hd% center = ival
        case ('MPN')    ! pressure (vert.location)
          call set_datum (s% p  ,rval ,spti% corme)
        case ('MPPP')   ! pressure
          call set_datum (s% p  ,rval ,spti% corme)
        case ('MTN')    ! TEMPERATURE/DRY BULB TEMPERATURE
          call set_datum (s% t  ,rval ,spti% corme)
        case ('NDNDN')  ! wind direction
          call set_datum (s% dd ,rval ,spti% corme)
        case ('NFNFN')  ! wind speed
          call set_datum (s% ff ,rval ,spti% corme)
        case ('MSCCF')  ! SATELLITE CHANNEL CENTRE FREQUENCY
          call set_datum (s% f ,rval  ,spti% corme)
!       case ('MPCCO')  ! per cent confidence
!         if (spti% hd% dbkz == 1704) then
!           if (spti% pcc < 0) then
!             spti% pcc = ival
!             s%ff% pcc = spti% pcc
!           endif
!         endif
!       case ('MPCCO1') ! per cent confidence
!         if (spti% hd% dbkz == 1705) then
!           if (spti% pcc < 0) then
!             spti% pcc = ival
!             s%ff% pcc = spti% pcc
!           endif
!         endif

        !--------------------------------------------------------------
        ! Interpretation of PCC (per cent confidence) flags :
        !
        ! pcc_conv = 1: as used for BUFR reports from Pegasus data base
        !
        ! pcc_conv = 2: to be used for SKY data base
        !               (compatible with NetCDF interface)
        !--------------------------------------------------------------
        case ('MPCCO')  ! per cent confidence

          if     (pcc_conv == 2) then
            select case (spti% hd% center)
            case (WMO0_NOAA, WMO0_NASA, WMO0_CIMSS, WMO0_JMA)
              if (spti% pcc < 0) then
                spti% pcc = ival
                s%ff% pcc = spti% pcc
              endif
            end select

          elseif (pcc_conv == 1) then
            select case (spti% hd% center)
            case (WMO0_EUMET, WMO0_JMA, WMO0_RSMC)
              if (spti% pcc < 0) then
                spti% pcc = ival
                s%ff% pcc = spti% pcc
              endif
            end select
          endif

        case ('MPCCO0')  ! per cent confidence

          if     (pcc_conv == 2) then
            select case (spti% hd% center)
            case (WMO0_EUMET, WMO0_RSMC)
              if (spti% pcc < 0) then
                spti% pcc = ival
                s%ff% pcc = spti% pcc
              endif
            end select
          endif

        case ('MPCCO1')  ! per cent confidence

          if     (pcc_conv == 1) then
            select case (spti% hd% center)
            case (WMO0_NOAA, WMO0_NASA, WMO0_CIMSS)
              if (spti% pcc < 0) then
                spti% pcc = ival
                s%ff% pcc = spti% pcc
              endif
            end select
          endif

        !------------------------------------
        ! exit loop in case of satellite FY_X
        !------------------------------------
        case ('MEDRE')  ! descr.repl.factor
!         if (fy_x) irend = ival
        case ('Lcnt000')
!         if (fy_x) then
!           if (ival == 0) spt0 = spti
!           if (ival >=  ir) exit
!         endif
#ifdef  CHECKCODES
        !-----------------------
        ! check for unused codes
        !-----------------------
        !-----------------------------------------
        ! ignore the following entries (in SATOBS)
        !-----------------------------------------
        case ('MOBITQ') ! OVERALL QUALITY BITS
        case ('MADDF')  ! ASSOCIATED FIELD SIGNIFICANCE
        case ('YCCCC')  ! ICAO LOCATION INDICATOR
        case ('YDDDD')  ! SHIP OR MOBILE LAND STATION IDENTIFIER
        case ('YXXNN')  ! AIRCRAFT FLIGHT NUMBER
        case ('MADDF0') ! ASSOCIATED FIELD SIGNIFICANCE
        case ('MVTSO')  ! VERTICAL SOUNDING SIGNIFICANCE
        case ('Loop000')! Start of Loop - 101000
        case ('MDREP')  ! DELAYED DESCRIPTOR REPLICATION FACTOR
        !------------------------------------
        ! additionally ignore in EUMETSAT AMV
        !------------------------------------
        case ('MSACL')  ! SATELLITE CLASSIFIKATION
        case ('MSSXD')  ! SEGMENT SIZE AT NADIR IN X DIRECTION
        case ('MSSYD')  ! SEGMENT SIZE AT NADIR IN Y DIRECTION
        case ('MSEC')   ! SECOND
        case ('MSCBW')  ! SATELLITE CHANNEL BAND WIDTH
        case ('MCOCT')  ! COLDEST CLUSTER TEMPERATURE
        case ('Lcnt001':'Lcnt011') ! Loop Counter
        case ('Loop001':'Loop011') ! Start of Loop
        case ('MOGC0'  :'MOGC8'  ) ! ORIGINATING/GENERATING CENTRE
        case ('MSGAPP0':'MSGAPP7') ! GENERATING APPLICATION
        case ('MMAQC0' :'MMAQC1' ) ! MANUAL/AUTOMATIC QUALITY CONTROL
        case ('MNOCT0' :'MNOCT1' ) ! NOMINAL CINFIDENCE THRESHOLD
        case ('MTISI1') ! TIME SIGNIFICANCE
        case ('MGG1')   ! HOUR
        case ('NGG0')   ! MINUTE
        case ('MSEC0')  ! SECOND
        case ('MTISI2') ! TIME SIGNIFICANCE
        case ('MGG2')   ! HOUR
        case ('NGG1')   ! MINUTE
        case ('MSEC1')  ! SECOND
        case ('NDNDN0') ! WIND DIRECTION
        case ('NFNFN0') ! WIND SPEED
        case ('MDPIN')  ! DATA PRESENT INDICATOR
        case ('MSGAPP') ! GENERATING APPLICATION
        case ('MMAQC')  ! MANUAL/AUTOMATIC QUALITY CONTROL
        case ('MNOCT')  ! NOMINAL CINFIDENCE THRESHOLD
        case ('MHEAM0') ! HEIGHT ASSIGNMENT METHOD
        case ('MPN0')   ! PRESSURE (VERT.LOCATION)
        case ('MTRCM')  ! TRACE CORRELATION METHOD
        !--------------------------------
        ! additionally ignore in GEOS AMV
        !--------------------------------
        case ('MTISI':'MTISI0')  ! TIME SIGNIFICANCE
        !---------------------------------
        ! additionally ignore in MODIS AMV
        !---------------------------------
        case ('YSUPL')           ! CHARACTERS
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
      end do
      !--------------------------------------
      ! in case of unknown codes:
      ! print suspicious BUFR record and exit
      !--------------------------------------
      if (unused) then
        call bufr_print_sections (bufr)
        call bufr_print_subset   (bufr, is)
        call finish ('read_amv_bufr','code(s) not implemented')
      endif

      !------------------------------------
      ! use station number for station name
      !------------------------------------
      spti% statid = satname (spti% ident)

      !----------------
      ! standard checks
      !----------------
      lkeep    = .true.
      spti% ps = s% p% o
      call init_time (spti% actual_time, yyyy, mo, dd, hh, mi)
      call check_report_1 (spti)
      call check_wind_method (spti, s)

      if (lkeep) then
        !-------------------------
        ! check for double entries
        !-------------------------
!       ....
!
        !--------------------------------------------------------
        ! insert in list if most important information is present
        !--------------------------------------------------------
        if (lkeep) call check_store_amv (s, spti, obs, lkeep)
        if (lkeep) then
          nkeep = nkeep + 1
        endif
      endif
      if (ir >= irend) exit
    end do
    lkeep = nkeep > 0

  end subroutine read_amv_bufr
!------------------------------------------------------------------------------
  subroutine check_amv_qual (gobs)
  type(t_obs) ,intent(inout) :: gobs
  !------------------------------------
  ! check AMV quality index (postponed)
  !------------------------------------
    integer :: is

    if (chk_amv_post) then
      do is = 1, gobs% n_spot
        if (gobs% spot(is)% hd% obstype == OT_SATOB) &
          call check_amv_indx (gobs% spot(is))
      end do
    endif

  end subroutine check_amv_qual
!------------------------------------------------------------------------------
  subroutine check_amv_indx (spt)
  type(t_spot) ,intent(inout) :: spt
  !---------------------------------------
  ! check AMV quality index (spot by spot)
  !---------------------------------------

    type(t_amv_select), pointer :: select
    integer                :: il         ! level index (1=low,2=med,3=high)
    integer                :: ih         ! hemisphere index (1=NH,2=TR,3=SH)

    select => amv_select
    do
      if (.not.associated(select)) then
        !------------------------------------------------------------
        ! don't use observation in case of no matching namelist group
        !------------------------------------------------------------
        spt% fwmodel = 0
        spt% params  = 0._wp
        call decr_rpt_use (spt ,CHK_NOTUSED)     ! No matching selection set
        exit
      endif

      if (all(select% meths  == -1)   .or. &
          any(select% meths  == spt% stret)) then

        if (all(select% satids == -1)   .or. &
            any(select% satids == spt% ident)) then

          il = 2
          if (spt% ps > select% p_low)  il = 1
          if (spt% ps < select% p_high) il = 3
          ih = 2
          if (spt% col% c% dlat >= select% lat_N) ih = 1
          if (spt% col% c% dlat <= select% lat_S) ih = 3

          if (spt% pcc       <= select% qimin (il,ih)) &
            call decr_rpt_use (spt ,CHK_QI)

          !---------------------------------------------------------------
          ! Set associated forward operator (in-situ or model layer winds)
          !---------------------------------------------------------------
          spt% fwmodel = select% fwmodel
          spt% params  = select% params((2*il-1):2*il)
          exit
        endif
      endif
      select => select% next
    end do

  end subroutine check_amv_indx
!------------------------------------------------------------------------------
  subroutine check_amv_surf (spt)
  type(t_spot) ,intent(inout) :: spt

    type(t_amv_select), pointer :: select
    integer                     :: ih       ! hemisphere index

    select => amv_select
    do
      if (.not.associated(select)) then
        exit
      endif
      if (spt% stret >= 1            .and.  &
         (all(select% meths  == -1)   .or.  &
          any(select% meths  == spt% stret))) then
        if (all(select% satids == -1)   .or. &
            any(select% satids == spt% ident)) then
          ih = 2
          if (spt% col% c% dlat >= select% lat_N) ih = 1
          if (spt% col% c% dlat <= select% lat_S) ih = 3
          if (select% seaonly(ih) .and. &
              spt% sl_bg > 0._wp) call decr_rpt_use (spt, CHK_SURF)
          exit
        endif
      endif
      select => select% next
    end do
  end subroutine check_amv_surf
!==============================================================================
  subroutine check_wind_method (spot, s)
    type(t_spot) ,intent(inout) :: spot ! observation meta data
    type(t_amv)  ,intent(in)    :: s    ! AMV observation read
    !-----------------------------------------------------------------------
    ! Code Figure 'MSDWCM' (satellite derived wind computation meth.)
    ! 0 Not used
    ! 1 Wind derived from cloud motion observed in the infrared channel
    ! 2 Wind derived from cloud motion observed in the visible channel
    ! 3 Wind derived from motion observed in the water vapour channel
    !   (only the above are used up to now)
    ! 4 Wind derived from motion observed in a combination of spectral chans
    ! 5 Wind derived from motion obs. in the water vapour ch. in clear air
    ! 6 Wind derived from motion observed in the ozon channel
    ! 7 Wind derived from cloud motion observed in water vapour channel
    !   (cloud or clear air not specified)
    !
    ! Add satellite/frequency digit (n*100) to spot% stret
    !-----------------------------------------------------------------------
    real(wp), parameter :: epsilon = 0.00001e+15_wp
    integer             :: digits

    if (spot% stret >= 1) then
      digits = 0
      select case (spot% ident)
      case (55:57,70)   ! Meteosat-8, Meteosat-9, Meteosat-10, Meteosat-11
        if     (abs (s%f%o - 0.4721142137D15) < epsilon) then
           digits = 100                 ! VIS1: 0.6 um
        elseif (abs (s%f%o - 0.399723D15)     < epsilon) then
           digits = 300                 ! VIS2 (high resolution): 0.75 um
        elseif (abs (s%f%o - 0.3701142098D15) < epsilon) then
           digits = 200                 ! VIS3: 0.8 um
        elseif (abs (s%f%o - 0.4796679874D14) < epsilon) then
           digits = 100                 ! WV1 : 6.2 um
        elseif (abs (s%f%o - 0.4078810064D14) < epsilon) then
           digits = 200                 ! WV2 : 7.3 um
        elseif (abs (s%f%o - 0.344589D14)     < epsilon) then
           digits = 100                 ! IR1 : 8.7 um
        elseif (abs (s%f%o - 0.310344D14)     < epsilon) then
           digits = 200                 ! IR2 : 9.7 um
        elseif (abs (s%f%o - 0.2775860013D14) < epsilon) then
           digits = 300                 ! IR3 : 10.8 um
        else
!         WRITE(0,*) '%E SATOBFREQ: channel frequency for',spot% ident,&
!                    ' has not been recognized',s%f%o
        endif
      case (254:259) ! GOES-10, GOES-11, GOES-12, GOES-13, GOES-14, GOES-15
        if     (abs (s%f%o - 0.461538D15) < epsilon) then
           digits = 100                 ! VIS1 (0.65um)
        elseif (abs (s%f%o - 0.405405D14) < epsilon) then
           digits = 100                 ! WV1 (7.4um)
        elseif (abs (s%f%o - 0.428571D14) < epsilon) then
           digits = 200                 ! WV2 (7.0um)
        elseif (s%f%o > 0.441176D14-epsilon .and. s%f%o < 0.461219D14+epsilon) then
           digits = 300                 ! WV3 (6.5 till 6.8um)
        elseif (abs (s%f%o - 0.280374D14) < epsilon) then
           digits = 100                 ! IR1 (10.7um)
        elseif (abs (s%f%o - 0.768699D14) < epsilon) then
           digits = 200                 ! IR2 (3.9um)
        else
           WRITE(0,*) '%E SATOBFREQ: channel frequency for',spot% ident,&
                      ' has not been recognized',s%f%o
        endif
      case (270:273) ! GOES-16, GOES-17, GOES-18, GOES-19
        if     (abs (s%f%o - 0.468426D15) < epsilon) then
           digits = 100                 ! VIS1 (0.64um)
        elseif (abs (s%f%o - 0.408437D14) < epsilon) then
           digits = 100                 ! WV1 (7.3um)
        elseif (abs (s%f%o - 0.431356D14) < epsilon) then
           digits = 200                 ! WV2 (7.0um)
        elseif (abs (s%f%o - 0.484317D14) < epsilon) then
           digits = 300                 ! WV3 (6.2um)
        elseif (abs (s%f%o - 0.267672D14) < epsilon) then
           digits = 100                 ! IR1 (11.1um)
        elseif (abs (s%f%o - 0.768699D14) < epsilon) then
           digits = 200                 ! IR2 (3.9um)
        else
           WRITE(0,*) '%E SATOBFREQ: channel frequency for',spot% ident,&
                      ' has not been recognized',s%f%o
        endif
      case (173, 174) ! Himawari-8/9
        if     (abs (s%f%o - 0.4684257D15) < epsilon) then
           digits = 100                 ! VIS1: 0.645 um
        elseif (abs (s%f%o - 0.483536D14) < epsilon) then
           digits = 100                 ! WV1 : 6.25 um
        elseif (abs (s%f%o - 0.434481D14) < epsilon) then
           digits = 200                 ! WV2 : 6.95 um
        elseif (abs (s%f%o - 0.410674D14) < epsilon) then
           digits = 300                 ! WV3 : 7.35 um
        elseif (abs (s%f%o - 0.288261D14) < epsilon) then
           digits = 100                 ! IR1 (10.45um)
        else
           WRITE(0,*) '%E SATOBFREQ: channel frequency for',spot% ident,&
                      ' has not been recognized',s%f%o
        endif
      case default
!       WRITE(0,*) '%E SATOBFREQ: channel frequency for',spot% ident,&
!                  ' has not been recognized',s% f% o
      end select
      spot% stret = spot% stret + digits
    else
!     WRITE(0,*) '%E SATOBFREQ: invalid value for "stret"',spot% stret
    endif
  end subroutine check_wind_method
!==============================================================================
  subroutine check_store_amv (amv, spot, obs, lkeep)
  type(t_amv)  ,intent(inout) :: amv   ! AMV level information
  type(t_spot) ,intent(inout) :: spot  ! meta data of this observation
  type(t_obs)  ,intent(inout) :: obs   ! data of all observations
  logical      ,intent(out)   :: lkeep ! flag: keep or reject
  !==================================================
  ! check data for consistency
  ! store accepted data in observation data structure
  !==================================================

    type (t_spot) ,pointer :: s
    integer ,parameter     :: n = 2      ! u&v
    integer                :: id

    lkeep = .true.
    if (spot% col% c% dlon == invalid .or. &
        spot% col% c% dlat == invalid .or. &
        amv%  p% o         <= 0       .or. &
        amv%  p% o         >= 110000  .or. &
        amv%  p% qc        /= 0       .or. &
        amv% ff% qc        /= 0       .or. &
        amv% dd% qc        /= 0            ) lkeep = .false.

    !---------------------------------------------------------------------
    ! Quality Index Check
    ! depends on: QI, satellite id, latitude, pressure, computation method
    !---------------------------------------------------------------------
    if (lkeep) then

      spot% ps = amv% p% o
      if (chk_amv_cdfin) call check_amv_indx (spot)

      if (spot% use% state <= STAT_DISMISS) lkeep = .false.
    end if

    !--------------------
    ! store data into OBS
    !--------------------
    if (lkeep) then
      call new_spot (obs, 1, set_id=.true.)
      s => obs% spot (obs% n_spot)
      id    = s% id
      s     = spot
      s% id = id
!     s% int_type  = ITY_ICOL
      s% col% nlev = 1
!     s% cost      = 1._wp
!     s% char      = CHR_ID
      call new_obs (obs, n, s)
!     s% nr        = n
!     call new_int (obs, s, n)
      obs%   body (s%o%i+1)          = amv% ff
      obs%   body (s%o%i+1) %o       = amv% ff% o * (-sin (d2r * amv% dd% o))
      obs%   body (s%o%i+1) %lev_typ = VN_P
      obs%   body (s%o%i+1) %plev    = amv% p% o
      obs%   olev (s%o%i+1)          = amv% p% o
      obs%  varno (s%o%i+1)          = VN_U
      obs%   body (s%o%i+2)          = amv% ff
      obs%   body (s%o%i+2) %o       = amv% ff% o * (-cos (d2r * amv% dd% o))
      obs%   body (s%o%i+2) %lev_typ = VN_P
      obs%   body (s%o%i+2) %plev    = amv% p% o
      obs%   olev (s%o%i+2)          = amv% p% o
      obs%  varno (s%o%i+2)          = VN_V
      !---------------------------------------------------------
      ! currently we set 'olev' to 'plev' as some parts
      ! of the code still assume 'olev' is a pressure coordinate
      ! only for z2p_am==2 this is done later in TSK_SETUP_FULL
      !---------------------------------------------------------
      if (z2p_amv     == 2       .and. &
          amv% p% src == SRC_DOK .and. &
          amv% z% qc  == 0             ) then
        obs% olev (s%o%i+1)          = amv% z% o
        obs% olev (s%o%i+2)          = amv% z% o
        obs% body (s%o%i+1) %lev_typ = VN_HEIGHT
        obs% body (s%o%i+2) %lev_typ = VN_HEIGHT
      endif
      call set_xuv (s)
    else
      call decr_rpt_use (spot ,CHK_INSDAT, STAT_DISMISS)
    endif
  end subroutine check_store_amv
!==============================================================================
  pure subroutine construct_amv (s)
  type (t_amv) ,intent(out) :: s
!#if defined(__SX__)
!   ! Default initialisation does not work with sxf90 rev.360 (SX-6)
!   s = t_amv (inv_datum,inv_datum,inv_datum,inv_datum,inv_datum,inv_datum)
!#endif
    s% total %mn =''
    s% p     %mn ='p'
    s% t     %mn ='t'
    s% ff    %mn ='ff'
    s% dd    %mn ='dd'
  end subroutine construct_amv
!==============================================================================
  subroutine read_nml_amv_obs
  !=====================================
  ! read AMV specific namelist /AMV_OBS/
  !=====================================
    character(len=8) :: satids (nsat)! satellite ids (mnemonics)
    integer          :: meths  (nm)  ! wind computation method
    integer          :: qimin  (3,3) ! minimum quality index
    real(wp)         :: p_low        ! lower  pressure bound [hPa]
    real(wp)         :: p_high       ! higher pressure bound [hPa]
    real(wp)         :: lat_N        ! latitude separating TR,NH
    real(wp)         :: lat_S        ! latitude separating TR,SH
    logical          :: seaonly(3)   ! use only over sea
    integer          :: fwmodel      ! forward model
    real(wp)         :: params (6)   ! parameters (offset,width)*(low,med,high)
    integer          :: pass         ! 1:apply during cdfin-read, 2:apply later

    namelist  /AMV_OBS/ satids, qimin, meths, &
                        p_low, p_high, lat_N, lat_S, seaonly, &
                        fwmodel, params, pass

    integer                     :: ierr
    logical                     :: first = .true.
    type(t_amv_select), pointer :: p
#if defined(__ibm__)
    integer                     :: ios
#endif

    !-----------------------------------------
    ! read namelist groups /AMV_OBS/ only once
    !-----------------------------------------
    if (.not. first) return

    !-------------
    ! set defaults
    !-------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') ' Reading namelist /AMV_OBS/'
      write(6,'()')
    endif

    do
      satids   = ''
      qimin    = 999    ! -99:all, 100:none else qi,lower bound
      meths    = -1
      p_low    = 700._wp
      p_high   = 400._wp
      lat_N    =  20._wp
      lat_S    = -20._wp
      seaonly  = .false.
      fwmodel  = 0      ! 0: original (in-situ obs.), n: n-point formula
      params   = -HUGE (0._wp)
      pass     = -1
      !---------------------------------
      ! read namelist, consistency check
      !---------------------------------
      if (dace% lpio) then
        call position_nml ('AMV_OBS' ,lrewind=first ,status=ierr)
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=AMV_OBS, iostat=ios)
          if(ios/=0)call finish('read_nml_amv_obs','ERROR in namelist /AMV_OBS/')
#else
          read (nnml ,nml=AMV_OBS)
#endif
        end select
      endif
      first = .false.
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !------------------------
      ! broadcast to other PE's
      !------------------------
      call p_bcast (satids,  dace% pio)
      call p_bcast (meths,   dace% pio)
      call p_bcast (qimin,   dace% pio)
      call p_bcast (p_low,   dace% pio)
      call p_bcast (p_high,  dace% pio)
      call p_bcast (lat_N,   dace% pio)
      call p_bcast (lat_S,   dace% pio)
      call p_bcast (seaonly, dace% pio)
      call p_bcast (fwmodel, dace% pio)
      call p_bcast (params,  dace% pio)
      call p_bcast (pass,    dace% pio)
      !--------------------------------
      ! Adjust default parameter values
      !--------------------------------
      if (all (params(3:4)== -HUGE (0._wp))) then
        params(3:4) = params(1:2)               ! med  same as for low
      end if
      if (all (params(5:6)== -HUGE (0._wp))) then
        params(5:6) = params(3:4)               ! high same as for med
      end if
      !------------------------------------------------
      ! Consistency checks for forward model parameters
      !------------------------------------------------
      select case (fwmodel)
      case (0)
        params = 0
      case (1)
        params(2::2) = 0        ! Enforce width == 0
      case (2:3)
        if (any (params(2::2) < 0)) &
        call finish ('read_nml_amv_obs','all width parameters must be >= 0!')
      case default
        call finish ('read_nml_amv_obs','fwmodel out of range')
      end select
      if (any (params(1::2) == -HUGE (0._wp))) &
        call finish ('read_nml_amv_obs','offset parameters must all be set!')
      if (any (params(1::2) - params(2::2)/2 < -100._wp)) &
        call finish ('read_nml_amv_obs','offset-width/2 outside allowed range!')

      p => amv_select
      allocate (amv_select)

      amv_select% satids   =  satid (satids)
      amv_select% meths    =  meths
      amv_select% qimin    =  qimin
      amv_select% p_low    =  p_low  * 100._wp
      amv_select% p_high   =  p_high * 100._wp
      amv_select% lat_N    =  lat_N
      amv_select% lat_S    =  lat_S
      amv_select% seaonly  =  seaonly
      amv_select% fwmodel  =  fwmodel
      amv_select% params   =  params * 100._wp
      amv_select% next     => p
      if (pass >= 0) chk_amv_pass = pass
      chk_amv_cdfin = iand (chk_amv_pass, 1) /= 0
      chk_amv_post  = iand (chk_amv_pass, 2) /= 0
      if (any(amv_select% satids == 0)) &
        call finish ('read_nml_amv_obs','invalid satellite name')
    end do
  end subroutine read_nml_amv_obs
!------------------------------------------------------------------------------
  subroutine clean_amv_nml
    type(t_amv_select), pointer :: p1, p2
    p1 => amv_select
    do
      if (.not. associated(p1)) exit
      p2 => p1% next
      deallocate (p1)
      p1 => p2
    end do
    nullify (amv_select)
    if (associated (amv_obserr)) deallocate (amv_obserr)
    nullify (amv_obserr)
  end subroutine clean_amv_nml
!==============================================================================
  subroutine read_amv_netcdf (ifile, i_source, obs, lkeep, nkeep,cc)
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in source-file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data type to set
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.
  integer       ,intent(in) ,optional :: cc        ! part of year ccyy

  !========================================
  ! Read AMV observations from  netCDF File
  !========================================
    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
  real    ,allocatable      :: mpn     (:)   ! pressure (vert.location) ; older mppp
  real    ,allocatable      :: nht     (:)   ! height of top of cloud (m) for misr data
  real    ,allocatable      :: ndndn   (:)   ! wind direction; integer in BUFR
  real    ,allocatable      :: nfnfn   (:)   ! wind speed
  real    ,allocatable      :: mcoct   (:)   ! coldest cluster temperature ; older mtn cloud top t
! real    ,allocatable      :: mtn     (:,:) ! (2-dim)         temperature
  integer ,allocatable      :: msidp   (:)   ! satellite instrument used in data proce.
  integer ,allocatable      :: msdwcm  (:)   ! satellite derived wind computation meth.
  real    ,allocatable      :: mszan   (:)   ! satellite zenith angle
  integer ,allocatable      :: mlase   (:)   ! land/sea qualifier
  integer ,allocatable      :: mheam   (:)   ! height assignment method
  integer ,allocatable      :: mogc    (:)   ! originating/generating centre
  integer ,allocatable      :: mcc     (:)   ! cloud type
  real    ,allocatable      :: msccf   (:)   ! satellite channel centre frequency
  integer ,allocatable      :: mdpin   (:,:) ! data present indicator for flag_table
  integer ,allocatable      :: msgapp  (:,:) ! generating application
  integer ,allocatable      :: mpcco   (:,:) ! per cent confidence: meteosat mpcco0; goes mpcco
  integer ,allocatable      :: mmaqc   (:,:) ! manual/automatic quality control
  integer ,allocatable      :: mnoct   (:,:) ! nominal confidence threshold
  !---------------------------------------
  ! NetCDF variable IDs for body   section
  !---------------------------------------
  ! variable ID's in NetCDF file for
  !    expansion  in NetCDF file BUFR- data section4
  integer              :: varid_MPN      ! pressure (vert.location)          [Pa]
  integer              :: varid_NHT      ! height of top of cloud            [m]
  integer              :: varid_NDNDN    ! wind direction                    [degree]
  integer              :: varid_NFNFN    ! wind speed                        [m/s] (real)
  integer              :: varid_MCOCT    ! coldest cluster temperature       [K]   (real)
  integer              :: varid_MTN      ! (2-dim)         temperature       [K]   (real)
  integer              :: varid_MSIDP    ! satellite instrument used in data proce.
  integer              :: varid_MSDWCM   ! satellite derived wind computation meth.
  integer              :: varid_MSZAN    ! satellite zenith angle
  integer              :: varid_MLASE    ! land/sea qualifier
  integer              :: varid_MHEAM    ! height assignment method
  integer              :: varid_MOGC     ! originating/generating centre
  integer              :: varid_MMIOGC   ! originating/generating centre (new)
  integer              :: varid_MSCCF    ! satellite channel centre frequency
  integer              :: varid_MDPIN    ! data present indicator for flag_table
  integer              :: varid_MSGAPP   ! generating application
  integer              :: varid_MGAPP    ! generating application
  integer              :: varid_MCC      ! cloud type (sequence in new Bufr template)
  integer              :: varid_MPCCO    ! per cent confidence
  integer              :: varid_MMAQC    ! manual/automatic quality control
  integer              :: varid_MNOCT    ! nominal confidence threshold

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
  integer              :: ncvars         ! NetCDF number of variables  defined in this NetCDF file
! integer              :: ncatts         ! NetCDF number of attributes defined in this NetCDF file
! integer              :: unlim_dimid    ! NetCDF id for of unlimited dimension defined in this NetCDF file
  integer              :: numDims        ! NetCDF number of dimensions for individual variable in NetCDF variable
  integer              :: numAtts        ! NetCDF number of attributes for individual variable in NetCDF variable
  integer              :: dimid_varqi    ! NetCDF dimension id for  variables with quality index
  integer              :: dimid_report   !                     for  Reports
  integer              :: dimid_vars     !                     for  Variables

  integer              :: len_varqi      ! number of variables with quality index
  integer              :: len_report     ! number of reports in NetCDF file
  integer              :: len_vars       ! number of variables in data present indicator field in NetCDF file
  integer              :: len_t          ! number of levels in temperature mtn field  in NetCDF file
! integer              :: lev1           ! number of used vertical levels
  integer              :: j, i           ! loop index
  integer              :: nc0            !        dimension for netcdf getvar 1-dimensional arrays
  integer              :: nc1            ! first  dimension for netcdf getvar mpcco in start / count
  integer              :: nc2            ! second dimension for netcdf getvar mpcco in start / count
  integer              :: nc3            ! first  dimension for netcdf getvar mdpin in start / count
  integer              :: nc4            ! second dimension for netcdf getvar mdpin in start / count
  integer              :: nc5            ! first  dimension for netcdf getvar mtn   in start / count
  integer              :: nc6            ! second dimension for netcdf getvar mtn   in start / count
  integer              :: entry1,entry   ! position in source file (subset)


  integer ,allocatable :: ifield1 (:)    !
  real    ,allocatable :: rfield1 (:)    !
  integer ,allocatable :: ifield  (:,:)  !
! real    ,allocatable :: rfield  (:,:)  !
  real    ,allocatable :: rfield2 (:,:)  !
  integer ,allocatable :: ifieldq (:,:)  !

  integer ,allocatable :: qi_cc_aq_t(:,:,:)! quality index array for
                                         ! 1.index : press.,dir.,force,temp.
                                         ! 2.index : qi information
                                         !           1. per cent confidence
                                         !           2. manual/automatic quality control
                                         !           3. nominal confidence threshold
                                         ! 3.index : reports
  integer ,allocatable :: qi_cc(:)       ! over all variable minimum confidence
  integer ,allocatable :: qi_nc(:)       ! over all variable nominal confidence level

  character(NF90_MAX_NAME)   :: yname_v  ! NetCDF variable name


  !----------------
  ! index variables
  !----------------
  integer             :: ir         ! record nr. (subset or loopcount)
! integer             :: is         ! sub-set    index
! integer             :: ie         ! entry in sub-set
! integer             :: id         ! descriptor index

  type (t_amv)        :: s          ! AMV observation read
  type (t_spot)       :: spt0, spti ! observation meta data
  type (t_spot) ,save :: empty      !
  type (t_use)        :: use        ! status variable

  logical             :: lpr_amv    ! amv   reports from netcdf are printed
  logical             :: lpr_extd   ! extended  printing of amv
  integer             :: npr_extd   ! number of extended  printing of aireps
! integer             :: unused     ! no. mnemonics not yet handled
  integer             :: jpos       ! order of defined variables with qi
  integer,dimension(4):: ord_pdft   ! index for press.,dir.,force,temp. in defined variables with qi
  integer             :: minm       ! minimum confidence
  integer             :: min_rep    ! minimum number rot processed reports
  integer             :: jrep       ! loop index over reports
  integer             :: iapp       ! selected generating application

  integer             :: dimids (dimids_max)

  ! logical for meteorological qc variables(over all reports)
  logical             :: l_press, l_wind,  l_mcoct, l_mtn   , l_msidp, l_msdwcm
  logical             :: l_mlase, l_mogc,  l_msccf, l_msgapp, l_mpcco, l_mmaqc
  logical             :: l_mheam, l_mdpin, l_mcc,   l_mmiogc, l_mszan, l_mnoct

  ! logical for reading in parts
  logical             :: l_part
  integer             :: i_met       ! index for begin of meteorological information
  integer             :: i_met_press ! index for pressure
  integer             :: i_met_ndndn ! index for wind direction
  integer             :: i_met_nfnfn ! index for wind force
  integer             :: i_met_mcoct ! index for temperature(coldest)
  integer,dimension(4):: i_m_pdft    ! index array for pres,dir.,for.,temp
  integer             :: sum_qi      ! sum of defined variables with qi
  integer             :: sum_qi_l    ! local sum of defined variables with qi
  integer             :: max_var_nmb ! maximum number of defined variables with qi
  logical             :: newbufr     ! new AMV BUFR template?

!------------------------------------------------------------------------------
  lpr_amv  = .false.; if (netcdf_verb > 1) lpr_amv = .true.
  lpr_extd = .true.
  lpr_extd = .false.
  npr_extd =   2
  if( lpr_extd) npr_extd  = 100
! default
  l_part = .false.
! l_part = .true.

  if( l_part ) then
     min_rep = 10
     min_rep = 100
     min_rep = 300
     min_rep = 600
  endif
  !------------------------------
  ! get default number of reports
  !------------------------------
  len_report =  size (s1date)
  nc0        = len_report
  if( l_part ) then
     nc0     = min(min_rep, len_report)
  endif

  lkeep = .false.
  nkeep = 0
  !---------------------------
  ! get variables ids and name
  !---------------------------
  status = nf90_Inquire (ncid, nVariables=ncvars)
  i_met = 0
  i_met_press = 0
  i_met_ndndn = 0
  i_met_nfnfn = 0
  i_met_mcoct = 0
  jpos     = 0
  ord_pdft = 0
  do j = 1 , ncvars
  status = nf90_Inquire_Variable(ncid, j,name=yname_v )
  if ( i_met == 0 ) then
    if (yname_v(1:7) == 'edition' ) then
      if ( lpr_amv .and. lpr_extd )                                      &
      write (6,'(a,i3,a,a16,a)') '1 nf90_Inquire_Variable(i) ',j,        &
                                    ' Variable name(o): ',trim(yname_v)
    else if ( yname_v(1:7) == 'section' ) then
      if ( lpr_amv .and. lpr_extd )                                      &
      write (6,'(a,i3,a,a16,a)') '2 nf90_Inquire_Variable(i) ',j,        &
                                    ' Variable name(o): ',trim(yname_v)
    else
      i_met = j
      if ( lpr_amv .and. lpr_extd )                                      &
      write (6,'(a,i3,a,a16,a)') '3 nf90_Inquire_Variable(i) ',j,        &
                                    ' Variable name(o): ',trim(yname_v)
    endif
  endif
  if (i_met /= 0 ) then
      if ( lpr_amv .and. lpr_extd )                                      &
      write (6,'(a,i3,a,a16,a)') '4 nf90_Inquire_Variable(i) ',j,        &
                                    ' Variable name(o): ',trim(yname_v)
   endif

  if (i_met /= 0  ) then
    if      ( (yname_v(1:3) == 'MPN'  .or. yname_v(1:3) == 'MPP')  .and. &
              i_met_press  == 0 ) then
      i_met_press = j
      jpos = jpos + 1
      ord_pdft(1) = jpos
    else if ( yname_v(1:5) == 'NDNDN'.and. i_met_ndndn == 0 ) then
      i_met_ndndn = j
      jpos = jpos + 1
      ord_pdft(2) = jpos
    else if ( yname_v(1:5) == 'NFNFN'.and. i_met_nfnfn == 0 ) then
      i_met_nfnfn = j
      jpos = jpos + 1
      ord_pdft(3) = jpos
    else if ( yname_v(1:5) == 'MCOCT'.and. i_met_mcoct == 0 ) then
      i_met_mcoct = j
      jpos = jpos + 1
      ord_pdft(4) = jpos
    endif

    if ( lpr_amv .and. lpr_extd )                                             &
    write (6,'(4(a,i3))') 'i_met',i_met,' i_met_press ',i_met_press,' i_met_ndndn ',i_met_ndndn, &
                          ' i_met_nfnfn ',i_met_nfnfn,' i_met_mcoct ',i_met_mcoct
  endif
  enddo
  i_m_pdft(1) = i_met_press
  i_m_pdft(2) = i_met_ndndn
  i_m_pdft(3) = i_met_nfnfn
  i_m_pdft(4) = i_met_mcoct

  !------------------------
  ! get dimension of fields
  !------------------------
  status = nf90_inq_varid (ncid, 'MPCCO' ,  varid_MPCCO)
  status = nf90_Inquire_Variable(ncid, varid_MPCCO, ndims=numDims, dimids=dimids, natts=numAtts)

  nc1 = 0
  nc2 = 0
  len_varqi = 0
  if(numDims >= 2) then
    dimid_varqi  = dimids(1)
    dimid_report = dimids(2)
    status = nf90_Inquire_Dimension(ncid, dimid_varqi,  len=len_varqi)
    status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
    nc1       = len_varqi
    nc2       = len_report
    if( l_part ) then
      nc1     = len_varqi
      nc2     = min(min_rep, len_report)
    endif
  else
    nc1       = len_report
    nc2       = 1
    len_varqi = max (1,len_varqi)
    if( l_part ) then
      nc1     = min(min_rep, len_report)
    endif
  endif

!----------------------------------------------------------------------

  nc3 = 0
  nc4 = 0
  len_vars = 0
  status = nf90_inq_varid (ncid, 'MDPIN' ,  varid_MDPIN)
  status = nf90_Inquire_Variable(ncid, varid_MDPIN, ndims=numDims, dimids=dimids, natts=numAtts)
! if(status /= nf90_noerr) then
!     print *,'read_amv_netcdf netcdf-error', TRIM(nf90_strerror(status)), &
!             ' no DATA PRESENT INDICATOR defined *** ABORT NOW ***'
!     call finish('read_obs_netcdf ','netcdf error ')
! endif
!
! Note: the new AMV BUFR template does not define MDPIN.  We rely on MCC, see further below.

  if(status == nf90_noerr) then
  if(numDims >= 2) then
    dimid_vars   = dimids(1)
    dimid_report = dimids(2)
    status = nf90_Inquire_Dimension(ncid, dimid_vars ,  len=len_vars )
    status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
    nc3       = len_vars
    nc4       = len_report
    if( l_part ) then
      nc3       = len_vars
      nc4     = min(min_rep, len_report)
    endif
  else
    nc3       = len_report
    nc4       = 1
    len_vars  = max (1,len_vars)
    if( l_part ) then
      nc3     = min(min_rep, len_report)
    endif
  endif
  endif

!----------------------------------------------------------------------

  status = nf90_inq_varid (ncid, 'MTN' ,  varid_MTN)

  if(status == nf90_noerr) then
    status = nf90_Inquire_Variable(ncid, varid_MTN, ndims=numDims, dimids=dimids, natts=numAtts)
  endif

  if(status /= nf90_noerr) then
!     print *,'read_amv_netcdf netcdf-error', TRIM(nf90_strerror(status)), &
!             'no MTN defined *** ABORT NOW ***'
      write(6,*) 'read_amv_netcdf netcdf-error', TRIM(nf90_strerror(status)), &
              'no MTN defined *** RETURN NOW ***'
      write(0,*) 'read_amv_netcdf netcdf-error', TRIM(nf90_strerror(status)), &
              'no MTN defined *** RETURN NOW ***'
      return

!     call finish('read_obs_netcdf ','netcdf error ')


  endif

  nc5 = 0
  nc6 = 0
  len_t = 0
  if(numDims >= 2) then
    dimid_vars   = dimids(1)
    dimid_report = dimids(2)
    status = nf90_Inquire_Dimension(ncid, dimid_vars ,  len=len_t     )
    status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
    nc5       = len_t
    nc6       = len_report
    if( l_part ) then
      nc5       = len_t
      nc6     = min(min_rep, len_report)
    endif
  else
    nc5       = len_report
    nc6       = 1
    len_t  = max (1,len_t)
    if( l_part ) then
      nc5     = min(min_rep, len_report)
    endif
  endif
  !----------------------------------------
  ! define number of reports in source-file
  !----------------------------------------
  i_source = len_report

  !---------------------------------------------
  ! allocate fields dependant from get_var count
  !---------------------------------------------
  allocate (mpn     (nc0))
  allocate (nht     (nc0))
  allocate (ndndn   (nc0))
  allocate (nfnfn   (nc0))
  allocate (mcoct   (nc0))
  allocate (msidp   (nc0))
  allocate (msdwcm  (nc0))
  allocate (mszan   (nc0))
  allocate (mlase   (nc0))
  allocate (mheam   (nc0))
  allocate (mogc    (nc0))
  allocate (msccf   (nc0))
  allocate (mcc     (nc0))

  allocate (mdpin   (len_vars ,nc0))
  allocate (mpcco   (len_varqi,nc0))
  allocate (mmaqc   (len_varqi,nc0))
  allocate (mnoct   (len_varqi,nc0))
  allocate (msgapp  (len_varqi,nc0))

  allocate (ifield  (len_varqi,nc0))
  allocate (ifield1 (          nc0))
  allocate (ifieldq (len_vars ,nc0))
! allocate (rfield  (len_varqi,nc0))
  allocate (rfield1 (          nc0))
  allocate (rfield2 (len_t    ,nc0))

  allocate (qi_cc_aq_t(4,3,nc0))
  allocate (qi_cc     (    nc0))
  allocate (qi_nc     (    nc0))


  !-------------------------------
  ! get meteorological information
  !-------------------------------

  !---------
  ! pressure
  !---------
  mpn = rvind
  l_press = .FALSE.
  status = nf90_inq_varid (ncid, 'MPN'  ,  varid_MPN )
  if (status == nf90_noerr) then
    l_press = .TRUE.

    status = nf90_get_var (ncid, varid_MPN, rfield1, start=(/1/), count=(/nc0/) )

!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_press = .FALSE.
    else
       mpn = rfield1
       where ( rfield1 == rmissing ) mpn = rvind
    endif

  else
     status = nf90_inq_varid (ncid, 'MPPP',  varid_mpn )
     if (status == nf90_noerr) then
       l_press   = .TRUE.
       status = nf90_get_var (ncid, varid_mpn, rfield1, start=(/1/), count=(/nc0/) )
   !   check for missing values
       if ( all( rfield1 == rmissing ))  then
          l_press   = .FALSE.
       else
          mpn   = rfield1
          where ( rfield1 == rmissing ) mpn   = rvind
     endif
    endif

  endif

  !-------
  ! height
  !-------
  nht = rvind
  status = nf90_inq_varid (ncid, 'NHT'  ,  varid_NHT )
  if (status == nf90_noerr) then
    status = nf90_get_var (ncid, varid_NHT, rfield1, start=(/1/), count=(/nc0/) )
    !-------------------------
    ! check for missing values
    !-------------------------
    if ( all( rfield1 == rmissing ))  then
    else
      nht   = rfield1
      where ( rfield1 == rmissing ) nht   = rvind
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
    status = nf90_get_var (ncid, varid_ndndn, ifield1, start=(/1/), count=(/nc0/) )
    if ( all( ifield1 == imissing ))  then
       l_wind = .FALSE.
    else
       where ( ifield1 == imissing ) ifield1 = iud
       rfield1 = ifield1
       where ( rfield1 == rud ) rfield1 = rvind
       ndndn = rfield1
    endif
    if ( l_wind ) then
      status = nf90_inq_varid (ncid, 'NFNFN',  varid_nfnfn)
      if (status == nf90_noerr) then
      status = nf90_get_var (ncid, varid_nfnfn, rfield1, start=(/1/), count=(/nc0/) )
        if ( all( rfield1 == rmissing ))  then
           l_wind = .FALSE.
        else
           where ( rfield1 == rmissing ) rfield1 = rvind
           nfnfn = rfield1
        endif
      endif
    endif
  endif

  !------------
  ! temperature
  !------------
  mcoct = rvind
  l_mcoct = .FALSE.
  l_mtn   = .FALSE.

  status = nf90_inq_varid (ncid, 'MCOCT'  ,  varid_MCOCT )
  if (status == nf90_noerr) then
    l_mcoct = .TRUE.
    status = nf90_get_var (ncid, varid_MCOCT , rfield1, start=(/1/), count=(/nc0/) )
!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_mcoct = .FALSE.
    else
       mcoct = rfield1
       where ( rfield1 == rmissing ) mcoct = rvind
    endif

  else
     status = nf90_inq_varid (ncid, 'MTN' , varid_MTN     )
     if (status == nf90_noerr) then
       l_mtn   = .TRUE.
       status = nf90_get_var (ncid, varid_MTN, rfield2, start=(/1, 1/), count=(/nc5, nc6/) )
   !   check for missing values
       if ( all( rfield2 == rmissing ))  then
          l_mtn   = .FALSE.
       else
          mcoct = rfield2(1,:)
          where ( rfield2(1,:) == rmissing ) mcoct = rvind
      endif
    endif

  endif

  !---------------------------------------------
  ! satellite instrument used in data processing
  !---------------------------------------------
  l_msidp = .FALSE.
  msidp   = -1
  status = nf90_inq_varid (ncid, 'MSIDP'  ,  varid_msidp )
  if (status == nf90_noerr) then
    l_msidp  = .TRUE.
    status = nf90_get_var (ncid, varid_MSIDP, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_msidp = .FALSE.
    else
      msidp = ifield1
      where ( ifield1 == imissing ) msidp = -1
      where ( ifield1 <  0        ) msidp = -1 ! enforce missing value
      where ( ifield1 >  15       ) msidp = -1
    endif
  endif

  !------------------------------------------
  ! satellite derived wind computation method
  !------------------------------------------
  l_msdwcm = .FALSE.
  msdwcm   = -1
  status = nf90_inq_varid (ncid, 'MSDWCM'  ,  varid_MSDWCM )
  if (status == nf90_noerr) then
    l_msdwcm  = .TRUE.
    status = nf90_get_var (ncid, varid_MSDWCM, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_msdwcm = .FALSE.
    else
      msdwcm = ifield1
      where ( ifield1 == imissing ) msdwcm = -1
    endif
  endif

  !-----------------------
  ! satellite zenith angle
  !-----------------------
  l_mszan = .FALSE.
  mszan   = rvind
  status = nf90_inq_varid (ncid, 'MSZAN'  ,  varid_MSZAN )
  if (status == nf90_noerr) then
    l_mszan  = .TRUE.
    status = nf90_get_var (ncid, varid_MSZAN, rfield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( rfield1 == rmissing ))  then
      l_mszan = .FALSE.
    else
      mszan = rfield1
      where ( rfield1 == rmissing ) mszan = rvind
    endif
  endif

  !-------------------
  ! Land/sea qualifier
  ! WMO table 0 08 012
  ! 0 - Land
  ! 1 - Sea
  ! 2 - Coast
  ! 3 - Missing value
  !-------------------
  l_mlase = .FALSE.
  mlase   = -1
  status = nf90_inq_varid (ncid, 'MLASE'  ,  varid_MLASE )
  if (status == nf90_noerr) then
    l_mlase  = .TRUE.
    status = nf90_get_var (ncid, varid_MLASE, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_mlase = .FALSE.
    else
      mlase = ifield1
      where ( ifield1 == imissing ) mlase = -1
    endif
  endif

  !-------------------------
  ! height assignment method
  !-------------------------
  l_mheam = .FALSE.
  mheam   = -1
  status = nf90_inq_varid (ncid, 'MHEAM'  ,  varid_MHEAM )
  if (status == nf90_noerr) then
    l_mheam  = .TRUE.
    status = nf90_get_var (ncid, varid_MHEAM, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_mheam = .FALSE.
    else
      mheam = ifield1
      where ( ifield1 == imissing ) mheam = -1
    endif
  endif

  !------------------------------
  ! originating/generating centre
  !------------------------------
  l_mogc   = .FALSE.
  l_mmiogc = .FALSE.
  mogc     = -1
  status   = nf90_inq_varid (ncid, 'MOGC'  ,  varid_MOGC )
  if (status == nf90_noerr) then
    l_mogc  = .TRUE.
    status = nf90_get_var (ncid, varid_MOGC, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_mogc = .FALSE.
    else
      mogc = ifield1
      where ( ifield1 == imissing ) mogc = -1
    endif
  else
    status = nf90_inq_varid (ncid, 'MMIOGC'  ,  varid_MMIOGC )
    if (status == nf90_noerr) then
    l_mmiogc  = .TRUE.
    status = nf90_get_var (ncid, varid_MMIOGC, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_mmiogc = .FALSE.
    else
      mogc = ifield1
      where ( ifield1 == imissing ) mogc = -1
    endif
    endif
  endif

  !-----------------------------------
  ! satellite channel centre frequency
  !-----------------------------------
  l_msccf = .FALSE.
  msccf   = rvind
  status = nf90_inq_varid (ncid, 'MSCCF'  ,  varid_MSCCF )
  if (status == nf90_noerr) then
    l_msccf  = .TRUE.
    status = nf90_get_var (ncid, varid_MSCCF, rfield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( rfield1 == rmissing ))  then
      l_msccf = .FALSE.
    else
      msccf = rfield1
      where ( rfield1 == rmissing ) msccf = rvind
    endif
  endif

  !--------------------
  ! per cent confidence
  !--------------------
  mpcco = -1
  l_mpcco = .FALSE.
!
! EUMETSAT AQC Scheme(Holmlund(1998) corresponds to second flag
! generating application     method
!     1                      standard method:  comparison  with forecast /with neighbours speed direction height
!     2                      special  method:  comparison  NO   forecast /with neighbours speed direction height
!     3                      NOAA/NESDIS recursive filter function(RFF) (METEOSAT 9 ; not 7)
!
! GOES     AQC Scheme
! generating application     method
!     1                      EUMETSAT QI without using forecast
!     2                      NOAA/NESDIS recursive filter function(RFF)
!     3                      EUMETSAT QI with    using forecast
!
!  consistent choice for QI: EUMETSAT    AQC generating application 2
!                          : NOAA/NESDIS AQC generating application 1
!
! as the japanese and chinese AQC are not available yet
!
! suppose japanese as NOAA/NESDIS (101 i.e 1; 102 2 ; 103 3)
!         chinese  as EUMETSAT
!

! per cent confidence usage flag
  if (pcc_conv == 1) then
    select case (mogc(1))
    case (WMO0_EUMET, WMO0_JMA, WMO0_RSMC)
      status = nf90_inq_varid (ncid, 'MPCCO'  ,  varid_MPCCO )
    case (WMO0_NOAA, WMO0_NASA, WMO0_CIMSS)
      status = nf90_inq_varid (ncid, 'MPCCO1' ,  varid_MPCCO )
    end select
  elseif (pcc_conv == 2) then
    select case (mogc(1))
    case (WMO0_NOAA, WMO0_NASA, WMO0_CIMSS, WMO0_JMA)
      status = nf90_inq_varid (ncid, 'MPCCO'  ,  varid_MPCCO )
    case (WMO0_EUMET, WMO0_RSMC)
      status = nf90_inq_varid (ncid, 'MPCCO0' ,  varid_MPCCO )
    end select
  endif
! special case NPP, NOAA-20, GOES 16, ..., GOES-19 satellite processed by NOAA
  if (mogc(1) == WMO0_NOAA .and. any (istidn(1) == [224, 225, 270, 271, 272, 273])) then
      status = nf90_inq_varid (ncid, 'MPCCO'  ,  varid_MPCCO )
  endif

  if (status == nf90_noerr) then
    l_mpcco = .TRUE.
    status = nf90_get_var (ncid, varid_MPCCO, ifield(:, :), start=(/1,1/), count=(/nc1, nc2/) )

    !   check for missing values
    if ( all( ifield == imissing ))  then
      l_mpcco = .FALSE.
    else
      mpcco = ifield
      where ( ifield == imissing ) mpcco = -1
    endif
  endif

  !---------------------------------
  ! manual/automatic quality control
  !---------------------------------
  mmaqc = -1
  l_mmaqc = .FALSE.
  select case (mogc(1))
  case (WMO0_EUMET, WMO0_RSMC)
    status = nf90_inq_varid (ncid, 'MMAQC0'  ,  varid_MMAQC )
  case (WMO0_NOAA, WMO0_NASA, WMO0_CIMSS, WMO0_JMA)
    status = nf90_inq_varid (ncid, 'MMAQC1'  ,  varid_MMAQC )
  end select
  status = nf90_inq_varid (ncid, 'MMAQC'  ,  varid_MMAQC )
  if (status == nf90_noerr) then
    l_mmaqc = .TRUE.
    status = nf90_get_var (ncid, varid_MMAQC, ifield(:, :), start=(/1,1/), count=(/nc1, nc2/) )

    !   check for missing values
    if ( all( ifield == imissing ))  then
      l_mmaqc = .FALSE.
    else
      mmaqc = ifield
      where ( ifield == imissing ) mmaqc = -1
    endif
  endif

  !-----------------------------
  ! nominal confidence threshold
  !-----------------------------
  mnoct = -1
  l_mnoct = .FALSE.
  select case (mogc(1))
  case (WMO0_EUMET, WMO0_RSMC)
    status = nf90_inq_varid (ncid, 'MNOCT0'  ,  varid_MNOCT )
  case (WMO0_NOAA, WMO0_NASA, WMO0_CIMSS, WMO0_JMA)
    status = nf90_inq_varid (ncid, 'MNOCT1'  ,  varid_MNOCT )
  end select
  status = nf90_inq_varid (ncid, 'MNOCT'  ,  varid_MNOCT )
  if (status == nf90_noerr) then
    l_mnoct = .TRUE.
    status = nf90_get_var (ncid, varid_MNOCT, ifield(:, :), start=(/1,1/), count=(/nc1, nc2/) )

    !   check for missing values
    if ( all( ifield == imissing ))  then
      l_mnoct = .FALSE.
    else
      mnoct = ifield
      where ( ifield == imissing ) mnoct = -1
    endif
  endif

  !--------------------------------------
  ! data present indicator for flag_table
  !--------------------------------------
  l_mdpin = .FALSE.
  mdpin   = -1
  status = nf90_inq_varid (ncid, 'MDPIN'  ,  varid_MDPIN )
  if (status == nf90_noerr) then
    l_mdpin  = .TRUE.
    status = nf90_get_var (ncid, varid_MDPIN, ifieldq, start=(/1, 1/), count=(/nc3, nc4/) )
    !   check for missing values
    if ( all( ifieldq /= 0))  then
!    if ( all( ifieldq == imissing))  then
      l_mdpin = .FALSE.
    else
      mdpin = ifieldq
      where ( ifieldq /= 0 ) mdpin = -1
!      where ( ifieldq == imissing ) mdpin = -1
    endif

  endif

  !--------------------------------------
  ! data present indicator for flag_table
  !--------------------------------------
  l_msgapp = .FALSE.
  msgapp   = -1
  status = nf90_inq_varid (ncid, 'MSGAPP' ,  varid_MSGAPP )
  if (status == nf90_noerr) then
    l_msgapp  = .TRUE.
    status = nf90_get_var (ncid, varid_MSGAPP, ifield, start=(/1, 1/), count=(/nc1, nc2/) )
    !   check for missing values
    if ( all( ifield == imissing ))  then
      l_msgapp = .FALSE.
    else
      msgapp = ifield
      where ( ifield == imissing ) msgapp = -1
    endif
  endif

  if (.not. l_msgapp) then
   status = nf90_inq_varid (ncid, 'MGAPP' ,  varid_MGAPP )
   if (status == nf90_noerr) then
    l_msgapp  = .TRUE.
    status = nf90_get_var (ncid, varid_MGAPP, ifield, start=(/1, 1/), count=(/nc1, nc2/) )
    !   check for missing values
    if ( all( ifield == imissing ))  then
      l_msgapp = .FALSE.
    else
      msgapp = ifield
      where ( ifield == imissing ) msgapp = -1
    endif
   endif
  endif

  !------------------------------------------------------------------
  ! The new AMV BUFR template 3 10 077 defines additional variables.
  ! We detect the new template via presence of one of these, e.g. MCC
  !------------------------------------------------------------------
  newbufr = .false.
  !--------------------------------------
  ! Cloud type
  !--------------------------------------
  l_mcc = .FALSE.
  mcc   = -2
  status = nf90_inq_varid (ncid, 'MCC'  ,  varid_MCC )
  if (status == nf90_noerr) then
    newbufr = .true.
    l_mcc  = .TRUE.
    status = nf90_get_var (ncid, varid_MCC, ifield1, start=(/1/), count=(/nc0/) )
    !   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_mcc = .FALSE.
      mcc   = -1
    else
      mcc = ifield1
      where ( ifield1 == imissing ) mcc = -1
    endif
  endif

  ! Old AMV BUFR templates:
  if (.not. newbufr) then

  !  get variables index for quality index
  sum_qi = 0
  !  get effective number of variables with qi from main obs.part
  max_var_nmb = maxval (i_m_pdft ) - (i_met - 1)
  do j = 1, len_vars
    if ( mdpin(j ,1)  /= -1 ) then
      if (lpr_amv .and. lpr_extd)                                             &
       write (6,'(i3,a,i3,a)') j,'.variable index for quality index',mdpin(j,1)
      if ( j <= max_var_nmb ) then
        sum_qi = sum_qi  + 1
      endif
    endif
  end do
  !----------------------------------
  ! combine qi information in array
  !----------------------------------
  qi_cc_aq_t( :, 1, : ) = -1
  qi_cc_aq_t( :, 2, : ) = -1
  qi_cc_aq_t( :, 3, : ) = -1

  sum_qi_l = 0
  do j = 1 , len_vars
   if ( mdpin(j ,1)  /= -1 ) then
     sum_qi_l = sum_qi_l + 1
     if ( sum_qi_l > sum_qi ) exit

     do i = 1 , 4
       if ( j == i_m_pdft(i) - (i_met - 1) ) then
!      variable with qi found
         if ( ord_pdft(i) /= 0 ) then
            if ( l_mpcco  ) then
              qi_cc_aq_t( i, 1, : ) = mpcco( sum_qi_l , : )
            endif
            if ( l_mmaqc  ) then
              qi_cc_aq_t( i, 2, : ) = mmaqc( sum_qi_l , : )
            endif
            if ( l_mnoct  ) then
              qi_cc_aq_t( i, 3, : ) = mnoct( sum_qi_l , : )
            endif
         endif
       endif
     end do
   endif
  end do

  do jrep = 1,min( nc0,len_report)
    minm = 10000
    do j    = 1,4
      if ( qi_cc_aq_t(j,1,jrep) > 0 ) then
        minm = min( minm, qi_cc_aq_t(j,1,jrep))
      endif
    end do
    if ( minm /= 10000 ) then
      qi_cc(jrep) = minm
    else
      qi_cc(jrep) = -1
    endif
  end do

  do jrep = 1,min( nc0,len_report)
    minm = 10000
    do j    = 1,4
      if ( qi_cc_aq_t(j,3,jrep) > 0 ) then
        minm = min( minm, qi_cc_aq_t(j,3,jrep))
      endif
    end do
    if ( minm /= 10000 ) then
      qi_nc(jrep) = minm
    else
      qi_nc(jrep) = -1
    endif
  end do

  else
! New AMV BUFR template
!
! BUFR Code Table 0 01 044 - Standard generating application
!   0   | Reserved
!   1   | Full weighted mixture of individual quality tests
!   2   | Weighted mixture of individual tests, but excluding forecast comparison
!   3   | Recursive filter function
!   4   | Common quality index (QI) without forecast
!   5   | QI without forecast
!   6   | QI with forecast
!   7   | Estimated error in m/s converted to a percent confidence
! 8-254 | Reserved
!  255  | Missing value
!
! IF Generatin Center == 160 THEN
!   0  first guess
!   1  QI derived from EUMETSAT without forecast fields consistency check
!   2  QI derived from NEDIS RFF (Recursive Filter Function) method
!   3  QI derived from EUMETSAT with forecast fields consistency check
!   4  QI derived from NESDIS EE (estimated error) methods
!   5  QI derived from EUMETSAT without forecast fields consistency check common to all satellite operators

  do jrep = 1,min( nc0,len_report)
     if (any (istidn(jrep) == [3,4,5,55,56,57,70,811,852])) then
        ! Meteosat, (Single-/Dual-)Metop wind products: use QI with forecast
        iapp = 6
     else if (any (istidn(jrep) == [225])) then
        ! NOAA 20, polar wind products: QI without forecast consistency
        iapp = 1
     else
        ! Other satellites: use QI without forecast
        iapp = 5
     end if
!NEC$ loop_count(4)
     do j = 1, nc1
        if (msgapp(j,jrep) == iapp) then
           qi_cc(jrep) = mpcco(j,jrep)
           exit
        end if
     end do
  end do

  endif

  !----------------------------------
  ! list of defined variables in file
  !----------------------------------
  ! default
!  if (l_mpcco) l_msgapp = .TRUE.

   if (netcdf_verb > 0) then
    write (6,'(a,i3,1x,a,a,a,/,3x,18l8)') 'pe=',dace% pe ,               &
     'l_press l_wind  l_mcoct l_mtn   l_msidp l_msdwcm l_mszan l_mlase ',&
     'l_mogc  l_msccf l_mheam l_mdpin l_msgapp l_mpcco l_mmaqc l_mnoct', &
     'l_mmiogc l_mcc',                                                   &
      l_press,l_wind ,l_mcoct,l_mtn  ,l_msidp, l_msdwcm,l_mszan,l_mlase, &
      l_mogc ,l_msccf,l_mheam,l_mdpin,l_msgapp,l_mpcco ,l_mmaqc,l_mnoct, &
      l_mmiogc,l_mcc
   end if

   !-------------------------------
   ! preset total number of reports
   !-------------------------------
   entry   = sum (source(1:ifile-1)% entries)


   !------------------
   ! loop over reports
   !------------------
! original: ir number of subsets
! netcdf  : ir number of report
!  do ir = 1, len_report
   do ir = 1,min( nc0,len_report)

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
   report_subt  = s2ikz(ir)
   bufr_type    = s1cat(ir)
   bufr_subtype = s1catls(ir)
   centre       = s1cent(ir)

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
   head% modtype     = rept_char(obstype)% mod
   head% buf_type    = bufr_type
   head% buf_subtype = bufr_subtype
   head% codetype    = obsid% codetype
   head% time        = stime(ir)
   head% db_time     = db_time(ir)
   head% idbk        = report_subti
   head% source      = ifile
   head% record      = ir
   head% id          = entry1
   head% center      = s1cent(ir)
   head% subcenter   = s1cents(ir)
   head% satid       = istidn(ir)

   if ( lpr_amv .and. ir < npr_extd ) then
     write (6,'()')
     write (6,'( 8(a16, i6  ,/),   &
               & 2(a16,2i8  ,/),   &
               & 7(a16, i6  ,/) )' )                       &
        'pe='         ,dace% pe,                           &
        'head ir='    ,ir,                                 &
        'obstype='    , head% obstype  ,                   &
        'dbkz='       , head% dbkz     ,                   &
        'modtype='    , head% modtype  ,                   &
        'buf_type='   , head% buf_type ,                   &
        'buf_subtype=', head% buf_subtype,                 &
        'codetype='   , head% codetype ,                   &
        'time='       , head% time     ,                   &
        'db_time='    , head% db_time  ,                   &
        'dbk='        , head% idbk     ,                   &
        'source='     , head% source   ,                   &
        'record='     , head% record   ,                   &
        'id='         , head% id       ,                   &
        'center='     , head% center   ,                   &
        'subcenter='  , head% subcenter,                   &
        'satid='      , head% satid
   endif

   !--------------------------------------------
   ! perform simple generic check on report type
   !--------------------------------------------
   call check_report_0 (use, head, 1)
   if (use% state <= STAT_DISMISS) then
     lkeep = .false.
     cycle
   endif

   !------------------
   ! create new report
   !------------------
   spt0       = empty
   spt0% use  = use
   spt0% hd   = head
   spti       = spt0
!  spti% hd% subset = ir
   !----------------------------------------------------------
   ! construct 'empty' Report data structure
   !----------------------------------------------------------
   call construct_amv (s)
   !------------------------------
   ! process the following entries
   !------------------------------
   spti% corme        =   max ( s1updat(ir), 0)
   spti% col% c% dlat =   mlah(ir)
   spti% col% c% dlon =   mloh(ir)
   spti% actual_time  =   obs_time(ir)
   spti% statid       =   ystidn(ir)
   spti% ident        =   istidn(ir)

!  satellite instrument used in data processing
   spti% sttyp  = msidp(ir)
!  satellite derived wind computation method
   spti% stret  = msdwcm(ir)
!  satellite zenith angle
   spti% stzen  = mszan(ir)
!  land/sea qualifier
   spti% stlsf  = mlase(ir)
!  height assignment method
   spti% stclf  = mheam(ir)
!  originating/generating centre
   spti% hd% center = mogc(ir)

!  overall per cent confidence
   spti% pcc  = qi_cc(ir)
!  overall nominal confidence level
   if (.not.newbufr) spti% ncot   = qi_nc(ir)

  !------------------------------------
  ! use station number for station name
  !------------------------------------
  spti% statid = satname (spti% ident)

  if ( lpr_amv .and. ir < npr_extd ) then
    write (6,'()')
    write (6,'(   a20, i6  ,a, / , &
              &   a20, i10 ,   / , &
              & 2(a20,f10.3,   /), &
              &   a20,2i10 ,   / , &
              &   a20, i10 ,   / , &
              &   a20, a   ,   / , &
              & 2(a20, i10 ,   /), &
              &   a20,f10.3,   / , &
              & 5(a20, i10 ,   /) )' )                   &
          'pe=',dace% pe,' amv   spti ',                 &
          'spti% corme        = ', spti% corme        ,  &
          'spti% col% c% dlat = ', spti% col% c% dlat ,  &
          'spti% col% c% dlon = ', spti% col% c% dlon ,  &
          'spti% actual_time  = ', spti% actual_time  ,  &
          'spti% ident        = ', spti% ident        ,  &
          'spti% statid       = ', spti% statid       ,  &
          'spti% sttyp        = ', spti% sttyp        ,  &
          'spti% stret        = ', spti% stret        ,  &
          'spti% stzen        = ', spti% stzen        ,  &
          'spti% stlsf        = ', spti% stlsf        ,  &
          'spti% stclf        = ', spti% stclf        ,  &
          'spti% hd% center   = ', spti% hd% center   ,  &
          'spti% pcc          = ', spti% pcc          ,  &
          'spti% ncot         = ', spti% ncot
  endif

! pressure (vert.location)
  call set_datum (s% p  ,mpn(ir)  ,spti% corme)

! height of cloud top
  call set_datum (s% z  ,nht(ir)  ,spti% corme)

! wind direction/ speed
  call set_datum (s% dd ,ndndn(ir) ,spti% corme)
  call set_datum (s% ff ,nfnfn(ir) ,spti% corme)

! coldest cluster temperature       [K]
  call set_datum (s% t  ,mcoct(ir) ,spti% corme)

! satellite channel centre frequency
  call set_datum (s% f  ,msccf(ir) ,spti% corme)

! per cent confidence
  if (mcc(ir) /= -2) then
    s% p%  pcc = 0
    s% dd% pcc = qi_cc(ir)
    s% ff% pcc = qi_cc(ir)
    s% t%  pcc = 0
  else
    s% p%  pcc = qi_cc_aq_t(1, 1, ir)
    s% dd% pcc = qi_cc_aq_t(2, 1, ir)
    s% ff% pcc = qi_cc_aq_t(3, 1, ir)
    s% t%  pcc = qi_cc_aq_t(4, 1, ir)
  end if

! if Manual/automatic quality control defined
  if ( l_mmaqc ) then
    do j = 1 , 4
    if (  qi_cc_aq_t(j, 2, ir) == 0  .or. qi_cc_aq_t(j, 2, ir) == 1  .or.  &
          qi_cc_aq_t(j, 2, ir) == 5  ) then
       continue
    else
       if ( j == 1 ) s%  p%  pcc = 0
       if ( j == 2 ) s% dd%  pcc = 0
       if ( j == 3 ) s% ff%  pcc = 0
       if ( j == 4 ) s%  t%  pcc = 0
    endif
    end do
  endif
! if nominal confidence threshold defined
  if ( l_mnoct ) then
    do j = 1 , 4
    if (  qi_cc_aq_t(j, 1, ir) <  qi_cc_aq_t(j, 3, ir) ) then
      if ( j == 1 ) s%  p%  pcc = 0
      if ( j == 2 ) s% dd%  pcc = 0
      if ( j == 3 ) s% ff%  pcc = 0
      if ( j == 4 ) s%  t%  pcc = 0
    endif
    end do
  endif

  !----------------------------------------------------------------
  ! derive pressure level if only height is given
  ! z2p_amv = 0 ! do not derive pressure level
  ! z2p_amv = 1 ! derive pressure level from US standard atmosphere
  ! z2p_amv = 2 ! as 1) later re-calculated from background profile
  !----------------------------------------------------------------
  s%  z% pcc = s%  p% pcc
  if (z2p_amv   >  0     .and. &
      s%  p% o  == rvind .and. &
      s%  z% qc == QC_OK       ) then
    s% p% qc = QC_OK
    s% p% src= SRC_DOK
    s% p% o  = p_h_usstd (real (s% z% o, wp))
  endif

  !------------------------------
  ! printout quality control bits
  !------------------------------
  if ( lpr_amv .and. ir < npr_extd ) then
     write (6,'()')
     write (6,'(   a20, i6  ,a, / ,                       &
               & 5(a20,f12.3 , i5,   /) ,                 &
               &   a20,e12.4            )' )              &
           'pe=',dace% pe,' amv   spti obs  pcc',         &
           's%  p% o  %pcc    = ', s%  p% o, s%  p% pcc,  &
           's%  z% o  %pcc    = ', s%  z% o, s%  z% pcc,  &
           's% dd% o  %pcc    = ', s% dd% o, s% dd% pcc,  &
           's% ff% o  %pcc    = ', s% ff% o, s% ff% pcc,  &
           's%  t% o  %pcc    = ', s%  t% o, s%  t% pcc,  &
           's%  f% o          = ', s%  f% o
  endif

  !-------------------------------------------------------------------------------------------
  ! WMO table (0 33 035)
  ! Manual/automatic quality control
  !
  ! Code
  ! figure
  !  0 Automatic quality control passed and not manually checked
  !  1 Automatic quality control passed and manually checked and passed
  !  2 Automatic quality control passed and manually checked and deleted
  !  3 Automatic quality control failed and manually not checked
  !  4 Automatic quality control failed and manually checked and failed
  !  5 Automatic quality control failed and manually checked and re-inserted
  !  6 Automatic quality control flagged data as questionable and not manually checked
  !  7 Automatic quality control flagged data as questionable and manually checked and failed
  !  8 Manually checked and failed
  !  9-14  Reserved
  !  15 Missing value
  !-------------------------------------------------------------------------------------------
  !----------------
  ! standard checks
  !----------------
  lkeep    = .true.
  spti% ps = s% p% o
  call check_report_1 (spti)
  call check_wind_method (spti, s)

  if (lkeep) then
    !-------------------------
    ! check for double entries
    !-------------------------
!   ....
!
    !--------------------------------------------------------
    ! insert in list if most important information is present
    !--------------------------------------------------------
    if (lkeep) call check_store_amv (s, spti, obs, lkeep)
    if (lkeep) then
      nkeep = nkeep + 1
    endif
  endif
  lkeep = nkeep > 0

  !------------------------
  ! end loop over reports
  !------------------------
  enddo

  !-----------
  ! deallocate
  !-----------
  if (allocated(mpn     )) deallocate (mpn     )
  if (allocated(nht     )) deallocate (nht     )
  if (allocated(ndndn   )) deallocate (ndndn   )
  if (allocated(nfnfn   )) deallocate (nfnfn   )
  if (allocated(mcoct   )) deallocate (mcoct   )
  if (allocated(msidp   )) deallocate (msidp   )
  if (allocated(msdwcm  )) deallocate (msdwcm  )

  if (allocated(mszan   )) deallocate (mszan   )
  if (allocated(mlase   )) deallocate (mlase   )
  if (allocated(mheam   )) deallocate (mheam   )
  if (allocated(mogc    )) deallocate (mogc    )
  if (allocated(msccf   )) deallocate (msccf   )
  if (allocated(mdpin   )) deallocate (mdpin   )
  if (allocated(msgapp  )) deallocate (msgapp  )
  if (allocated(mcc     )) deallocate (mcc     )
  if (allocated(mpcco   )) deallocate (mpcco   )
  if (allocated(mmaqc   )) deallocate (mmaqc   )
  if (allocated(mnoct   )) deallocate (mnoct   )

  if (allocated(ifield  )) deallocate (ifield  )
  if (allocated(ifield1 )) deallocate (ifield1 )
  if (allocated(ifieldq )) deallocate (ifieldq )
! if (allocated(rfield  )) deallocate (rfield  )
  if (allocated(rfield1 )) deallocate (rfield1 )
  if (allocated(rfield2 )) deallocate (rfield2 )

  if (allocated(qi_cc_aq_t )) deallocate (qi_cc_aq_t  )
  if (allocated(qi_cc      )) deallocate (qi_cc       )
  if (allocated(qi_nc      )) deallocate (qi_nc       )

  if ( lpr_amv ) then
     write (6,'()')
     write (6,'(a,i3,a)') 'pe=',dace% pe,' ENDE    MO_AMV xxxxxxxxxxxxxxxxxxxxxxx'
     write (6,'()')
  endif

  end subroutine read_amv_netcdf
!------------------------------------------------------------------------------
  subroutine read_nml_obserr_amv
    !========================================
    ! read AMV specific namelist /OBSERR_AMV/
    !========================================
    character(len=8) :: satids(nsat) ! Satellite ids (mnemonics)
    integer          :: meths (nm)   ! Wind computation method
    real(wp)         :: lat_N        ! Latitude separating TR,NH
    real(wp)         :: lat_S        ! Latitude separating TR,SH
    real(wp)         :: err   (nlev) ! Externally specified obs. errors
    real(wp)         :: scale (3)    ! Regional scaling factor

    namelist /OBSERR_AMV/          satids, meths, lat_N, lat_S, err, scale

    integer                     :: i, ierr
    logical                     :: first
    type(t_obserr_amv), pointer :: p(:)
#if defined(__ibm__)
    integer                     :: ios
#endif
    !-------------
    ! set defaults
    !-------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') ' Reading namelist /OBSERR_AMV/'
      write(6,'()')
    endif
    first = .true.
    do
      satids   = ''
      meths    = -1
      lat_N    =  20._wp
      lat_S    = -20._wp
      scale    =   1._wp
      err      = -HUGE (0._wp)
      !---------------------------------
      ! read namelist, consistency check
      !---------------------------------
      if (dace% lpio) then
        call position_nml ('OBSERR_AMV' ,lrewind=first ,status=ierr)
        first=.false.
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=OBSERR_AMV, iostat=ios)
          if(ios/=0)call finish('read_nml_obserr_amv',&
                                'ERROR in namelist /OBSERR_AMV/')
#else
          read (nnml ,nml=OBSERR_AMV)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !------------------------
      ! broadcast to other PE's
      !------------------------
      call p_bcast (satids,   dace% pio)
      call p_bcast (meths,    dace% pio)
      call p_bcast (lat_N,    dace% pio)
      call p_bcast (lat_S,    dace% pio)
      call p_bcast (err,      dace% pio)
      call p_bcast (scale,    dace% pio)

      do i = 2, nlev
         if (err(i) < 0._wp) err(i) = err(i-1)
      end do
      if (any (err <= 0._wp) .or. any (scale <= 0._wp)) then
         call finish ('read_nml_obserr_amv','invalid error specification')
      end if
      !------------------------------------------------
      ! Extend list of observation error specifications
      !------------------------------------------------
      nerr =  nerr + 1
      p    => amv_obserr
      allocate (amv_obserr(nerr))
      if (nerr > 1) then
         amv_obserr(:nerr-1) = p
         deallocate (p)
      end if

      amv_obserr(nerr)% satids = satid (satids)
      amv_obserr(nerr)% meths  = meths
      amv_obserr(nerr)% lat_N  = lat_N
      amv_obserr(nerr)% lat_S  = lat_S
      amv_obserr(nerr)% err    = err
      amv_obserr(nerr)% scale  = scale

      if (any (amv_obserr(nerr)% satids == 0)) then
         if (dace% lpio) &
            write(0,*) "Unknown Satellite: ", &
                       pack (satids, (amv_obserr(nerr)% satids == 0))
         call finish ('read_nml_obserr_amv','invalid satellite name')
      end if

      if (dace% lpio) then
         write(6,'(A,99(:,1x,A))') ' Satellite:',            &
              satname (pack (amv_obserr(nerr)% satids,      &
                             amv_obserr(nerr)% satids /= -1))
         write(6,'(A,20I10)')      '   meths =', pack (meths, meths /= -1)
         write(6,'(A,20F10.6)')    '     err =', err
         write(6,'(A,3F10.6,A)')   '   scale =', scale, "  (NH/TR/SH)"
         write(6,'(A,F10.3)')      '   lat_N =', lat_N
         write(6,'(A,F10.3)')      '   lat_S =', lat_S
         write(6,'()')
      end if
    end do
  end subroutine read_nml_obserr_amv
!------------------------------------------------------------------------------
  function obserr (satid, retrtype, lat, p) result (err)
    real(wp)                   :: err
    integer        ,intent(in) :: satid      ! Satellite id
    integer        ,intent(in) :: retrtype   ! Retrieval type
    real(wp)       ,intent(in) :: lat        ! Latitude [deg]
    real(wp)       ,intent(in) :: p          ! Height   [Pa]

    integer  :: i, k
    real(wp) :: w1, w2

    err = -HUGE (0._wp)
    if (nerr == 0) return
    do i = nerr, 1, -1
      if (     amv_obserr(i)% satids(1) >  -1 .and. &
          all (amv_obserr(i)% satids    /= satid   )) cycle
      if (     amv_obserr(i)% meths (1) >  -1 .and. &
          all (amv_obserr(i)% meths     /= retrtype)) cycle
      if (p >= pt(1)) then
         err = amv_obserr(i)% err (1)
      else if (p <= pt(nlev)) then
         err = amv_obserr(i)% err (nlev)
      else
         do k = 2, nlev
            if (pt(k) < p) then
               w1  = (p-pt(k))/(pt(k-1)-pt(k))
               w2  = 1._wp - w1
               err = w2 * amv_obserr(i)% err(k)   &
                   + w1 * amv_obserr(i)% err(k-1)
               exit
            endif
         end do
      endif
      !--------------------------------------
      ! Apply latitude-dependent scale factor
      !--------------------------------------
      if      (lat >= amv_obserr(i)% lat_N) then
         k = 1
      else if (lat <= amv_obserr(i)% lat_S) then
         k = 3
      else
         k = 2
      end if
      err = err * amv_obserr(i)% scale(k)
      exit
    end do
  end function obserr
!------------------------------------------------------------------------------
  subroutine read_nml_amv_heightcorr ()
    !============================================
    ! read AMV specific namelist /AMV_HEIGHTCORR/
    !============================================
    character(255) :: hgtcor_files(nsat*3)  ! Height correction files

    namelist /AMV_HEIGHTCORR/ hgtcor_files, verbose, p_const_hcor

    integer                     :: i, ierr, j, k, nl
    integer                     :: unit     ! IO unit
    logical                     :: first
    character(len=256)          :: file     ! Height correction file name
    character(len=255)          :: line
    character(len=16)           :: sat_nam
    character(len=8)            :: channel  ! IR, VIS, WV
    integer                     :: satid    ! Numeric satellite identifier
    integer                     :: method   ! Retrieval method
    real(wp)                    :: lev, bias
    real(wp)                    :: lat_nb   ! Northern bound
    real(wp)                    :: lat_sb   ! Southern bound
    integer                     :: stat     ! Parser status
    logical                     :: lstat    ! Pointer status
    type(t_amv_hgtcor), pointer :: p
#if defined(__ibm__)
    integer                     :: ios
#endif
    integer                     :: nhgtcor  ! Number of entries
    character(len=*), parameter :: methods(3) = [ "IR ","VIS","WV " ]
    !-------------
    ! set defaults
    !-------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') ' Namelist /AMV_HEIGHTCORR/'
      write(6,'()')
    endif
    nhgtcor = 0
    do
      hgtcor_files = ''
      verbose      = 2
      p_const_hcor = 0._wp    ! Constant correction above (hPa)
      !---------------------------------
      ! read namelist, consistency check
      !---------------------------------
      if (dace% lpio) then
        call position_nml ('AMV_HEIGHTCORR',lrewind=.true.,status=ierr)
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=AMV_HEIGHTCORR, iostat=ios)
          if(ios/=0)call finish('read_nml_amv_heightcorr',&
                                'ERROR in namelist /AMV_HEIGHTCORR/')
#else
          read (nnml ,nml=AMV_HEIGHTCORR)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !------------------------
      ! broadcast to other PE's
      !------------------------
      call p_bcast (hgtcor_files, dace% pio)
      call p_bcast (verbose,      dace% pio)
      call p_bcast (p_const_hcor, dace% pio)

      p_const_hcor =  p_const_hcor * 100._wp    ! hPa -> Pa

      k = count (hgtcor_files /= "")
      if (dace% lpio) then
         write(6,'(A,A,99(:,/17x,A))') '  hgtcor_files = ',      &
                                   (trim (hgtcor_files(i)), i=1,k)
         write(6,'(A,i6)'            ) '  verbose      = ', verbose
         write(6,'(A,f9.2," hPa")'   ) '  p_const_hcor = ', p_const_hcor/100
         write(6,'()')
      end if

      do i = 1, size (hgtcor_files)
         if (hgtcor_files(i) == "") cycle
         file = path_file (input, hgtcor_files(i))
         if (dace% lpio) then
            unit = get_unit_number ()
            open (unit, file=file, action="read", iostat=ierr)
         end if
         call p_bcast (ierr, dace% pio)
         if (ierr /= 0) call finish ('read_nml_amv_heightcorr', &
                                     'cannot open '//trim (file))
         !---------------------------------------
         ! Read height bias information from file
         !---------------------------------------
         sat_nam = ""
         channel = ""
         ierr    = -1
         if (dace% lpio) then
            !------------------------
            ! Scan for satellite name
            !------------------------
            do
               read (unit,'(A)',iostat=ierr) line
               if (ierr /= 0) exit
               if (line(1:1) /= "#") cycle
               j = index (line, "Satellite")
               if (j == 0) cycle
               j = index (line, ":")
               if (j == 0) cycle
               sat_nam = adjustl (line(j+1:))
               if (sat_nam /= "") ierr = 0
               exit
            end do
         end if
         call p_bcast (ierr,    dace% pio)
         call p_bcast (sat_nam, dace% pio)
         if (ierr /= 0) call finish ("read_nml_amv_heightcorr",        &
                                     "satname missing in "//trim (file))
         satid = get_satid (sat_nam)
         if (satid <= 0) call finish ("read_nml_amv_heightcorr",          &
                                      "invalid satellite "//trim (sat_nam))
         !------------------------
         ! Scan for latitude range
         ! Scan for channels
         !------------------------
         lat_nb =  90
         lat_sb = -90
         ierr   = -1
         stat   = 0
         first  = .true.
         if (dace% lpio) then
            do
               line = ""
               read (unit,'(A)',iostat=ierr) line
               if (ierr /= 0) exit
               if (line(1:1) /= "#") cycle
               j = index (line, "Latitude")
               if (j > 0) then
                  k = index (line, ":")
                  if      (k > 0 .and. index (line, "Max") > 0) then
                     read(line(k+1:),*) lat_nb
                  else if (k > 0 .and. index (line, "Min") > 0) then
                     read(line(k+1:),*) lat_sb
                  else
                     call finish ("read_nml_amv_heightcorr",   &
                                  "invalid line: "//trim (line))
                  end if
                  cycle
               end if
               j = index (line, "Channel:")
               if (j == 0) cycle
               j = index (line, ":")
               channel = adjustl (line(j+1:))
               if (channel == "" .and. first) exit
               select case (channel)
               case ("IR")
                  method  = 1
               case ("VIS")
                  method  = 2
               case ("WV")
                  method  = 3
               case default
                  method  = -1
               end select
               first = .false.
               stat  = 1
               if (.not. associated (amv_hgtcor)) then
                  allocate (amv_hgtcor)         ! Start new list
                  p => amv_hgtcor
               else
                  allocate (p% next)            ! Append to list
                  p => p% next
               end if
               p% satid  = satid
               p% method = method
               p% lat_nb = lat_nb
               p% lat_sb = lat_sb
               k = 0
               do
                  line = ""
                  read (unit,'(A)',iostat=ierr) line
                  if (ierr /= 0) exit
                  if (line(1:1) == "#" .or. line == "") then
                     if (stat == 1) cycle       ! Scan for bias data
                     if (line(1:1) == "#") then
                        backspace (unit)
                     end if
                     exit
                  end if
                  stat = 2
                  read(line,*,iostat=ierr) lev, bias
                  if (ierr /= 0) exit
                  k = k + 1
                  if (k > nhlev) then
                     call finish ("read_nml_amv_heightcorr",            &
                                  "too many levels for "//trim (sat_nam))
                  end if
                  p% lev (k) = lev  * 100  ! level -> Pa
                  p% bias(k) = bias * 100  ! bias  -> Pa
               end do
               stat = 3
               if (k == 0) then
                  stat = -1
                  exit
               end if
               nhgtcor = nhgtcor + 1
            end do
         end if
         call p_bcast (stat,    dace% pio)
         if (stat <= 0) call finish ("read_nml_amv_heightcorr",              &
                                     "no valid channels for "//trim (sat_nam))
         if (dace% lpio) then
            call return_unit_number (unit)
         end if
      end do
      exit
    end do

    if (verbose > 2 .and. dace% lpio) then
       print *, "nhgtcor =", nhgtcor
    end if
    call p_bcast (nhgtcor, dace% pio)
    !----------------------
    ! Broadcast linked list
    !----------------------
    p => amv_hgtcor
    first = .true.
    do
       lstat = associated (p)
       call p_bcast (lstat, dace% pio)
       if (.not. lstat) exit
       if (.not.dace% lpio) then
          if (first) then
             first = .false.
             allocate (amv_hgtcor)         ! Start new list
             p => amv_hgtcor
          else
             allocate (p% next)            ! Append to list
             p => p% next
          end if
       end if
       call p_bcast (p% satid,  dace% pio)
       call p_bcast (p% method, dace% pio)
       call p_bcast (p% lat_nb, dace% pio)
       call p_bcast (p% lat_sb, dace% pio)
       call p_bcast (p% lev,    dace% pio)
       call p_bcast (p% bias,   dace% pio)
       if (dace% lpio) then
          p => p% next
       end if
    end do
    !-----------------
    ! Print final list
    !-----------------
    p => amv_hgtcor
    do
       if (.not. associated (p)) exit
       if (dace% lpio.and. verbose > 0) then
          write(6,'(A,A)'  )     '    Satellite: ', satname (p% satid)
          if (p% method > 0)   &
          write(6,'(A,I10,A10)') '       method: ', p% method,              &
                                                    trim (methods(p% method))
          write(6,'(A,2f10.3)')  '    lat.range: ', p% lat_sb, p% lat_nb
          if (verbose > 1) then
             nl = count (p% lev /= 0._wp)
             write(6,'(A,2F10.3,99(:,/15x,2F10.3))')                   &
                                 '  level, bias: ',                    &
                                 (p% lev(k)/100, p% bias(k)/100, k=1,nl)
          end if
          write(6,'()')
       end if
       if (p% lat_sb >= p% lat_nb) &
          call finish ("read_nml_amv_heightcorr","invalid: lat_min >= lat_max")
       if (any (p% lev(1:nhlev-1) <= p% lev(2:nhlev) .and. &
                p% lev(1:nhlev-1) /= 0._wp                )) then
          call finish ("read_nml_amv_heightcorr",        &
                       "level values must be descending!")
       end if
       p => p% next
    end do

  contains
    !------------------------------------------------------------
    ! Canonicalize satellite names to convention used in mo_satid
    !------------------------------------------------------------
    integer function get_satid (name)
      character(len=*), intent(in) :: name
      !----------------
      ! Local variables
      !----------------
      character(len(name))         :: uname     ! Upper case version of name
      character(len=16)            :: longname  ! Long name for matching

      uname    = toupper (name)
      longname = uname
      if      (uname(1:8) == "METEOSAT" .and. uname(9:9) /= " ") then
         longname(9:) = " "//adjustl (uname(9:))        ! Insert space here
      else if (uname(1:4) == "GOES"     .and. uname(5:5) /= " ") then
         longname(5:) = " "//adjustl (uname(5:))        ! Insert space here
      else if (uname(1:5) == "MTSAT"    .and. uname(6:6) /= "-") then
         longname(6:) = "-"//adjustl (uname(6:))        ! Insert minus here
      else if (uname(1:8) == "HIMAWARI" .and. uname(9:9) /= "-") then
         longname(9:) = "-"//adjustl (uname(9:))        ! Insert minus here
      end if
      get_satid = satid_longname (longname)

    end function get_satid
  end subroutine read_nml_amv_heightcorr
!------------------------------------------------------------------------------
  subroutine amv_heightcorr (obs)
    type(t_obs), intent(inout) :: obs       ! observation container
    !-------------------------------------------
    ! Apply height correction to satellite winds
    !-------------------------------------------
    integer                    :: i, is, k  ! indices
    integer                    :: satid     ! satellite id
    integer                    :: meth      ! wind computation method
    real(wp)                   :: pobs      ! "observed" pressure
    real(wp)                   :: pref      ! reference  pressure for correction
    real(wp)                   :: pnew      ! corrected  pressure
    real(wp)                   :: bias      ! height bias
    real(wp)                   :: w         ! interpolation weight
    type(t_spot)      ,pointer :: s         ! pointer to report meta data
    type(t_amv_hgtcor),pointer :: p         ! pointer to linked list
    integer                    :: np        ! processed winds
    integer                    :: nc        ! corrected winds
    integer                    :: nr        ! rejected  winds

    if (.not. associated (amv_hgtcor)) return

    np = 0
    nc = 0
    nr = 0

    do is = 1, obs% n_spot
       if (obs% spot(is)% hd% obstype /= OT_SATOB    ) cycle
       if (obs% spot(is)% use% state  <= STAT_DISMISS) cycle
       s => obs% spot(is)
       satid =      s% ident
       meth  = mod (s% stret, 100)      ! Remove frequency part
       pobs  =      s% ps               ! "observed" pressure
       pref  = max (pobs, p_const_hcor) ! reference  pressure for correction
       !print *, "is,meth,p=", is,meth,real(pobs)
       np    = np + 1
       !--------------------
       ! Look up height bias
       !--------------------
       p => amv_hgtcor
       do while (associated (p))
          if ( p% satid  == satid                      .and. &
              (p% method == meth .or. p% method == -1) .and. &
               p% lat_nb >= s% col% c% dlat            .and. &
               p% lat_sb <= s% col% c% dlat                  ) then
             bias = p% bias(1)
             if (pref <= p% lev(nhlev)) then
                bias = p% bias(nhlev)
             else if (pref < p% lev(1)) then
                do k = 2, nhlev
                   if (pref > p% lev(k)) then
                      if (p% lev(k) == 0._wp) pref = p% lev(k-1)
                      w    = (p% lev(k) - pref) / (p% lev(k) - p% lev(k-1))
                      bias = (1._wp-w) * p% bias(k) + w * p% bias(k-1)
                      exit
                   end if
                end do
             end if
             pnew = pobs - bias
             do i = s%o%i + 1, s%o%i + s%o%n, 2
                if (obs% varno(i) == VN_U) then
                   obs% body(i)  % bc   = bias
                   obs% body(i+1)% bc   = bias
                   obs% body(i)  % plev = pnew
                   obs% body(i+1)% plev = pnew
                   if (obs% body(i)% lev_typ == VN_P) then
                     obs% olev(i)       = pnew
                     obs% olev(i+1)     = pnew
                   endif
                end if
                !print *, "  pobs,bias,pnew=",real ([pobs,bias,pnew])
             end do
             !-------------------------
             ! Check observation status
             !-------------------------
             if (pnew > 0._wp) then
                nc = nc + 1
             else
                call decr_rpt_use (s, CHK_HEIGHT, comment="height corr.")
                nr = nr + 1
                if (verbose > 1) then
                   write(*,*) "amv_heightcorr: lat,lon,satid,method =", &
                              real(s% col% c% dlat),                 &
                              real(s% col% c% dlon), satid, meth
                   write(*,*) "amv_heightcorr: pobs,bias,pnew =", &
                              pobs,bias,pnew
                end if
             end if
             exit
          end if
          p => p% next
       end do
    end do

    np = p_sum (np)
    nc = p_sum (nc)
    nr = p_sum (nr)

    if (dace% lpio) then
      write(6,'(a,i9)') '    reports processed = ',np
      write(6,'(a,i9)') '    reports corrected = ',nc
      write(6,'(a,i9)') '    reports dismissed = ',nr
      write(6,'()')
      write(6,'(a)')    repeat('-',79)
    endif

  end subroutine amv_heightcorr
!------------------------------------------------------------------------------
end module mo_amv
