!
!+ Interface 3D-Var/LETKF <-> Zenith Total Delay observation operator
!
MODULE mo_std
!
! Description:
!   Interface 3D-Var/LETKF <-> Zenith Total Delay observation operator
!
! Current Maintainer: DWD, Harald Anlauf, Michael Bender
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_22        2013-02-13 Andreas Rhodin
!  Interface 3D-Var/LETKF <-> Zenith Total Delay observation operator
! V1_26        2013/06/27 Michael Bender
!  Update
! V1_28        2014/02/26 Andreas Rhodin
!  changed interface to new_int
! V1_29        2014/04/02 Andreas Rhodin
!  read and process bufr2netcdf ZTD data
! V1_35        2014-11-07 Harald Anlauf
!  implement ICVTYPE_ICON
! V1_42        2015-06-08 Andreas Rhodin
!  bias correction for ground based GNSS data; temporal interpolation for MEC
! V1_43        2015-08-19 Michael Bender
!  First running version of the GNSS bias correction within the 3dvar
! V1_44        2015-09-30 Andreas Rhodin
!  move type t_bufr_inv to mo_t_obs to simplify dependencies
! V1_45        2015-12-15 Michael Bender
!  adjoint code fixed; remove MSIS initialisation; encapsulate read_std_nml
! V1_46        2016-02-05 Harald Anlauf
!  process_std: add pieces for adjoint
! V1_47        2016-06-06 Harald Anlauf
!  finish implementation of adjoint code
! V1_50        2017-01-09 Andreas Rhodin
!  namelist /gpsgb_select/ for selecting GPSGB observations
! V1_51        2017-02-24 Andreas Rhodin
!  use generalised humidity for GPSGB; option to read netcdf input in parallel
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================
#include "tr15581.incf"
!=============
! Modules used
!=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception,   only: finish        ! abort routine
  use mo_kind,        only: wp,          &! working precision kind parameter
                            sp,          &! single  precision kind parameter
                            dp,          &! double  precision kind parameter
                            i8            ! 8-byte integer
  use mo_dace_string, only: char1, char3  ! convert integer -> character
  use mo_namelist,    only: position_nml,&! position namelist
                            nnml,        &! namelist Fortran unit number
                            POSITIONED    ! ok    code from position_nml
  use mo_mpi_dace,    only: dace,        &! MPI group info
                            p_bcast,     &! generic broadcast routine
                            p_min         ! minimum over PEs
  use mo_run_params,  only: path_file,   &! concatenate: path / file . suffix
                            data,        &! default path constatt input
                            obsinput      ! default path for observations
  USE mo_fortran_units, only:            &
                            get_unit_number,   &! get free Fortran unit
                            return_unit_number  ! release Fortran unit
  use mo_time,        only: t_time,      &! date&time data type
                            time_mjd,    &! derive time from modified Julian date
                            mjd,         &! modified Julian date from t_time
                            imm,         &! derive month from t_time
                            cyyyymmddhhmmss
!                           init_time,        &! initialise time data type
!                           iyyyymmdd,ihhmmss,&! conversion routines
!                           cyyyymmdd,chhmmss,&!
!                           i_time             ! t_time->yyyy,mm,dd,hh,mi,ss
  use mo_physics,     only: r2d,              &! factor radians -> degree
                            d2r,              &! factor degree  -> radians
                            gacc,             &! gravity acceleration
                            pi,               &! 3.1415....
                            rh_q,             &! relative <- specific humid.
                            q_rh,             &! relative -> specific hum.
                            q_rh_adj           ! adjoint routine
  use mo_cntrlvar,    only: gh_rh,            &! generalised humidity from rh
                            rh_gh              ! relative    humidity from gh
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,         only: nf90_Inquire,          &!
                            nf90_Inquire_Dimension,&!
                            nf90_inq_dimid,        &!
                            nf90_Inquire_Variable, &!
                            nf90_inq_varid,        &!
                            nf90_get_var,          &!
!                           nf90_strerror,         &!
!                           NF90_FLOAT,            &!
!                           NF90_INT,              &!
                            NF90_MAX_NAME,         &!
                            NF90_NOERR              !
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only:ncid,           &! NetCDF file id
                          dimids_max,     &! max number of NetCDF dimension ids
!                         imissing,       &! NetCDF _FillValue for integer
                          rmissing,       &! NetCDF _FillValue for reals
                          ymissing,       &! NetCDF _FillValue for character
                          s1cent,         &! data centre
                          stime,          &! header observation time (section1)
                          db_time,        &! data bank time
                          s1cents,        &! data sub centre
                          istidn,         &! WMO numeric station number combined
                          s1updat,        &! update sequence no.
                          mlah,           &! latitude
                          mloh,           &! longitude
                          obs_time,       &! body observation time
                          ystidn,         &! station identifier as variable
                          get_int,        &! read integer variable from NetCDF file
                          get_real         ! read real    variable from NetCDF file
  !------------------------
  ! access matrix data type
  !------------------------
  use mo_dec_matrix,  only: t_vector_segm   ! vector segment data type
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use Errors,         only: error_status,  &! status data type definition
                            enter_callee,  &!
            error_callee => error,         &!
                            display_status,&!
                            clear_status    !
  use MSIS,           only: MSIS_init       ! initialise MSIS atmosphere
  use mo_atm_state,   only: t_atm           ! atm. state data type
  USE mo_atm_transp,  only: scatter_level, &! scatter one horizontal grid level
                            gather_multi    ! gather  multi level field
  use mo_atm_grid,    only: t_grid,        &! atm. grid  data type
!                           VCT_P_HYB,     &! GME/HRM vertical coordinate
!                           VCT_Z_HYB,     &! COSMO/ICON vertical coordinate
                            VCT_Z_GEN       !
  use mo_grid_intpol, only: Grid_Indices,  &! get indices of neighbor gridpoints
                            mp,            &! max. number of interpolation points
                            alloc_imcol,   &! (re)allocate t_imcol
                            add_index       ! utility to set model column indices
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_obs_set,   only:  t_obs_block     ! obs data type
  use mo_t_obs,     only:                 &!
                           t_obs,         &! observation data type (container)
                           t_spot,        &! report data type
                           t_head,        &! component of t_spot
!                          t_imcol,       &! model column meta data type
!                          derive_dbkz,   &! derive DBKZ if not present
                           new_spot,      &! reserve memory
                           new_obs,       &! reserve memory
                           new_int,       &! reserve memory for interp.space
                           new_par,       &! reserve memory
                           set_xuv,       &! set unit vectors, zenith angle
!                          set_vqc_insitu,&! subroutine to set VQC bounds
                           shrink_report, &! release unused obs. in report
                           source,        &! list   of Report source files
!                          corme_conv,    &! corme handling (TEMP/BUFR)
                           TSK_INIT,      &! FLAGS: initialize module
                           TSK_READ,      &!  read observations
                           TSK_SET_CHR,   &!  set observation characteristics
                           TSK_SETUP_COLS,&!  determine model columns required
                           TSK_SETUP_FULL,&!  setup description of PSAS-space
                           TSK_SETUP_FUL0,&!  setup description of PSAS-space
                           TSK_SHRINK,    &!  release unused obs. in report
                           TSK_K,         &!          evaluate linear operator
                           TSK_Y,         &!          evaluate nonlinear oper.
                           TSK_R,         &!   setup observational errors
!                          CHR_ID,        &! H is the identity operator
                           CHR_NONL,      &! H is nonlinear
!                          CHR_INV,       &! H is invertible
                           CHR_EXP,       &! H is 'expensive'
                           OBS_TV,        &! interp. observation type: Tv
                           OBS_RH,        &! interp. observation type: rh
!                          OBS_HS,        &! interp. observation type: geop.
                           OBS_H,         &! interp. observation type: geop.
!                          ITY_ICOL,      &! interpolation type: column
                           source,        &! list   of source files
                           n_source,      &! number of source files
!                          m_source,      &! max. number of source files
                           add_source,    &! add source file to list
                           bufr_inv,      &! observation file inventory
                           p_bcast,       &!
                           netcdf_verb,   &! verbosity of NetCDF decoding
                           GPSGB,         &! module type
                           ITY_MCOLS,     &! interpolation type: column
                           FT_SPECIFIC     ! NetCDF file type flag
  use mo_t_col,      only: t_cols,        &! model columns data type
                           t_col,         &! 1 column of the model
                           COL_T,         &! column flag: temperature
                           COL_Q,         &!              spec.humidity
!                          COL_X,         &!              water + ice load
!                          COL_TV,        &! column flag: virt.temp.
!                          COL_RH,        &!              rel. hum .
!                          COL_GEOH,      &!              geop.h. half levs.
                           COL_GEO,       &!              geop.h. full levs.
                           COL_P           !              press.  full levs.
!                          COL_PH          !              press.  half levs.
  use mo_t_use,      only: STAT_FORGET,   &!
                           STAT_ACTIVE,   &!
                           STAT_REJECTED, &!
                           STAT_DISMISS,  &!
                           STAT_OBS_ONLY, &!
                           CHK_NONE,      &!
                           CHK_INSDAT,    &!
                           CHK_DOMAIN,    &!
                           CHK_HEIGHT,    &!
                           CHK_OPERATOR,  &!
                           CHK_WHITELIST, &!
                           t_use,         &! data type to hold state
                           use_0,         &! default values of type use
                           decr_use,      &! decrease state of datum
                           reverse_code,  &! convert code back (feedback to DACE)
                           stats           ! status translation table
  use mo_t_datum,    only: t_datum         ! date+time derived type
!                          inv_datum       ! default t_datum variable
  use mo_obs_tables, only: rept_use,      &! report type usage table
                           decr_rpt_use,  &! degrade status of report
!                          idb_dbk,       &! index in table rept_stat
                           check_report_1,&! basic checks on reports
                           check_report_0  ! basic checks on reports
  use mo_obs_rules,  only: get_rule,      &! routine to get a rule
                           t_set,         &! result data type
                           iud,           &! undefined integer value
                           rud             ! undefined real value
  use mo_vqc,        only: svqc,          &! default var.qual.cntrl stdev.
                           vqc_form        ! formulation of obs-costfunction
  use mo_fdbk_tables,only: OT_GPSGB,      &! observation type:GNSS ground based
                           OC_GPS,        &! observation code type: ZTD
                           OC_GPSGB,      &! observation code type: STD
                           VN_SPD,        &! observed quantity:slant path delay
                           VN_ZPD,        &! observed quant.: zenith path delay
                           VN_ELEV,       &! level type     : elevation
                           FL_PRACTICE     ! insufficient data flag
  use mo_wmo_tables, only: WMO0_DWD,      &! DWD   Offenbach
                           WMO6_LATLON,   &! Latitude/Longitude Grid
                           WMO6_ROTLL,    &! Rotated latitude/longitude grid
                           WMO6_GAUSSIAN, &! Gaussian grid
                           DWD6_ICON,     &! ICON unstructured grid id
                           DWD6_ICOSAHEDRON! Icosahedral based triangular grid
  use mo_gnss_bc    ,only: biascor_mode,  &! mode used for bias correction
                           t_decay,       &! biasc. accum. decay time (days)
                           n_required,    &! number of entries required for bc
                           bc_fallback     ! action if biasc-file not present
  !--------------------------------------------------------
  ! access generic (DA-sytem or COSMO) observation operator
  !--------------------------------------------------------
  use mo_std_selection, only:                 &!
                          std_status,         &! returns initial status
                          name2spec,          &! derive processing specifications
                          print_tables,       &! print GNSS stations & centers
                          read_nml_gpsgb_select! read namelist /STD_SELECT/
  use mo_std_coord, only: SlantHeader,        &! observation file header type
                          GNSSStation,        &! station data     derived type
                          SlantTotalDelay,    &! observation data derived type
                          k1, k2, k3,         &! constants for Thayer refrac.
                          RDRD,               &! R(dry)/R(vapor)
                          invalid,            &! invalid data
                          invalsp,            &! invalid data
                          Ellips2Cart,        &! calculate cartesian coordinate
                          WGS84Param           ! WGS84 reference ellipse
  use mo_std_operator, only:                  &!
         SlantData,          &! Slant data       derived type
         p_column,           &! model column derived type
         !---------
         ! Routines
         !---------
         read_std_nml,       &! read namelist, module initialisation
         read_gnss_stations, &! read GNSS station information
         scan_std_obs,       &! read observation file headers
         read_std_obs,       &! read the observations
         std_path_coord,     &! compute points on the slant path
         destruct,           &! deallocate components
         std_delay,          &! compute model slant total delay
         !--------------------
         ! namelist parameters
         !--------------------
         read_ascii,         &! flag to read ascii files
         verbose,            &! verbosity level
         NStepVertMod,       &! Number of vertical points inside the model
         NStepVertTop,       &! Number of vertical points above the model
         NStepOpt,           &! option for scaling points
         Hmax,               &! Maximum height for STD integration in m
         HScaleP,            &! pressure scale height, model
         HScaleP2,           &! pressure scale hgt, above model
         Hlevel,             &! model levels per interval
         Heights,            &! heights of model levels at upper interval bound.
         Hpoints,            &! supporting points per interval
         Href,               &! reference height in m
         UseRaytracer,       &! use raytracer for low elev.
         ZTDminUse,          &! minimum ZTD used
         ZTDmaxUse,          &! maximum ZTD used
         StatBelowSurface,   &! reject stations below model surface
         StatAboveSurface,   &! reject stations above model surface
         StatBelowColumn,    &! reject stations below max. column
         ZTDerror,           &! ZTD observation error, m
         HORIFile,           &! STD station file
         std_obs_file,       &! STD observation input files
         ztd_col,            &! /= 0 : use 1 model column for ZTDs
         pl_method,          &! method for pressure level estimation
         !----------------------
         ! direct access to data
         !----------------------
         CoordList,          &! list of stations
         GNSShash,           &! hash table: station id -> index
         STD                  ! observational data

  !------------------------------------
  ! Interpolation routines and adjoints
  !------------------------------------
  use Interpolation,only: Init_Spline,        &! (from GPSRO operator)
                          Spline

  !--------------------------------------
  ! Debugging of implemention of adjoint?
  !#define DEBUG_STD_AD
  !--------------------------------------
#undef  DEBUG_STD_AD
#ifdef  DEBUG_STD_AD
  use Interpolation_adj,  only:               &! (from GPSRO operator)
                          Init_Spline_adj,    &!
                          Spline_adj
#endif

implicit none

!================
! Public entities
!================
private
public :: process_std       ! perform specific tasks for the STD operator
public :: read_std_netcdf   ! read (ZTD) from bufrx2netcdf converted files
public :: read_std_nml_dace ! read namelist /STD_OBS/ and broadcast variables
public :: std_col2xi        ! convert t,q ,gp  to  t,rh,gh
public :: read_fdbk_gpsgb   ! restore SlantData from t_spot, body
public :: pl_method         ! method for pressure level estimation
public :: set_std           ! initialize auxiliary data used by STD operator
!------------------------------------------------------------------------------
  interface p_bcast
    module procedure p_bcast_GNSSStation
  end interface p_bcast
!------------------------------------------------------------------------------

  integer ::  sd_size = 0  ! size of type SlantData
  integer :: int_size = 0  ! size of type integer
  integer ::  wp_size = 0  ! size of type real(wp)

  real(wp) ,parameter :: Eps = 1._wp/RDRD - 1._wp       ! R(vapor)/R(dry) - 1

  ! Character set for GNSS station names
  ! Station names with other characters will be rejected.
  ! Valid GNSS station names are defined by the RINEX standard and
  ! "must contain capital ASCII letters or numbers".
  ! See https://igs.org/wg/rinex
  character (len=38), parameter   ::                           &
            GNSSset = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_'

!==============================================================================
contains
!==============================================================================
  subroutine read_std_nml_dace
  !------------------------------------------------
  ! read namelist /STD_OBS/ and broadcast variables
  !------------------------------------------------
    integer :: ierr ! returned from position_nml

    if (dace% lpio) then
      !-------------------------------
      ! read namelist on I/O processor
      !-------------------------------
      call position_nml ('STD_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
        call read_std_nml (nnml)
      case default
        write (6,*) ' read_std_nml_dace: namelist /STD_OBS/ not present'
        call read_std_nml (nnml, initonly=.true.)
      end select
    endif
    !-------------------
    ! broadcast namelist
    !-------------------
    call p_bcast (ierr,              dace% pio)
    call p_bcast (read_ascii,        dace% pio)
    call p_bcast (NStepVertMod,      dace% pio)
    call p_bcast (NStepVertTop,      dace% pio)
    call p_bcast (NStepOpt,          dace% pio)
    call p_bcast (Hmax,              dace% pio)
    call p_bcast (HScaleP,           dace% pio)
    call p_bcast (HScaleP2,          dace% pio)
    call p_bcast (Hlevel,            dace% pio)
    call p_bcast (Heights,           dace% pio)
    call p_bcast (Hpoints,           dace% pio)
    call p_bcast (Href,              dace% pio)
    call p_bcast (ZTDminUse,         dace% pio)
    call p_bcast (ZTDmaxUse,         dace% pio)
    call p_bcast (StatBelowSurface,  dace% pio)
    call p_bcast (StatAboveSurface,  dace% pio)
    call p_bcast (StatBelowColumn,   dace% pio)
    call p_bcast (k1,                dace% pio)
    call p_bcast (k2,                dace% pio)
    call p_bcast (k3,                dace% pio)
    call p_bcast (UseRaytracer,      dace% pio)
    call p_bcast (ZTDerror,          dace% pio)
    call p_bcast (std_obs_file,      dace% pio)
    call p_bcast (HORIFile,          dace% pio)
    call p_bcast (ztd_col,           dace% pio)
    call p_bcast (pl_method,         dace% pio)

    call p_bcast (biascor_mode,      dace% pio)
    call p_bcast (t_decay,           dace% pio)
    call p_bcast (n_required,        dace% pio)
    call p_bcast (bc_fallback,       dace% pio)
    call p_bcast (verbose,           dace% pio)

  end subroutine read_std_nml_dace
!==============================================================================
  subroutine process_std (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
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
  !-----------------------------------------------------------------------
  ! This subroutine is called from various points in the assimilation code
  ! in order to perform specific tasks for the Zenith Total Delay
  ! operator.
  !------------------------------------------------------------------------

    integer                   :: tsk          ! task, local copy
    integer                   :: if           ! file  index variable
    integer                   :: is           ! slant index variable
    integer                   :: ip           ! point index along slant
    integer                   :: ig           ! point index along neighbours
    type(SlantHeader)         :: head         ! observation file header
    character(len=256)        :: fname        ! file name
    integer                   :: n_rec        ! number of records in file
    integer                   :: n_stat = 0   ! size of list of stations
    type(SlantData)  ,pointer :: slants(:)    ! slant data for one station
    real(wp)                  :: htop = -1._wp! model top geoid corrected
    integer                   :: nnghb= -1    ! # neighbours in interpolation
    integer                   :: nngh3= -1    ! # neighbours in DACE interface
    integer                   :: nngho= -1    ! # columns in operator
    integer                   :: order        ! parameter to grid_indices
    type (t_grid)    ,pointer :: grid         ! atmospheric grid
!   real(wp)         ,pointer :: geoid(:,:,:)
    type(SlantData)  ,pointer :: s            ! pointer to actual slant
    type(p_column),allocatable:: mc(:,:)      ! pointer array to model columns
    type(p_column)   ,pointer :: h (:,:)      ! adjoint
    real(wp)                  :: delay        ! model slant total delay
    logical                   :: adjoint      ! flag for TSK_K set
    logical                   :: change       ! true for observations changed
    logical          ,pointer :: msk(:)       ! true for observations kept
    type(t_set)               :: set          ! set of characteristics for channels
    real(wp)                  :: z            ! model surface elevation

    integer                   :: idx (mp,4)
    integer                   :: np
    real(wp)                  :: w     (mp)
    integer                   :: ixmcol(mp)   ! indices to imcol
    integer                   :: mimcol = 200 ! estimated #of points in slant
    integer                   :: nimcol       ! # of model columns required
    integer(i8)               :: iatm         ! model columns required
    integer                   :: natm         ! # of model columns required
    integer                   :: m            ! number of model columns
    integer                   :: i, j, k, l, n! indices
    integer                   :: ic           ! index
    type(t_col)      ,pointer :: c            ! pointer to model column
    integer                   :: j1           ! index
    type (GNSSStation),target :: station
    logical                   :: tmp_stat
!   integer                   :: month
!   type(Error_Status),pointer:: errstat
    integer                   :: err          ! return code std_path_coord/delay
    integer                   :: colerr       ! return code set_std_data
    integer                   :: Nmiss        ! number of slants set to DISMISS
    integer                   :: ix           ! column index
    integer                   :: ke           ! number of model levels
    integer                   :: nn           ! number of variables / column
    integer                   :: nx, ny
    real(wp)                  :: dx, dy
    real(wp)                  :: dz           ! work
    real(wp)                  :: a1, a2, a3   ! work
    real(wp)     ,allocatable :: rh  (:), q  (:), t  (:), p  (:) ! variables
    real(wp)     ,allocatable :: rh_a(:), q_a(:), t_a(:), p_a(:) ! adjoints
    real(wp)     ,allocatable :: gh_a(:), gh (:)                 ! gen.humidity
    real(wp)     ,allocatable :: lnp (:)      ! ln(pf)
    real(wp)     ,allocatable :: geo (:)      ! geopotential
!   real(wp)     ,allocatable :: tv  (:)      ! aux.variables
    real(wp)     ,allocatable :: Hstd (:)     ! Jacobi matrix (part of)
    real(wp)     ,allocatable :: std_t(:,:)   ! Adjoint of slant total delay
    real(wp)     ,allocatable :: std_p(:,:)   ! Adjoint of slant total delay
    real(wp)     ,allocatable :: std_q(:,:)   ! Adjoint of slant total delay
    real(wp)     ,allocatable :: delay_t(:)   ! Adjoint of std w.r.t. t,q,geo
    real(wp)     ,allocatable :: delay_q(:)   ! Adjoint of std w.r.t. t,q,geo
    real(wp)     ,allocatable :: delay_z(:)   ! Adjoint of std w.r.t. t,q,geo
    real(wp)     ,allocatable :: d2t   (:)    ! spline, 2nd derivative
    real(wp)     ,allocatable :: d2q   (:)    ! spline, 2nd derivative
    real(wp)     ,allocatable :: d2lnp (:)    ! spline, 2nd derivative
    real(wp)     ,allocatable :: t_z   (:)    ! 1st derivative of spline
    real(wp)     ,allocatable :: q_z   (:)    !   w.r.t. support points
    real(wp)     ,allocatable :: p_z   (:)    !
#ifdef  DEBUG_STD_AD
    real(wp)     ,allocatable :: d2t_z (:,:)  ! spline, 2nd derivatives,
    real(wp)     ,allocatable :: d2q_z (:,:)  !   derivatives w.r.t.
    real(wp)     ,allocatable :: d2p_z (:,:)  !   input parameters
    real(wp)     ,allocatable :: d2t_t (:,:)  !
    real(wp)     ,allocatable :: d2q_q (:,:)  !
    real(wp)     ,allocatable :: d2p_p (:,:)  !
    real(wp)                  :: tmp          ! work
    real(wp)     ,allocatable :: tmp1  (:)    ! work
    real(wp)     ,allocatable :: tmp2  (:)    ! work
    real(wp)     ,allocatable :: t_int (:)    ! interpolated profile: t
    real(wp)     ,allocatable :: q_int (:)    ! interpolated profile: q
    real(wp)     ,allocatable :: p_int (:)    ! interpolated profile: ln(p)
    real(wp)     ,allocatable :: t_by_z(:,:), t_by_t(:,:) ! adjoints
    real(wp)     ,allocatable :: q_by_z(:,:), q_by_q(:,:) ! adjoints
    real(wp)     ,allocatable :: p_by_z(:,:), p_by_p(:,:) ! adjoints
#endif

    !==============================
    ! observation non_specific part
    !==============================
    if (rept_use(OT_GPSGB)% use (CHK_NONE) <= STAT_FORGET) return
    tsk = task

    !----------------
    ! tsk == TSK_INIT
    !----------------
    if (iand (TSK_INIT,tsk) /= 0) then
      tsk=tsk-TSK_INIT
      !----------------------
      ! module initialisation
      !----------------------

      !---------------
      ! read namelists
      !---------------
      call read_std_nml_dace      ! /STD_OBS/
      call read_nml_gpsgb_select  ! /GPSGB_SELECT/

!     !----------------
!     ! initialise MSIS
!     !----------------
!     if (ierr==POSITIONED) then
!       month = imm (atm% time)
!       nullify (errstat)
!       call enter_callee ('process_std',errstat)
!       call MSIS_init (month, errstat)
!       if (error (errstat)) then
!         call display_status (errstat)
!         call finish ('process_std','error in MSIS_init')
!       endif
!       call clear_status (errstat)
!     endif

      n_stat = 0
      if (read_ascii) then
        if (dace% lpio) then
          !------------------------------------------------
          ! append default path names to station input file
          ! read station coordinates on I/O processor
          !------------------------------------------------
          HORIFile = path_file (data, HORIFile)
          call read_gnss_stations
          n_stat = size (CoordList)
        endif
        !------------------------------
        ! broadcast station coordinates
        !------------------------------
        call p_bcast (n_stat,       dace% pio)
        if (.not.dace% lpio) allocate (CoordList (n_stat))
        call p_bcast (CoordList,    dace% pio)
        call p_bcast (GNSShash,     dace% pio)
        !----------------------------------------------------
        ! scan observation input files, update file inventory
        !----------------------------------------------------
        do if = 1, size (std_obs_file)
          if (std_obs_file(if) == '') cycle
          fname = path_file (obsinput, std_obs_file(if))
          if (dace% lpio) call scan_std_obs (fname, head)
          call p_bcast (head% NStat,   dace% pio)
          call p_bcast (head% NSlants, dace% pio)
          n_rec = head% NStat
          if (head% NSlants > 0) then
            call add_source (obsinput, std_obs_file(if),&
                             filetype= FT_SPECIFIC     ,&
                              obstype= OT_GPSGB        ,&
                             complete= .true.          ,&
                              entries= n_rec            )
            if (dace% lpio) then
              if(n_source>1) bufr_inv% subseto(n_source) =   &
                bufr_inv% subseto(n_source-1)
              bufr_inv   (OT_GPSGB)% file(n_source)    = .true.
              bufr_inv   (OT_GPSGB)% nrec              =            &
                bufr_inv (OT_GPSGB)% nrec              + n_rec
              bufr_inv   (OT_GPSGB)% nsubset           =            &
                bufr_inv (OT_GPSGB)% nsubset           + n_rec
              bufr_inv   (OT_GPSGB)% subseto(n_source) =     &
                bufr_inv (OT_GPSGB)% subseto(n_source) + n_rec
            endif
          endif
          call p_bcast (bufr_inv, dace% pio)
        end do
      endif
    endif
    if (tsk==0) return


    !-----------------------
    ! tsk == TSK_READ:
    ! read data
    !-----------------------
    if (iand (TSK_READ,tsk) /= 0) then
      if (read_ascii) then
        call read_std_ascii (obs% o)
      endif
      tsk=tsk-TSK_READ
      call print_tables
      if (tsk==0) return
    endif


    !-------------------------------------------------
    ! TSK_SETUP_COLS:
    ! determine model columns required by the operator
    !-------------------------------------------------
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      tsk=tsk-TSK_SETUP_COLS

      !--------------
      ! set model top
      !--------------
      grid => atm% grid
!     geoid => grid% geoid (grid% lb(1):grid% ub(1),  &
!                           grid% lb(2):grid% ub(2),1,&
!                           grid% lb(4):grid% ub(4)   )
      htop = grid% htopf
      if (htop == 0._wp) call finish ('process_std','GPSGB: htopf not defined')

      !----------------------------------------------------
      ! estimate horizontal grid spacing (from mo_veri_obs)
      !----------------------------------------------------
      select case (atm% grid% gridtype)
      case (WMO6_GAUSSIAN)
         dx = atm% grid% d_deg
         dy = atm% grid% d_deg
      case (WMO6_LATLON, WMO6_ROTLL)
         dx = atm% grid% di
         dy = atm% grid% dj
      case (DWD6_ICOSAHEDRON)
         nx = atm% grid% ni
         ny = atm% grid% ni
         dx = 63.486_sp/nx
         dy = 63.486_sp/ny
      case (DWD6_ICON)
         nx = atm% grid% ni
         ny = atm% grid% ni
         dx = 45.416_sp/nx   ! horizontal grid spacing, degrees
         dy = 45.416_sp/ny   !     "
      case default
         dx = -1.0
         dy = -1.0
         write(0,*) "gridtype =", atm% grid% gridtype
         !call finish('add_veri','unsupported gridtype')
      end select

      !--------------------------------
      ! set number of neighbour columns
      !--------------------------------
      select case (grid% gridtype)
      case (DWD6_ICOSAHEDRON)      ! GME
        nnghb = 3
      case (DWD6_ICON)             ! ICON
        nnghb = 3
      case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
        nnghb = 4
      case default
        call finish('process_std',                                     &
                    'GPSGB: gridtype not implemented: '//char3(grid% gridtype))
      end select

      call load_std (obs% o, spot, slants)
!     grid => atm% grid
      iatm = COL_T+COL_Q+COL_P+COL_GEO ! parameters required (T+Q)
      natm = 4                         ! number of parameters required
!     iatm = iatm + COL_TV             ! virtual temperature
!     natm = natm + 1                  ! (needed by adjoint)???
      !-------------------------------------------------------------
      ! set number of gridpoints already associated with this report
      !-------------------------------------------------------------
!     if (associated(imcol)) then
!       nimcol = size(imcol)
!     else
        nimcol = 0
        call alloc_imcol (spot% imcol, mimcol)
!     endif

      !--------------------------------------------------------------
      ! derive (temporary) observation info if not present (BUFR ZTD)
      !--------------------------------------------------------------
      tmp_stat = .not. associated (slants (1)% Station)
      if (tmp_stat) call spot2station (spot, grid, slants, station)

      !-----------
      ! slant loop
      !-----------
      Nmiss = 0
      do is = 1, size (slants)
        s => slants (is)
        if (tmp_stat) slants(is)% Station => station

        !-----------------------------------------
        ! get the coordinates along the slant path
        !-----------------------------------------
        call std_path_coord (s, htop, err, grid%model, grid%nz, dx)

        !------------------------------------------------------------------
        ! nngh3 : number of model columns used in DACE
        ! nnghb : number of model columns required for interpolation
        ! nngho : number of model columns passed to the the operator
        ! for ZTD nngh3 may differ (=1) from nnghb to speed up computations
        ! in this case the nearest neighbour column is passed nnghb times
        !------------------------------------------------------------------
        nngh3   = nnghb
        nngho   = nnghb
        order   = 2
        if (abs (s% obs% elevation * r2d - 90._wp) < 1.e-10_wp .and. ztd_col/=0) then
          nngh3                = 1
          order                = 1
          if (ztd_col>0) then
            nngho   = 1
            s% Nmod = 1
            s% Ntot = 1
          endif
        endif

        allocate (s% wi             (nngho, s% Nmod))
        allocate (s% idxi           (nngho, s% Nmod))
        allocate (s% SlantPointsGeo (       s% Ntot, 3))
        s% Nnghb = nngho
        s% Nhor  = s% Nmod
        s% Naloc = s% Nmod
        s% htop  = htop
        s% dx    = dx

        !---------------------------------------------------------
        ! loop over points along the slant within the model domain
        !---------------------------------------------------------
        points: do ip = s% Nmod, 1, -1

          call Grid_Indices                         &
            (s% SlantPointsEll(ip,1) * r2d,& ! <-- geodetic longitude
             s% SlantPointsEll(ip,2) * r2d,& ! <-- geodetic latitude
             grid,                         & ! <-- grid data type
             idx,                          & ! --> Grid point indices [Point, index]
             w,                            & ! --> Weight
             np                            ) ! --> number of points returned
          if (np == 0) then
            !---------------------------------------------
            ! abort if intermediate point is out of domain
            !---------------------------------------------
            ! => Does not work for ICON-D2
            !    For slants near the border some intermediate points can
            !    be out of the domain!
            ! => Check grid type ?


!!$            if (ip /= s% Nmod) then
!!$               write(*,*) 'ZYX ip, s% Nmod, s% SlantPointsEll(ip,:), s% obs%slant = ', &
!!$                    ip, s% Nmod, s% SlantPointsEll(ip,1:2)*r2d, s% SlantPointsEll(ip,3), &
!!$                    s% obs%slant
!!$               do ig=1, s% Nmod
!!$                  write(*,*) 'ZYX s% SlantPointsEll(ig,:) = ',                 &
!!$                       s% SlantPointsEll(ig,1:2)*r2d, s% SlantPointsEll(ig,3)
!!$               end do
!!$               call finish('process_std','SETUP_COLS: np==0,ip/=Nmod')
!!$            endif
            !--------------------------------------
            ! exit loop if out of horizontal domain
            !--------------------------------------
            Nmiss = Nmiss + 1
            s% Nhor                   = ip -1
            s% Naloc                  = 0
            s% idxi        (:,ip:)    = 0
            s% wi          (:,ip:)    = 0._wp
            s% SlantPointsGeo (ip:,:) = invalid
            call decr_use (obs% o% body(spot% o% i + is)% use, &
                           STAT_DISMISS, check=CHK_DOMAIN     )
            exit points
          else if (np/=nnghb) then
            !----------------------------------------------------
            ! abort if number of neighbour points is inconsistent
            !----------------------------------------------------
            write(0,*) dace% pe, is, ip,  'SETUP_COLS: np/=nnghb', np, nnghb, mp
            call finish('process_std','GPSGB: SETUP_COLS: np/=nnghb')
          endif

          !------------------------
          ! one model column only ?
          !------------------------
          if (order==1) then
            select case (abs(ztd_col))
            case (1)
              !------------------
              ! nearest neighbour
              !------------------
              i = maxloc (w(1:np), dim=1)
            case (2)
              !-------------------------------
              ! lowest model surface elevation
              !-------------------------------
              z = huge(z)
              do ig = 1, nnghb
                if (grid% hsurf (idx(ig,1),idx(ig,2),1,idx(ig,3)) < z) then
                  z = grid% hsurf (idx(ig,1),idx(ig,2),1,idx(ig,3))
                  i = ig
                endif
              end do
            case default
              call finish('process_std','GPSGB: SETUP_COLS, not implemented: ztd_col=')
            end select
            !-------------------------
            ! Keep single model column
            !-------------------------
            idx(1,:) = idx(i,:)     ! Column index and pe
            w  (:)   = 0.0_wp       ! Trivial weights
            w  (1)   = 1.0_wp
            np       = 1
          end if
          !------------------------------------------------------
          ! store model column info in the 3dvar/LETKF structures
          !------------------------------------------------------
          call add_index (idx, nngh3, obs% o, iatm, natm, spot% imcol,    &
                          nimcol, ixmcol, grid, spot% i_time, spot% w_time)
          if (nngh3==1) then
            s% idxi         (:,ip) = ixmcol (1)      ! pass nnghb times nearest neighbour
          else
            s% idxi         (:,ip) = ixmcol (:nngho)
          endif
          s% wi             (:,ip) = w      (:nngho)
        end do points

        spot% mke   = grid% nz
        !-----------------
        ! some diagnostics
        !-----------------
        if (tmp_stat) nullify (slants(is)% Station)
      end do   ! slant loop: do is = 1, size (slants)

      !--------------------
      ! memory deallocation
      !--------------------
      n = Nmiss - size(slants)
      call alloc_imcol (spot% imcol, nimcol)
      spot% n_spt = size (spot% imcol)
      call store_std (obs% o, spot, slants)
      call destruct                (slants)
      deallocate (slants)

      if (size(spot% imcol) == 0 .and. n == 0) then
         ! All slants leave the horizontal domain and
         ! have been set to STAT_OBS_ONLY => reject whole spot
         call decr_rpt_use (spot, CHK_DOMAIN, STAT_OBS_ONLY)
      else if (size(spot% imcol) == 0) then
         ! No grid columns found, this should not happen:
         ! Set state to DISMISS for whole spot, i.e. all slants and print warning
         write(*,'(3a,2(f8.2,tr2))')                                            &
              'mo_std, process_std, SETUP_COLS: No columns found for station ', &
              trim(spot%statid), ', at lon/lat ',                               &
              spot% col% c% dlon, spot% col% c% dlat
         call decr_rpt_use (spot, CHK_DOMAIN, STAT_DISMISS)
      end if

      if (tsk==0) return
    endif ! TSK_SETUP_COLS

    !--------------------------------
    ! set observation characteristics
    !--------------------------------
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      tsk=tsk-TSK_SET_CHR
      !------------------------
      ! already set in TSK_READ
      !------------------------
      if (tsk==0) return
    endif

    !--------------------------------
    ! setup description of PSAS-space
    !--------------------------------
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      tsk=tsk-TSK_SETUP_FUL0
      if (dace% pe==obs% o% pe) then
        m = 3 * atm% grid% nz
        call new_int (obs% o, spot, m * size(spot% imcol))
      endif
      if (tsk==0) return
    endif

    !========================================================
    ! set interpolation space: observed quantities and levels
    !========================================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      tsk=tsk-TSK_SETUP_FULL
      if (dace% pe==obs% o% pe) then
        m = 3 * atm% grid% nz
        do k=0,size(spot% imcol)-1
          i = spot% imcol(k+1)% imc(1)
          n = k*m
          obs% o% t_int (spot%i%i+n+1:spot%i%i+n+m:3) = OBS_TV
          obs% o% t_int (spot%i%i+n+2:spot%i%i+n+m:3) = OBS_RH
          obs% o% t_int (spot%i%i+n+3:spot%i%i+n+m:3) = OBS_H
          obs% o% lev   (spot%i%i+n+1:spot%i%i+n+m:3) = cols% col(i)% p
          obs% o% lev   (spot%i%i+n+2:spot%i%i+n+m:3) = cols% col(i)% p
          obs% o% lev   (spot%i%i+n+3:spot%i%i+n+m:3) = cols% col(i)% p
        end do
      endif
      if(tsk==0) return
    endif

    !----------------------------------------------
    ! TSK_Y: run the nonlinear observation operator
    ! TSK_K: in addition provide the Jacobi matrix
    !----------------------------------------------
    if (iand   (TSK_Y+TSK_K, tsk) /= 0) then
      adjoint = iand (TSK_K, tsk) /= 0
      grid => atm% grid
      if (adjoint) then
         select case (grid% vct)
         case (VCT_Z_GEN) ! ICON
         case default
            call finish('process_std',                                  &
                        'GPSGB: TSK_K not implemented, vct: '//char1(grid% vct))
         end select
      end if
      if (spot% pe_eval == dace% pe) then

        !--------------------------------------
        ! load the STD specific slant meta data
        !--------------------------------------
        call load_std (obs% o, spot, slants)

        !--------------------------------------------------------------
        ! derive (temporary) observation info if not present (BUFR ZTD)
        !--------------------------------------------------------------
        tmp_stat = .not. associated (slants (1)% Station)

        if (tmp_stat) then
           call spot2station (spot, atm% grid, slants, station)
           !----------------------------------------------------------------
           ! interpolate geoid correction and compute height above ellipsoid
           !----------------------------------------------------------------
           if (slants(1)%GeoidCorr == invalid) then
              ! The geoid undulation was not available and spot2station could
              ! not compute correct ellipsoidal and cartesion station
              ! coordinates. This is now done using the geoid.

              if (slants(1)%Nnghb == 1) then
                 ! ZTD with nearest neighbour, no interpolation
                 i  = slants(1)%idxi(1,1)
                 ic =  spot% imcol(i)% imc(1)
                 c  => cols% col(ic)
                 station% CoordEll(3) = station% CoordEll(3) + c%s%geoid
              else
                 ! ZTD or STD, interpolate several columns to GNSS station
                 do i = 1, slants(1)%Nnghb
                    k  = slants(1)%idxi(i,1)
                    ic =  spot% imcol(k)% imc(1)
                    c  => cols% col(ic)
                    station% CoordEll(3) = station% CoordEll(3) +         &
                                            slants(1)%wi(i,1) * c%s%geoid
                 end do
              end if

              !----------------------------------------------------------------
              ! transform ellipsoidal coordinates to ECEF Cartesian coordinates
              !----------------------------------------------------------------
              call  Ellips2Cart( station% CoordEll (1),& ! in : longitude
                                 station% CoordEll (2),& ! in : latitude
                                 station% CoordEll (3),& ! in : altitude
                                 station% CoordCart(1),& ! out: X
                                 station% CoordCart(2),& ! out: Y
                                 station% CoordCart(3),& ! out: Z
                                 WGS84Param            ) ! in : ref. ellipsoid

           end if     ! if (slants(1)%GeoidCorr == -999._sp) then
        end if        ! if (tmp_stat) then

        ke = cols% ke
        if (adjoint) then
           allocate (gh_a(ke), rh_a(ke), q_a(ke), t_a(ke), p_a(ke))
           allocate (gh  (ke), rh  (ke), q  (ke), t  (ke), p  (ke))
!          allocate (tv  (ke))
           allocate (geo (ke), lnp (ke))
           allocate (d2t  (ke),    d2q  (ke),    d2lnp (ke)   )
           allocate (t_z  (ke),    q_z  (ke),    p_z   (ke))
#ifdef DEBUG_STD_AD
           allocate (d2t_z(ke,ke), d2q_z(ke,ke), d2p_z (ke,ke))
           allocate (d2t_t(ke,ke), d2q_q(ke,ke), d2p_p (ke,ke))
           allocate (tmp1  (ke),    tmp2  (ke))
           allocate (t_int (ke),    q_int (ke),    p_int  (ke)   )
           allocate (t_by_z(ke,ke), q_by_z(ke,ke), p_by_z (ke,ke))
           allocate (t_by_t(ke,ke), q_by_q(ke,ke), p_by_p (ke,ke))
#endif
        end if
        !-----------------
        ! loop over slants
        !-----------------
        do i = 1, size (slants)
          s => slants (i)
          if (s% Naloc == 0) then
            !-----------------------------------------------------
            ! part of slant is out of domain, return invalid value
            !-----------------------------------------------------
            y%x (spot%o%i+i) = invalid
          else
            !-----------------------------------------------------
            ! gather the model columns from the 3dvar derived type
            !-----------------------------------------------------
            if (tmp_stat) slants(i)% Station => station
            allocate   (mc(s% Nnghb, s% Naloc))
            if (adjoint) then
              allocate (h (s% Nnghb, s% Naloc))
            else
              nullify (h)
            endif
            call set_std_data (mc, s, atm% grid, spot, cols, xi, colerr)
            if (adjoint) then
              do k = 1, s% Naloc
                do j = 1, s% Nnghb
                  l         = s% idxi (j,k)
                  ic        = spot% imcol(l)% imc(1)
                  allocate (h(j,k)% t (size(cols% col(ic)% geo)))
                  allocate (h(j,k)% p (size(cols% col(ic)% p  )))
                  allocate (h(j,k)% q (size(cols% col(ic)% q  )))
                  h(j,k)% t = 0._wp
                  h(j,k)% p = 0._wp
                  h(j,k)% q = 0._wp
                end do
              end do
            endif

            if (colerr == 0) then
               !-------------------------------------------------------------
               ! interpolate geoid correction and compute height above geooid
               !-------------------------------------------------------------
               if (slants(1)%GeoidCorr == invalid) then
                  ! Re-compute signal path with correct height above ellipsoid
                  call std_path_coord (s, s%htop, err, grid%model, grid%nz, s%dx, 10)
               end if

               if (.not. associated(s%SlantPointsGeo)) &
                    allocate(s%SlantPointsGeo(1:s% Ntot,1:3))
               s%SlantPointsGeo = s%SlantPointsEll

               points2: do ip = 1, s% Nmod
                  do j = 1, s% Nnghb
                     l = s% idxi(j,ip)
                     ic =  spot% imcol(l)% imc(1)
                     c  => cols% col(ic)
                     s%SlantPointsGeo(ip,3) = s%SlantPointsGeo(ip,3)     &
                                              - s%wi(j,ip) * c%s%geoid
                  end do
               end do points2

               !------------------------------------
               ! Call the Slant Total Delay operator
               !------------------------------------
               call std_delay (s, mc, delay, adjoint, h, err)

            else  ! colerr = 1
               ! Invalid columns, operator could not be called, reject slant
               delay = invalid
               err   = 10
            end if
            !---------------------------------
            ! Error codes:
            !   1: wrong number of columns
            !   2: out of domain
            !   3: station below model surface
            !   4: station too high
            !   5: ray path cuts model surface
            !   6: computed delay out of range
            !   ...
            !  10: invalid column data in set_std_data
            !---------------------------------
            select case (err)
            case (0)
               ! OK
            case (1,2)
               call decr_use (obs% o% body(spot% o% i + i)% use, check=CHK_DOMAIN)
            case (3)
               call decr_rpt_use (spot, CHK_HEIGHT, STAT_REJECTED, &
                                  comment='below model surface'    )
            case (4)
               call decr_rpt_use (spot, CHK_HEIGHT, STAT_REJECTED, &
                                  comment='above model surface'    )
            case (5)
               call decr_use (obs% o% body(spot% o% i + i)% use, check=CHK_HEIGHT)
            case (6)
               call decr_use (obs% o% body(spot% o% i + i)% use, check=CHK_OPERATOR)
            case (10)
               call decr_use (obs% o% body(spot% o% i + i)% use, check=CHK_DOMAIN)
            end select
            !-----------------
            ! store the result
            !-----------------
            if (iand (TSK_Y,tsk) /= 0) then
               y%x (spot%o%i+i) = delay
               if (pl_method == 20 .and. obs% o% body(spot%o%i+i)% plev < 0) then
                  if (delay /= invalid) then
                     obs% o% body(spot%o%i+i)% plev = s%Obs%Press
                     obs% o% body(spot%o%i+i)% lon  = s%Obs%PosEll(1)*r2d
                     obs% o% body(spot%o%i+i)% lat  = s%Obs%PosEll(2)*r2d
                     ! -180 < lon < 180
                     if (obs% o% body(spot%o%i+i)% lon > 180.0_sp) then
                        obs% o% body(spot%o%i+i)% lon = &
                             obs% o% body(spot%o%i+i)% lon - 360.0_sp
                     end if
                  end if
               end if
            end if
            !------------------
            ! store the adjoint
            !------------------
            if (adjoint) then

               obs% yi% x (spot% o%i+i) = delay
               allocate (std_t(ke,size (spot% imcol)))
               allocate (std_p(ke,size (spot% imcol)))
               allocate (std_q(ke,size (spot% imcol)))
               std_t = 0._wp
               std_p = 0._wp
               std_q = 0._wp
               !----------------------------------------------------------
               ! Accumulate contributions to adjoint for each model column
               !----------------------------------------------------------
               do k = 1, s% Naloc
                  do j = 1, s% Nnghb
                     l          = s% idxi (j,k)
                     std_t(:,l) = std_t(:,l) + h(j,k)% t
                     std_p(:,l) = std_p(:,l) + h(j,k)% p
                     std_q(:,l) = std_q(:,l) + h(j,k)% q
                  end do
               end do
!!$               print *, "### ____________________________"
!!$               print *, "### spot% z    =", real (spot% z)
!!$               print *, "### geoid( 1 ) =", real (mc(1,1)% geoid)
!!$               print *, "### gpm  ( 1 ) =", real (mc(1,1)% gpm)
!!$               print *, "### t    ( 1 ) =", real (mc(1,1)% t)
!!$               print *, "### max(std_t) =", maxval (abs (std_t))
!!$               print *, "### max(std_p) =", maxval (abs (std_p))
!!$               print *, "### max(std_q) =", maxval (abs (std_q))
!!$               print *, "### std_t(:,1) =", real (std_t(:,1))
!!$               print *, "### std_t(:,2) =", real (std_t(:,2))
!!$               print *, "### std_t(:,3) =", real (std_t(:,3))
!!$               print *, "### std_p(:,1) =", real (std_p(:,1))
!!$               print *, "### std_p(:,2) =", real (std_p(:,2))
!!$               print *, "### std_p(:,3) =", real (std_p(:,3))
!!$               print *, "### std_q(:,1) =", real (std_q(:,1))
!!$               print *, "### std_q(:,2) =", real (std_q(:,2))
!!$               print *, "### std_q(:,3) =", real (std_q(:,3))

               !--------------------------------------------
               ! The following code is inspired by mo_occ_1d
               !--------------------------------------------
               nn = 3 * ke
               !---------------
               ! Jacobi matrix:
               !---------------
               allocate (Hstd(nn))
               Hstd = 0._wp
               allocate (delay_t(ke))
               allocate (delay_q(ke))
               allocate (delay_z(ke))

               do l = 1, size (spot% imcol)
                  n  = (l-1) * nn
                  ix = spot% imcol(l)% imc(1)
                  q  =     cols% col(ix)% q
                  t  =     cols% col(ix)% t
                  geo=     cols% col(ix)% geo
                  lnp=     cols% col(ix)% p
                  p  = exp(cols% col(ix)% p)
                  rh = rh_q (q, t, p)
!                 tv = t * (1 + Eps * q)
!!                tv = t * (1 + Eps * q - qx)   ! incl. condensates
                  !----------------------
                  ! temporary: set t,q,ps
                  !----------------------
                  q_a = 1._wp; rh_a=0._wp; t_a=0._wp; p_a=0._wp
                  call q_rh_adj (q_a, rh_a, t_a, p_a, rh, t, p)
                  !---------------------------------
                  ! account for generalised humidity
                  !---------------------------------
                  gh   =  gh_rh (rh, .false.)
                  call rh_gh (rh, gh, gh_a)
                  rh_a = rh_a * gh_a
                  !====================================================
                  ! Interpolation to domain of STD operator (fixed geo)
                  !====================================================
                  call init_spline (geo, t,   d2t  )
                  call init_spline (geo, q,   d2q  )
                  call init_spline (geo, lnp, d2lnp)
                  !--------------------------------------------------------
                  ! Optimized version of special case for adjoint of spline
                  !--------------------------------------------------------
                  do j = 1, ke
                     j1 = j+1; if (j == ke) j1 = j-1
                     dz      = (geo(j1) - geo(j))
                     a2      =   0.5_wp * d2t(j)
                     a3      = (d2t(j1) - d2t(j)) / (6*dz)
                     a1      = (t  (j1) - t  (j)) / dz - dz*(a2 + dz*a3)
                     t_z (j) = - a1
                     a2      =   0.5_wp * d2q(j)
                     a3      = (d2q(j1) - d2q(j)) / (6*dz)
                     a1      = (q  (j1) - q  (j)) / dz - dz*(a2 + dz*a3)
                     q_z (j) = - a1
                     a2      =     0.5_wp * d2lnp(j)
                     a3      = (d2lnp(j1) - d2lnp(j)) / (6*dz)
                     a1      = (lnp  (j1) - lnp  (j)) / dz - dz*(a2 + dz*a3)
                     p_z (j) = - a1 * p(j)
                  end do
#ifdef DEBUG_STD_AD
                  !------------------------
                  ! "Full" (slower) version
                  !------------------------
                  call init_spline_adj (geo, t,   d2t,   d2t_z, d2t_t)
                  call init_spline_adj (geo, q,   d2q,   d2q_z, d2q_q)
                  call init_spline_adj (geo, lnp, d2lnp, d2p_z, d2p_p)
                  do j = 1, ke
                     call spline_adj (geo, t,   d2t,   d2t_z, d2t_t, geo(j), &
                                      t_int(j), t_by_z(:,j), t_by_t(:,j),    &
                                      tmp, tmp1, tmp2)
                     call spline_adj (geo, q,   d2q,   d2q_z, d2q_q, geo(j), &
                                      q_int(j), q_by_z(:,j), q_by_q(:,j),    &
                                      tmp, tmp1, tmp2)
                     call spline_adj (geo, lnp, d2lnp, d2p_z, d2p_p, geo(j), &
                                      p_int(j), p_by_z(:,j), p_by_p(:,j),    &
                                      tmp, tmp1, tmp2)
                     !-------------------------------
                     ! dp(j)/dz = p(j)*d(ln(p(j)))/dz
                     !-------------------------------
                     p_by_z(:,j) = p_by_z(:,j) * exp (p_int(j))
                     !-------------------------------------------------------
                     ! Verify that t_by_z, q_by_z, and p_by_z are diagonal,
                     ! and that t_by_t, q_by_q, and p_by_p are unit matrices:
                     !-------------------------------------------------------
                     tmp1    = 0._wp
                     tmp1(j) = 1._wp
                     if (maxval (abs (t_by_t(:,j) - tmp1)) > 1.e-15_wp) then
                        print *, "### j,t_by_t=",j,real(t_by_z(:,j))
                     end if
                     if (maxval (abs (t_by_z(:,j) - t_z(j)*tmp1)) > abs (t_z(j))*1.e-15_wp) then
                        print *, "### j,t_z   =",j,real(t_z (j)), abs (t_z (j)-t_by_z(j,j))
                        print *, "### j,t_by_z=",j,real(t_by_z(:,j))
                     end if
                     if (maxval (abs (q_by_q(:,j) - tmp1)) > 1.e-15_wp) then
                        print *, "### j,q_by_q=",j,real(q_by_q(:,j))
                     end if
                     if (maxval (abs (q_by_z(:,j) - q_z(j)*tmp1)) > abs (q_z(j))*1.e-15_wp) then
                        print *, "### j,q_z   =",j,real(q_z (j)), abs (q_z (j)-q_by_z(j,j))
                        print *, "### j,q_by_z=",j,real(q_by_z(:,j))
                     end if
                     if (maxval (abs (p_by_p(:,j) - tmp1)) > 1.e-15_wp) then
                        print *, "### j,p_by_p=",j,real(p_by_p(:,j))
                     end if
                     if (maxval (abs (p_by_z(:,j) - p_z(j)*tmp1)) > abs (p_z(j))*1.e-15_wp) then
                        print *, "### j,p_z   =",j,real(p_z (j)), abs (p_z (j)-p_by_z(j,j))
                        print *, "### j,p_by_z=",j,real(p_by_z(:,j))
                     end if
                  end do
                  !========================================================
                  ! Jacobian for change of basis via vertical interpolation
                  !========================================================
!                 delay_t  = matmul (t_by_t, std_t(:,l))
!                 delay_q  = matmul (q_by_q, std_q(:,l))
!                 delay_z  = matmul (t_by_z, std_t(:,l)) + &
!                            matmul (q_by_z, std_q(:,l)) + &
!                            matmul (p_by_z, std_p(:,l))
#endif
                  !========================================================
                  ! Jacobian for change of basis via vertical interpolation
                  !========================================================
                  delay_t  = std_t(:,l)
                  delay_q  = std_q(:,l)
                  delay_z  = std_t(:,l) * t_z(:) + &
                             std_q(:,l) * q_z(:) + &
                             std_p(:,l) * p_z(:)
                  !
                  Hstd(1:3*ke:3) = delay_t(:) +   &
                                   delay_q(:) * t_a
                  Hstd(2:3*ke:3) = delay_q(:) * rh_a
                  Hstd(3:3*ke:3) = delay_z(:)
                  !
                  k = obs% H% ia (spot% i% i + n + 1)
                  do j=1,nn
                     obs% H% ia (spot% i% i + n + j) = k
                     ! Note: spot% o% n == 1
                     obs% H% packed (k) = Hstd (j)
                     obs% H% ja     (k) = spot% o% i + 1
                     k = k + 1
                  end do
                  obs% H% ia (spot% i% i + n + nn + 1) = k
               end do
               deallocate (Hstd)
               deallocate (delay_t, delay_q, delay_z)
               deallocate (std_t, std_p, std_q)
            end if
            !---------
            ! clean up
            !---------
            do k = 1, s% Naloc
              do j = 1, s% Nnghb
                deallocate (mc(j,k)% gpm)
                deallocate (mc(j,k)% p)
                deallocate (mc(j,k)% q)
                deallocate (mc(j,k)% t)
                if (adjoint) then
                  deallocate (h(j,k)% t)
                  deallocate (h(j,k)% p)
                  deallocate (h(j,k)% q)
                endif
              end do
            end do
            deallocate (mc)
            if (adjoint ) deallocate (h)
            if (tmp_stat) nullify (slants(i)% Station)
          endif
        end do
        call destruct (slants)
        deallocate    (slants)
        if (adjoint) then
           deallocate (gh  , rh  , q  , t  , p  )
           deallocate (gh_a, rh_a, q_a, t_a, p_a)
           deallocate (lnp , geo)
!          deallocate (tv)
           deallocate (d2t,   d2q,   d2lnp)
           deallocate (t_z,   q_z,   p_z  )
#ifdef DEBUG_STD_AD
           deallocate (d2t_z, d2q_z, d2p_z)
           deallocate (d2t_t, d2q_q, d2p_p)
           deallocate (tmp1, tmp2)
           deallocate (t_int, q_int, p_int, t_by_z, q_by_z, p_by_z, t_by_t, q_by_q, p_by_p)
#endif
        endif
      endif
      if (iand (TSK_Y, tsk) /= 0) tsk = tsk - TSK_Y
      if (iand (TSK_K, tsk) /= 0) tsk = tsk - TSK_K

      if (tsk==0) return
    endif

    !-----------------------------------------------
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !-----------------------------------------------
    if (iand (TSK_R,tsk) /= 0) then
      if (obs% o% pe == dace% pe) then
        call load_std (obs% o, spot, slants)
        !------------------------------------------
        ! setup observation error covariance matrix
        ! setup variational quality control bounds
        !------------------------------------------
        call get_rule (type     = spot% hd% modtype,  &! <- module      type
                       obstype  = spot% hd% obstype,  &! <- observation type
                       codetype = spot% hd% codetype, &! <- code type
                       bf_type  = iud,                &! <- no BUFR     type
                       bf_subt  = iud,                &! <- no BUFR  subtype
                       db_kz    = iud,                &! <- no Datenbankkennzahl
                       stname   = '',                 &! <- no Station Name
                       lat      = spot% col% c% dlat, &! <- latitude
                       lon      = spot% col% c% dlon, &! <- longitude
                       o        = set            )     ! -> channel information
        if (.not. associated (obs% o% s_vqc)) then
          allocate (obs% o% s_vqc (obs% o% n_obs))
          obs% o% s_vqc = svqc
        endif
        if (.not. associated (obs% o% f_vqc)) then
          allocate (obs% o% f_vqc (obs% o% n_obs))
          obs% o% f_vqc = vqc_form
        endif
        if (set% sgm_vq /= rud) then
          obs% o% s_vqc (spot% o% i+1 : spot% o% i + spot% o% n) = set% sgm_vq
        endif
        if (set% frm_vq /= iud) then
          obs% o% f_vqc (spot% o% i+1 : spot% o% i + spot% o% n) = set% frm_vq
        endif
        n = spot% o% n
        i = spot% o% i
        k = obs% R% ia (i + 1)
        do j=1,n
          obs% R% ia (i + j) = k
          obs% R% packed (k) = ZTDerror ** 2
          obs% R% ja (k) = i + j
          k = k + 1
        end do
        obs% R% ia (i + n + 1) = k
        call destruct (slants)
        deallocate (slants)
      endif
      tsk = tsk - TSK_R
      if (tsk == 0) return
    endif

    !------------------------------------------
    ! TSK_SHRINK:
    ! release unused observations in the report
    !------------------------------------------
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change, mask=msk)
      if (change) then
        ! currently not yet implemented...
!       call load_std (obs% o, spot, slants)
!       ....
!       call store_std (obs% o, spot, slants)
!       call destruct                (slants)
!       deallocate (slants)
        deallocate (msk)
      endif
      tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !==========================
    ! abort if any task is left
    !==========================
    if (tsk /= 0) then
      write(0,*)  'process_std:  TSK =',tsk
      call finish('process_std','GPSGB: TSK_?')
    endif

  end subroutine process_std
!------------------------------------------------------------------------------
  subroutine spot2station (spot, grid, slants, station)
  type (t_spot)      ,intent(in)  :: spot      ! 3dvar report meta data
  type (t_grid)      ,intent(in)  :: grid      ! model grid meta data
  type(SlantData)    ,intent(in)  :: slants(:) ! slant data for one station
  type (GNSSStation) ,intent(out) :: station   ! STD station meta data
  !---------------------------------------------------------
  ! derive STD station meta data from 3dvar report meta data
  !---------------------------------------------------------
    integer  :: idx (mp,4)
    integer  :: np
    real(wp) :: w     (mp)
    integer  :: i
    !-----------
    ! station id
    !-----------
    station% name        = spot% statid
    station% SName       = spot% statid
    station% ID          = 0
    station% Country     = ''
    !-------------------------
    ! geographical coordinates
    !-------------------------
    station% CoordGeo(1) = spot% col% c% dlon * d2r
    station% CoordGeo(2) = spot% col% c% dlat * d2r
    station% CoordGeo(3) = spot% z
    !-----------------------------------
    ! initialize ellipsoidal coordinates
    !-----------------------------------
    station% CoordEll(1) = station% CoordGeo(1)
    station% CoordEll(2) = station% CoordGeo(2)
    station% CoordEll(3) = station% CoordGeo(3)
    !------------------------------------------------------
    ! compute height above ellipsoid using geoid undulation
    !------------------------------------------------------
    if (slants(1)%GeoidCorr /= invalid) then
       station% CoordEll(3) = station% CoordEll(3) + slants(1)%GeoidCorr
    end if

    ! Warning
    ! Ellipsoidal coordinates can be computed only if the geoid undulation
    ! is available (slants(1)%GeoidCorr /= -999._sp).
    ! Without the geoid undulation the geographical coordinates will
    ! be used, i.e. the ellipsoidal and cartesian coordinates are
    ! somewhat wrong. They are still needed in order to find an
    ! approximation of the signal path and to get the corresponding
    ! model columns.

    !----------------------------------------------------------------
    ! transform ellipsoidal coordinates to ECEF Cartesian coordinates
    !----------------------------------------------------------------
     call  Ellips2Cart( station% CoordEll (1),& ! in : longitude
                        station% CoordEll (2),& ! in : latitude
                        station% CoordEll (3),& ! in : altitude
                        station% CoordCart(1),& ! out: X
                        station% CoordCart(2),& ! out: Y
                        station% CoordCart(3),& ! out: Z
                        WGS84Param            ) ! in : reference ellipsoid
  end subroutine spot2station
!==============================================================================
  subroutine read_std_netcdf (ifile, i_source, obs, head, rec1, recl, &
                              lkeep, nkeep,cc)
  !===================================================================
  ! Read ZTD and/or STD observations from netCDF (converted BUFR data)
  ! BUFR data are encoded using table D sequence 307022 or 307024
  !===================================================================
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data to set
  type (t_head) ,intent(in)           :: head      ! header data already encoded
  integer       ,intent(in)           :: rec1      ! first record to read
  integer       ,intent(in)           :: recl      ! last  record to read
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.
  integer       ,intent(in) ,optional :: cc        ! part of year ccyy

  integer :: status, varid

  ! Find table D sequence
  ! => Access variable defined only in 307024: Geoid undulation, NGEUN
  status = nf90_inq_varid(ncid, 'NGEUN', varid)
  if (status == NF90_NOERR) then
     ! NGEUN exists => 307024 (new STD/ZTD BUFR)
     call read_netcdf_307024 (ifile, i_source, obs, head, rec1, recl, &
                              lkeep, nkeep,cc)
  else
     ! NGEUN does not exist => 307022 (old ZTD BUFR used by E-GVAP)
     call read_netcdf_307022 (ifile, i_source, obs, head, rec1, recl, &
                              lkeep, nkeep,cc)
  end if

end subroutine read_std_netcdf


!==============================================================================
  subroutine read_netcdf_307024 (ifile, i_source, obs, head, rec1, recl, &
                                 lkeep, nkeep,cc)
  !===============================================================
  ! Read STD or ZTD observations from netCDF (converted BUFR data)
  ! BUFR data are encoded using table D sequence 307024
  ! => New ZTD/STD BUFR
  !===============================================================
  ! Warning
  ! Only a subset of the BUFR descriptors is read.
  ! The current implementation can replace the old ZTD BUFR and
  ! the ASCII STD data provided by the GFZ. New features are
  ! not yet implemented!
  !===============================================================
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data to set
  type (t_head) ,intent(in)           :: head      ! header data already encoded
  integer       ,intent(in)           :: rec1      ! first record to read
  integer       ,intent(in)           :: recl      ! last  record to read
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.
  integer       ,intent(in) ,optional :: cc        ! part of year ccyy

    !================
    ! local variables
    !================
    type (t_use)             :: use             ! status variable
    type (t_head)            :: hd              ! report header
    type (t_spot)            :: spt0, spti      ! report meta data
    type (t_spot) ,save      :: empty           !
    type (SlantData), allocatable :: slants (:)    ! slant observations
    integer                  :: status          ! NetCDF status variable
    integer                  :: ncvars          ! NetCDF number of variables  defined in this NetCDF file
    integer                  :: len_report      ! number of reports in NetCDF file
    integer                  :: len_slants      ! number of reports in NetCDF file
    integer                  :: i, j, k         ! loop index
    integer                  :: is              ! report loop index
!   integer                  :: nc0             ! dimension for netcdf getvar 1-dimensional arrays
!   integer                  :: nc1             ! first  dimension in start / count
!   integer                  :: nc2             ! second dimension in start / count
    character(NF90_MAX_NAME) :: yname_v         ! NetCDF variable name
    logical                  :: lpr_std         ! std   reports from netcdf are printed
    logical                  :: lpr_extd        ! extended  printing of std
    integer                  :: npr_extd        ! number of extended  printing of temps
    !------------ GNSS station ------------------------------------------------
    integer                  :: i_met_YSOSN = 0 ! index for station id
    integer                  :: i_met_DATE  = 0 !   date
    integer                  :: i_met_TIME  = 0 !   time
    integer                  :: i_met_MLAH  = 0 !   latitude
    integer                  :: i_met_MLOH  = 0 !   longitude
    integer                  :: i_met_MHP   = 0 !   station height
    integer                  :: i_met_NGEUN = 0 !   geoid undulation
    !------------ ZTD observation ---------------------------------------------
    integer                  :: i_met_NZPDNA= 0 !   atmospheric path delay
    integer                  :: i_met_NEENAZ= 0 !   estimated error in path delay
    integer                  :: i_met_MEQCGD= 0 !   quality flags
    !------------ STD observation ---------------------------------------------
    integer                  :: i_met_MPTID = 0 !   satellite PRN
    integer                  :: i_met_MDA   = 0 !   azimuth
    integer                  :: i_met_MDE   = 0 !   elevation
    integer                  :: i_met_NPDNA = 0 !   STD
    !                                 NZPDNA    !   STD mapped into zenith => NZPDNA
    !--------------------------------------------------------------------------
    integer                  :: numDims         ! number of dimensions for NetCDF variable
    integer                  :: numAtts         ! number of attributes for NetCDF variable
    integer                  :: n               ! number of reports to read
    integer                  :: dimid_slants    !
    integer                  :: dimid_report    !
    integer                  :: varid_NADES     ! variable id
    integer                  :: entry1,entry    ! position in source file (subset)
    ! data from netCDF
    character (12), allocatable :: ypcid(:)      ! GNSS product string
    character (1),  allocatable :: ypcidarr(:,:) ! GNSS product string for neCDF read
    real(sp), allocatable    :: mhp    (:)       ! station height
    real(sp), allocatable    :: ngeun  (:)       ! geoid undulation
    real(sp), allocatable    :: nzpdna (:,:)     ! ZTDs
    real(sp), allocatable    :: neenaz (:,:)     ! ZTD errors
    integer,  allocatable    :: meqcgd (:,:)     ! Quality flags GNSS
    integer,  allocatable    :: mptid  (:,:)     ! satellite PRN
    integer,  allocatable    :: msacl  (:,:)     ! satellite system specification
    real(sp), allocatable    :: mda    (:,:)     ! bearing or azimuth
    real(sp), allocatable    :: mde    (:,:)     ! elevation
    real(sp), allocatable    :: npdna  (:,:)     ! STDs

!!$    integer,  allocatable    :: nqfgd (:)       ! quality flags
!!$    real(sp), allocatable    :: mhp   (:)       ! station height
!!$    real(sp), allocatable    :: mda   (:)       ! bearing or azimuth
!!$    real(sp), allocatable    :: mde   (:)       ! elevation
!!$    real(sp), allocatable    :: nades (:)       ! atmospheric path delay
!!$    real(sp), allocatable    :: neerr (:)       ! estimated error in path delay

    logical                  :: lk              ! flag to keep this report
    character(len=4)         :: cnt_prd         ! code for center/product
    character(len=10)        :: station         ! possibly modified station name
    integer                  :: center          ! center id within DACE (not used so far)
    integer                  :: dim_loop_id     ! index of dimension
    integer                  :: dim_records= -1 ! BUFR_records in netCDF file
    integer                  :: dim_YPCID= -1   ! sring length of YPCID
    integer                  :: dim_loop_0 = -1 ! = 1 if ZTDs available
    integer                  :: dim_loop_1 = -1 ! Number of ZTDs
    integer                  :: dim_loop_2 = -1 ! = 2 for gradients
    integer                  :: dim_loop_3 = -1 ! = 1 if STDs available
    integer                  :: dim_loop_4 = -1 ! Number of STDs
    integer                  :: dim_loop_5 = -1 ! = 2 for gradients
    integer                  :: varid
    integer                  :: iascii_ymissing ! ASCII code of ymissing

    integer                  :: dimids (dimids_max) ! NetCDF dimension ids
    integer                  :: Nstd            ! number of valid slants in report

    type (SlantTotalDelay), allocatable  :: STD(:)
    integer, allocatable                 :: UseSTD(:)

    lpr_std  = .false.
    if (netcdf_verb > 1) lpr_std = .true.
    ! debug
    !lpr_extd = .TRUE.
    !lpr_std  = .TRUE.

    !------------------------------
    ! get default number of reports
    !------------------------------
    len_report =  recl - rec1 + 1
!   nc0        = len_report

    lkeep = .false.
    nkeep = 0

    ! Compare ASCII code instead of strings
    iascii_ymissing = iachar(ymissing)

    !-----------------------------------------
    ! get dimension of record and replications
    !-----------------------------------------
    ! BUFR_records    = number of records in netCDF file
    ! Loop_000_maxlen = 1 if ZTDs available
    ! Loop_001_maxlen = max. Number of ZTDs in any record
    ! Loop_002_maxlen = 2 for gradients
    ! Loop_003_maxlen = 1 if STDs available
    ! Loop_004_maxlen = max. Number of STDs in any record
    ! Loop_005_maxlen ! = 2 for gradients

    status = nf90_inq_dimid(ncid, 'BUFR_records', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_records)
    end if
    status = nf90_inq_dimid(ncid, 'YPCID_strlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_YPCID)
    end if
    status = nf90_inq_dimid(ncid, 'Loop_000_maxlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_loop_0)
    end if
    status = nf90_inq_dimid(ncid, 'Loop_001_maxlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_loop_1)
    end if
    status = nf90_inq_dimid(ncid, 'Loop_002_maxlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_loop_2)
    end if
    status = nf90_inq_dimid(ncid, 'Loop_003_maxlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_loop_3)
    end if
    status = nf90_inq_dimid(ncid, 'Loop_004_maxlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_loop_4)
    end if
    status = nf90_inq_dimid(ncid, 'Loop_005_maxlen', dim_loop_id)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dim_loop_id, len=dim_loop_5)
    end if

    if ( lpr_std .and. lpr_extd ) then
       write (6,'(a,i3)') 'Number records   (BUFR_records)   : ', dim_records
       write (6,'(a,i3)') 'ZTD replication  (Loop_000_maxlen): ', dim_loop_0
       write (6,'(a,i3)') 'ZTD replications (Loop_001_maxlen): ', dim_loop_1
       write (6,'(a,i3)') 'ZTD Gradients    (Loop_002_maxlen): ', dim_loop_2
       write (6,'(a,i3)') 'STD replication  (Loop_003_maxlen): ', dim_loop_3
       write (6,'(a,i3)') 'STD replications (Loop_004_maxlen): ', dim_loop_4
       write (6,'(a,i3)') 'STD Gradients    (Loop_005_maxlen): ', dim_loop_5
    end if

    !------------------
    ! Step 1: Read ZTDs
    !------------------
    if (dim_loop_0 == 1) then

       ! WARNING
       ! Reading ZTDs from the new BUFR is not yet supported!
       ! => The code is incomplete and has not been tested.
       ! => The dbkz is related to the observation code type, i.e.
       !    a given dbkz is related to ZTDs or STDs but not to both.
       ! => dbkz = 95 is used only for STDs, a dbkz for  ZTDs and the
       !    new BUFR is not defined.

       !------------------------
       ! get dimension of fields
       !------------------------
       ! ZTD related data with three dimensions:
       ! NZPDNA(BUFR_records, Loop_000_maxlen, Loop_001_maxlen)
       ! => BUFR_records        - number of records in file
       ! => Loop_000_maxlen = 1 - ZTD data available
       ! => Loop_001_maxlen     - max. number of ZTDs per record

       len_report  = dim_records
       len_slants  = dim_loop_1

       !----------------------------------------
       ! define number of reports in source-file
       !----------------------------------------
       i_source = len_report

       !----------------
       ! allocate arrays
       !----------------
       n = min (len_report, rept_use(OT_GPSGB)% max_proc)

       allocate( mhp    (n)            )   ! Station height
       allocate( ngeun  (n)            )   ! Geoid undulation
       allocate( nzpdna (n,len_slants) )   ! ZTD
       allocate( neenaz (n,len_slants) )   ! ZTD error

       !allocate( (n,len_slants) )

       !--------------------------------
       ! read variables from NetCDF file
       !--------------------------------
       call get_real (mhp,   'MHP'  ,  invalsp, (/  n/),            (/  rec1/))
       call get_real (ngeun, 'NGEUN',  invalsp, (/  n/),            (/  rec1/))
       call get_real (nzpdna,'NZPDNA', invalsp, (/len_slants,1,n/), (/1,1,rec1/))
       call get_real (neenaz,'NEENAZ', invalsp, (/len_slants,1,n/), (/1,1,rec1/))

       !-------------------------------
       ! preset total number of reports
       !-------------------------------
       entry = sum (source(1:ifile-1)% entries) + rec1 - 1

       !------------------
       ! loop over reports
       !------------------
       do is = 1, n

          entry1  = entry    + 1
          entry   = entry    + 1

          !---------------
          ! initialize use
          !---------------
          use = use_0

          !--------------------
          ! define head section
          !--------------------
          hd            = head
          hd% time      = stime   (is)
          hd% db_time   = db_time (is)
          hd% source    = ifile
          hd% record    = is
          hd% id        = entry1
          hd% center    = s1cent  (is)
          hd% subcenter = s1cents (is)
          if ( lpr_std .and. is < npr_extd ) then
             write (6,'()')
             write (6,'( 8(a16, i8,/),   &
                  & 2(a16, a ,/),   &
                  & 6(a16, i8,/) )' )                      &
                  'pe='         ,dace% pe,                         &
                  'head is='    ,is,                               &
                  'obstype='    , hd% obstype  ,                   &
                  'dbkz='       , hd% dbkz     ,                   &
                  'modtype='    , hd% modtype  ,                   &
                  'buf_type='   , hd% buf_type ,                   &
                  'buf_subtype=', hd% buf_subtype,                 &
                  'codetype='   , hd% codetype ,                   &
                  'time='       , cyyyymmddhhmmss (hd% time)   ,   &
                  'db_time='    , cyyyymmddhhmmss (hd% db_time),   &
                  'dbk='        , hd% idbk     ,                   &
                  'source='     , hd% source   ,                   &
                  'record='     , hd% record   ,                   &
                  'id='         , hd% id       ,                   &
                  'center='     , hd% center   ,                   &
                  'subcenter='  , hd% subcenter
          endif

          !--------------------------------------------
          ! perform simple generic check on report type
          !--------------------------------------------
          call check_report_0 (use, hd, 1)
          if (use% state <= STAT_DISMISS) cycle

          !------------------
          ! create new report
          !------------------
          spt0             = empty
          spt0% use        = use
          spt0% hd         = hd
          spti             = spt0

          !----------------------------------
          ! check for valid GNSS station name
          !----------------------------------
          if (ystidn(is) == '') then
             ! Station name empty, reject report
             call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
                  comment='read_netcdf_307024: empty station name')
          end if

          if (verify(trim(ystidn(is)), GNSSset) > 0) then
             ! Station name with invalid characters (not in GNSSset)
             ! => reject station
             call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
                  comment='read_netcdf_307024: station name with invalid character')
          end if

          if (lpr_std) then
             if (spti%use%state <= STAT_DISMISS) then
                write(*,*) '307024: ZTD rejected because of ',     &
                                    'invalid station name : >>>',  &
                                    trim(ystidn(is)), '<<<'
             end if
          end if

          !--------------------------
          ! check for sufficient data
          !--------------------------
!!$          if (mhp   (is) == -999._sp .or. &
!!$              mptid (is) == -1       .or. &
!!$              mda   (is) == -999._sp .or. &
!!$              mde   (is) == -999._sp .or. &
!!$              npdna (is) == -999._sp .or. &
!!$              mlah  (is) == rmissing .or. &
!!$              mloh  (is) == rmissing      ) then
!!$             call decr_rpt_use (spti, CHK_INSDAT, comment=&
!!$                     'read_std_netcdf: insufficient data')
!!$
!!$             if ( lpr_std ) then
!!$                print *,'pe=',dace% pe,   &
!!$                     ' read_std_netcdf after decr_rpt_use CHK_INSDAT is=',is
!!$             endif
!!$             cycle
!!$          endif

          !------------------------------
          ! process report header entries
          !------------------------------
          spti% corme        = max ( s1updat(is), 0)
          spti% col% c% dlat = mlah     (is)
          spti% col% c% dlon = mloh     (is)
          spti% z            = mhp      (is)
          spti% actual_time  = obs_time (is)
          spti% ident        = istidn   (is)
          spti% col% nlev    = 1
          call set_xuv (spti)

          !------------------------
          ! set center / processing
          !------------------------
          call name2spec (                            &
                 name_center   = ystidn(is),          & !<= station name
                 wmo_center    = spti% hd% center,    & !<= BUFR WMO center
                 wmo_subcenter = spti% hd% subcenter, & !<= BUFR WMO subcenter
                 cenpro        = ypcid(is),           & !<= center/product
                 lat           = spti% col% c% dlat,  & !<= station latitude
                 lon           = spti% col% c% dlon,  & !<= station longitude
                 station       = station,             & !=> station ID
                 processing    = spti% stret,         & !=> product center ID
                 center        = spti% center_id     )  !=> processing center ID

          if (spti% stret < 0) then
             ! Station could not be associated with any product => reject
             call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
                  comment='read_netcdf_307024: station without product info')
          end if

          if (lpr_std) then
             if (spti%use%state <= STAT_DISMISS) then
                write(*,*) '307024: ZTD rejected because of missing product',  &
                       'station name: >>>', trim(ystidn(is)),                  &
                       '<<<, product ID = ', spti% stret,                      &
                       ', processing center ID = ', spti% center_id
             end if
          end if

          cnt_prd      = ystidn (is) (6:9)
          if (index(ystidn(is), '-') > 0) then
             ! E-GVAP station name with center/product: remove '-'
             ystidn(is) = ystidn(is)(:index(ystidn(is), '-')-1) // &
                          ystidn(is)(index(ystidn(is), '-')+1:)
          end if
          spti% statid = ystidn (is)          ! currently keep old station name
          !     spti% statid = station              ! possibly new station name

          if ( lpr_std  .and. is < npr_extd ) then
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
                   'spti% col% nlev    = ', spti% col% nlev
          endif

          !----------------
          ! standard checks
          !----------------
          call check_report_1 (spti)
          lk = spti% use% state > STAT_DISMISS
          if (lk) then
             !----------------------------------------------
             ! set slants meta data (currently for ZTD only)
             !----------------------------------------------
             slants% Ntot  = 1 ! number of supporting points along the slant path
             slants% Nmod  = 1 ! number of supporting points below top model layer
             slants% Nup   = 0 ! number of supporting points above top model layer
             slants% Nhor  = 1 ! number of supporting points inside the model grid
             slants% Naloc = 1 ! number of points used (either Nmod or 0 if Nhor<Nmod)
             slants% obs% time      = mjd(obs_time(is))! Modified Julian date of data
             slants% obs% site      = 0                ! GPS station ID
             slants% obs% satellite = 0                ! site ID, last 2 digits of STA
             !slants% obs% elevation = mde  (is) *d2r   ! elevation in radian
             !slants% obs% azimuth   = mda  (is) *d2r   ! azimuth   in radian
             !slants% obs% slant     = nades(is) ! slant delay in m
             !slants% obs% zslant    = nades(is) ! slant delay mapped into zenith direction

             slants% GeoidCorr       = invalid        ! geoid undulation not available

             !call check_store_std (slants, spti, obs, nqfgd(is), lk, VN_ZPD, cnt_prd, station)
             if (lk) nkeep = nkeep + 1
          endif

       end do

    end if  ! End of step 1 - read ZTDs

    !------------------
    ! Step 2: Read STDs
    !------------------
    if (dim_loop_3 == 1) then

       !------------------------
       ! get dimension of fields
       !------------------------
       ! Slant related data with three dimensions:
       ! NPDNA(BUFR_records, Loop_003_maxlen, Loop_004_maxlen)
       ! => BUFR_records        - number of records in file
       ! => Loop_003_maxlen = 1 - STD data available
       ! => Loop_004_maxlen     - max. number of STDs per record

       len_report  = dim_records
       len_slants  = dim_loop_4

       !----------------------------------------
       ! define number of reports in source-file
       !----------------------------------------
       i_source = len_report

       !----------------
       ! allocate arrays
       !----------------
       n = min (len_report, rept_use(OT_GPSGB)% max_proc)

       allocate( ypcid  (n)            )  ! GNSS product string
       allocate( ypcidarr (dim_YPCID,n))  ! GNSS product string for netCDF read
       allocate( mhp    (n)            )  ! Station height
       allocate( ngeun  (n)            )  ! Geoid undulation
       allocate( mptid  (len_slants,n) )  ! Satellite PRN
       allocate( msacl  (len_slants,n) )  ! Satellite system specification
       allocate( mda    (len_slants,n) )  ! Azimuth
       allocate( mde    (len_slants,n) )  ! Elevation
       allocate( npdna  (len_slants,n) )  ! STD
       allocate( nzpdna (len_slants,n) )  ! ZTD = mapped STD

       allocate( UseSTD (len_slants)   )  ! Use current slant: Quality check
       !allocate( (n,len_slants) )

       !--------------------------------
       ! read variables from NetCDF file
       !--------------------------------

       ! Character array GNSS product: YPCID
       ypcid = ''
       status = nf90_inq_varid (ncid, 'YPCID', varid)
       if (status == NF90_NOERR) then
          status = nf90_get_var (ncid, varid, ypcidarr,                  &
                                 count=(/dim_YPCID,n/), start=(/1,rec1/))
          if (status == NF90_NOERR) then
             do is = 1, n
                do j=1, dim_YPCID
                   if (iachar(ypcidarr(j,n)) == iascii_ymissing) then
                      ypcid(is)(j:) = ' '
                      exit
                   end if
                   ypcid(is)(j:j) = ypcidarr(j,is)
                end do
             end do
          end if
       end if
       if (status /= NF90_NOERR) then
          ! GNSS product info not available
          call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS,                   &
                    comment='read_netcdf_307024: No product info (YPCID) available.')
          return
       end if

       ! Real and integer arrays
       call get_real (mhp,    'MHP'  ,   invalsp, (/  n/), (/  rec1/))
       call get_real (ngeun,  'NGEUN',   invalsp, (/  n/), (/  rec1/))
       call get_int  (msacl,  'MSACL',   -1      , (/len_slants,1,n/), (/1,1,rec1/))
       call get_int  (mptid,  'MPTID',   -1      , (/len_slants,1,n/), (/1,1,rec1/))
       call get_real (mda,    'MDA'  ,   invalsp, (/len_slants,1,n/), (/1,1,rec1/))
       call get_real (mde,    'MDE'  ,   invalsp, (/len_slants,1,n/), (/1,1,rec1/))
       call get_real (npdna,  'NPDNA',   invalsp, (/len_slants,1,n/), (/1,1,rec1/))
       call get_real (nzpdna, 'NZPDNA0', invalsp, (/len_slants,1,n/), (/1,1,rec1/))

       if (.false.) then
       !if (.true.) then
          write(*,*)
          write(*,*) 'n, len_slants, rec1 = ', n, len_slants, rec1

          do is = 1, n
             write(*,*) 'Record, is    = ', is
             write(*,*) 'ystidn(is)    = ', ystidn(is)
             write(*,*) 'stime   (is)  = ', stime   (is), cyyyymmddhhmmss (stime(is))
             write(*,*) 'obs_time (is) = ', obs_time (is), cyyyymmddhhmmss (obs_time(is))
             write(*,*) 'mlah  (is)    = ', mlah  (is), rmissing
             write(*,*) 'mloh  (is)    = ', mloh  (is), rmissing
             write(*,*) 'mhp   (is)    = ', mhp   (is)

             do k=1, len_slants
                write(*,*) 'Slant No,  k            = ', k
                write(*,*) 'Elevation: mde  (is,k)  = ', mde  (k,is)   ! elevation
                write(*,*) 'Azimuth:   mda  (is,k)  = ', mda  (k,is)   ! azimuth
                write(*,*) 'STD:       npdna(is,k)  = ', npdna(k,is)   ! slant delay in m
                write(*,*) 'ZTD:       nzpdna(is,k) = ', nzpdna(k,is)  ! STD mapped to zenith
             end do
             write(*,*)
          end do

          call finish('read_netcdf_307024','GPSGB: Debug output: Slant records')
       end if

       !-------------------------------
       ! preset total number of reports
       !-------------------------------
       entry = sum (source(1:ifile-1)% entries) + rec1 - 1

       !------------------
       ! loop over reports
       !------------------
       do is = 1, n

          entry1  = entry    + 1
          entry   = entry    + 1

          !---------------
          ! initialize use
          !---------------
          use = use_0

          !--------------------
          ! define head section
          !--------------------
          hd            = head
          hd% time      = stime   (is)
          hd% db_time   = db_time (is)
          hd% source    = ifile
          hd% record    = is
          hd% id        = entry1
          hd% center    = s1cent  (is)
          hd% subcenter = s1cents (is)

          if ( lpr_std .and. is < npr_extd ) then
             write (6,'()')
             write (6,'( 8(a16, i8,/), &
                  & 2(a16, a ,/), &
                  & 6(a16, i8,/) )' )                              &
                  'pe='         ,dace% pe,                         &
                  'head is='    ,is,                               &
                  'obstype='    , hd% obstype  ,                   &
                  'dbkz='       , hd% dbkz     ,                   &
                  'modtype='    , hd% modtype  ,                   &
                  'buf_type='   , hd% buf_type ,                   &
                  'buf_subtype=', hd% buf_subtype,                 &
                  'codetype='   , hd% codetype ,                   &
                  'time='       , cyyyymmddhhmmss (hd% time)   ,   &
                  'db_time='    , cyyyymmddhhmmss (hd% db_time),   &
                  'dbk='        , hd% idbk     ,                   &
                  'source='     , hd% source   ,                   &
                  'record='     , hd% record   ,                   &
                  'id='         , hd% id       ,                   &
                  'BUFRcenter=' , hd% center   ,                   &
                  'BUFRsubcenter=', hd% subcenter
          endif

          !--------------------------------------------
          ! perform simple generic check on report type
          !--------------------------------------------
          call check_report_0 (use, hd, 1)
          !if (use% state <= STAT_DISMISS) write (6,*) is, 'check_report_0: DISMISS', use%state
          if (use% state <= STAT_DISMISS) cycle

          !------------------
          ! create new report
          !------------------
          spt0             = empty
          spt0% use        = use
          spt0% hd         = hd
          spti             = spt0

          !----------------------------------
          ! check for valid GNSS station name
          !----------------------------------
          if (ystidn(is) == '') then
             ! Station name empty, reject report
             call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS,            &
                          comment='read_netcdf_307024: empty station name')
          end if

          if (verify(trim(ystidn(is)), GNSSset) > 0) then
             ! Station name with invalid characters (not in GNSSset)
             ! => reject station
             call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
                  comment='read_netcdf_307024: station name with invalid character')
          end if

          if (lpr_std) then
             if (spti%use%state <= STAT_DISMISS) then
                write(*,*) '307024: STD rejected because of ',     &
                                    'invalid station name : >>>',  &
                                    trim(ystidn(is)), '<<<'
             end if
          end if

          !--------------------------
          ! check for sufficient data
          !--------------------------
          ! => No need to check GNSS product (ypcid), already done after reading
          ! => Geoid undulation (ngeun) might be missing, can be obtained
          !    using geoid implemented in dace
          if (mhp   (is)         == invalid                 .or. &
              mlah  (is)         == rmissing                .or. &
              mloh  (is)         == rmissing                .or. &
              count(mptid (:,is) == -1)       == len_slants .or. &
              count(msacl (:,is) == -1)       == len_slants .or. &
              count(mda   (:,is) == invalid)  == len_slants .or. &
              count(mde   (:,is) == invalid)  == len_slants .or. &
              count(npdna (:,is) == invalid)  == len_slants      ) then

             call decr_rpt_use (spti, CHK_INSDAT, comment=&
                     'read_std_netcdf: insufficient data')

             if ( lpr_std ) then
                print *,'pe=',dace% pe,   &
                     ' read_std_netcdf after decr_rpt_use CHK_INSDAT is=',is
             endif

             cycle
          endif

          !------------------------------
          ! process report header entries
          !------------------------------
          spti% corme        = max ( s1updat(is), 0)
          spti% col% c% dlat = mlah     (is)
          spti% col% c% dlon = mloh     (is)
          spti% z            = mhp      (is)
          spti% actual_time  = obs_time (is)
          spti% ident        = istidn   (is)
          spti% col% nlev    = 1
          call set_xuv (spti)

          !------------------------
          ! set center / processing
          !------------------------

          ! Keep old code until E-GVAP decides about new BUFR
          !
          !call name2spec (ystidn(is),                              &
          !                spti% hd% center,   spti% hd% subcenter, &
          !                spti% col% c% dlat, spti% col% c% dlon,  &
          !                station,            spti% stret,         &
          !                center                                   )
          !cnt_prd      = ypcid(is)  !ystidn (is) (6:9)
          !if (index(ystidn(is), '-') > 0) then
          !   ! E-GVAP station name with center/product: remove '-'
          !   ystidn(is) = ystidn(is)(:index(ystidn(is), '-')-1) // &
          !                ystidn(is)(index(ystidn(is), '-')+1:)
          !end if

          call name2spec (                            &
                 name_center   = ystidn(is),          & !<= station name
                 wmo_center    = spti% hd% center,    & !<= BUFR WMO center
                 wmo_subcenter = spti% hd% subcenter, & !<= BUFR WMO subcenter
                 cenpro        = ypcid(is),           & !<= center/product
                 lat           = spti% col% c% dlat,  & !<= station latitude
                 lon           = spti% col% c% dlon,  & !<= station longitude
                 station       = station,             & !=> station name
                 processing    = spti% stret,         & !=> DWD product ID
                 center        = spti% center_id     )  !=> DWD processing center ID

          if (spti% stret < 0) then
             ! Station could not be associated with any product => reject
             call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
                    comment='read_netcdf_307024: station without product info')
          end if
          if (lpr_std) then
             if (spti%use%state <= STAT_DISMISS) then
                write(*,*) '307024: STD rejected because of missing product',  &
                       'station name: >>>', trim(ystidn(is)),                  &
                       '<<<, product ID = ', spti% stret,                      &
                       ', processing center ID = ', spti% center_id
             end if
          end if

          spti% statid = ystidn (is)          ! GNSS station name
          cnt_prd      = ypcid(is)            ! GNSS product

          if ( lpr_std  .and. is < npr_extd ) then
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
                   'spti% col% nlev    = ', spti% col% nlev
          endif

          !----------------
          ! standard checks
          !----------------
          call check_report_1 (spti)
          lk = spti% use% state > STAT_DISMISS

          if (lk) then
             ! Find number of valid slants in record

             ! Check each single slant for invalid values
             ! => Only the values used by the operator need to be checked
             UseSTD = 1
             where (mde(:,is)   == invalid) UseSTD = -1
             where (mda(:,is)   == invalid) UseSTD = -1
             where (npdna(:,is) == invalid) UseSTD = -1
             Nstd = count(UseSTD == 1)
             allocate( slants(Nstd) )

             !----------------------------------------------
             ! set slants meta data for STDs
             !----------------------------------------------

             ! Do not change the settings below!
             ! These settings make no sense for STDs but allow to assimilate
             ! STDs with the ZTD operator. Only some near-zenith STDs will
             ! be assimilated but the assimilation runs without errors and
             ! writes STDs to the feedback files. => Nice for testing ...
             ! In case of real STD assimilation the settings have no effect
             ! and will automatically be replaced by meanigful values.
             slants(:)% Ntot  = 1 ! number of supporting points along the slant path
             slants(:)% Nmod  = 1 ! number of supporting points below top model layer
             slants(:)% Nup   = 0 ! number of supporting points above top model layer
             slants(:)% Nhor  = 1 ! number of supporting points inside the model grid
             slants(:)% Naloc = 1 ! number of points used (either Nmod or 0 if Nhor<Nmod)

             ! Satellite PRN
             ! Combine satellite system specification and satellite PRN
             ! BUFR Satellite system specification, descriptor 0 02 020:
             ! GPS      401   =>   0 < NewPRN < 100, New PRN = PRN = mptid(k,is)
             ! GLONASS  402   => 100 < NewPRN < 200, New PRN = 100 + PRN
             ! Galileo  403   => 200 < NewPRN < 300, New PRN = 200 + PRN
             ! BeiDou   404   => 300 < NewPRN < 400, New PRN = 300 + PRN
             ! => NewPRN = (msacl(k,is)-401)*100 + mptid(k,is)

             j = 0
             do k=1, len_slants
                if (UseSTD(k) < 0) cycle
                j = j + 1
                ! Slant observation
                slants(j)% obs% time      = mjd(obs_time(is)) ! MJD of slant
                slants(j)% obs% station   = 0                 ! GPS station ID, not used
                slants(j)% obs% site      = 0                 ! GPS site ID, not used
                slants(j)% obs% satellite =                  &! satellite PRN
                           (msacl(k,is)-401)*100+mptid(k,is)  !   and GNSS
                slants(j)% obs% elevation = mde  (k,is) *d2r  ! elevation in radian
                slants(j)% obs% azimuth   = mda  (k,is) *d2r  ! azimuth   in radian
                slants(j)% obs% slant     = npdna(k,is)       ! slant delay in m
                slants(j)% obs% zslant    = nzpdna(k,is)      ! STD mapped to zenith

                ! GNSS station coordinates
                ! Only the geoid undulation is stored here, the setup of the coordinates
                ! and any coordinate transformations will be done later => spot2station
                slants(j)% GeoidCorr      = ngeun(is)         ! geoid undulation
             end do
             ! Anstelle der 8 das neue Quality flag übergeben ???
             ! => geht nach obs % body % pcc  ????
             ! => deallocates slants
             call check_store_std (slants, spti, obs, 8, lk, VN_SPD, &
                                   cnt_prd, spti% statid)
          endif

       end do    ! do is = 1, n,  loop over reports

    end if  ! End of step 2 - read STDs

  end subroutine read_netcdf_307024

!==============================================================================
  subroutine read_netcdf_307022 (ifile, i_source, obs, head, rec1, recl, &
                                 lkeep, nkeep,cc)
  !========================================================
  ! Read ZTD observations from netCDF (converted BUFR data)
  ! BUFR data are encoded using table D sequence 307022
  ! => E-GVAP ZTD observations
  !========================================================
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data to set
  type (t_head) ,intent(in)           :: head      ! header data already encoded
  integer       ,intent(in)           :: rec1      ! first record to read
  integer       ,intent(in)           :: recl      ! last  record to read
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.
  integer       ,intent(in) ,optional :: cc        ! part of year ccyy


    !================
    ! local variables
    !================
    type (t_use)             :: use             ! status variable
    type (t_head)            :: hd              ! report header
    type (t_spot)            :: spt0, spti      ! report meta data
    type (t_spot) ,save      :: empty           !
    type (SlantData)         :: slants (1:1)    !
    integer                  :: status          ! NetCDF status variable
    integer                  :: ncvars          ! NetCDF number of variables  defined in this NetCDF file
    integer                  :: len_report      ! number of reports in NetCDF file
    integer                  :: len_slants      ! number of reports in NetCDF file
    integer                  :: i, j            ! loop index
    integer                  :: is              ! report loop index
!   integer                  :: nc0             ! dimension for netcdf getvar 1-dimensional arrays
!   integer                  :: nc1             ! first  dimension in start / count
!   integer                  :: nc2             ! second dimension in start / count
    character(NF90_MAX_NAME) :: yname_v         ! NetCDF variable name
    logical                  :: lpr_std         ! std   reports from netcdf are printed
    logical                  :: lpr_extd        ! extended  printing of std
    integer                  :: npr_extd        ! number of extended  printing of temps
    integer                  :: i_met_YSOSN = 0 ! index for station id
    integer                  :: i_met_DATE  = 0 !   date
    integer                  :: i_met_TIME  = 0 !   time
    integer                  :: i_met_MLAH  = 0 !   latitude
    integer                  :: i_met_MLOH  = 0 !   longitude
    integer                  :: i_met_MHP   = 0 !   station height
    integer                  :: i_met_NQFGD = 0 !   quality flags
    integer                  :: i_met_MDA   = 0 !   bearing or azimuth
    integer                  :: i_met_MDE   = 0 !   elevation
    integer                  :: i_met_NADES = 0 !   atmospheric path delay
    integer                  :: i_met_NEERR = 0 !   estimated error in path delay
    integer                  :: numDims         ! number of dimensions for NetCDF variable
    integer                  :: numAtts         ! number of attributes for NetCDF variable
    integer                  :: n               ! number of reports to read
    integer                  :: dimid_slants    !
    integer                  :: dimid_report    !
    integer                  :: varid_NADES     ! variable id
    integer                  :: entry1,entry    ! position in source file (subset)
    integer,  allocatable    :: nqfgd (:)       ! quality flags
    real(sp), allocatable    :: mhp   (:)       ! station height
    real(sp), allocatable    :: mda   (:)       ! bearing or azimuth
    real(sp), allocatable    :: mde   (:)       ! elevation
    real(sp), allocatable    :: nades (:)       ! atmospheric path delay
    real(sp), allocatable    :: neerr (:)       ! estimated error in path delay
    logical                  :: lk              ! flag to keep this report
    character(len=4)         :: cnt_prd         ! code for center/product
    character(len=10)        :: station         ! possibly modified station name
    integer                  :: center          ! center id within DACE (not used so far)
    integer                  :: dimids (dimids_max) ! NetCDF dimension ids

    lpr_std  = .false.
    if (netcdf_verb > 1) lpr_std = .true.

    ! debug
    !lpr_extd = .TRUE.
    !lpr_std  = .TRUE.

    !------------------------------
    ! get default number of reports
    !------------------------------
    len_report =  recl - rec1 + 1
!   nc0        = len_report

    lkeep = .false.
    nkeep = 0

    !---------------------------
    ! get variables ids and name
    !---------------------------
    status = nf90_Inquire (ncid, nVariables=ncvars)
    do j = 1 , ncvars
      status = nf90_Inquire_Variable(ncid, j,name=yname_v )
      !-------------------------------------------------
      ! dimensions:
      !  BUFR_records = UNLIMITED ; // (20670 currently)
      !  Loop_000_maxlen = 25 ;
      !  section1_length = 18 ;
      !  section2_length = 18 ;
      !  YSOSN_strlen = 20 ;
      ! variables to skip:
      !  int edition_number(BUFR_records) ;
      !  byte section1(BUFR_records, section1_length) ;
      !  byte section2(BUFR_records, section2_length) ;
      !  int section....
      !-------------------------------------------------
      if (yname_v(1:7) == 'edition' ) cycle
      if (yname_v(1:7) == 'section' ) cycle
      if (lpr_std .and. lpr_extd )                                 &
      write (6,'(a,i3,a,a16,a)') '3 nf90_Inquire_Variable(i) ',j,  &
                                 ' Variable name(o): ',trim(yname_v)
      !---------------------------------------------------------------------------
      ! char YSOSN (BUFR_records, YSOSN_strlen) ; ""      ; "STATION OR SITE NAME"
      !---------------------------------------------------------------------------
      if (yname_v == 'YSOSN') i_met_YSOSN = j
      !---------------------------------------------------------
      ! int MJJJ (BUFR_records); -2147483647; "YEAR"
      ! int MMM  (BUFR_records); -2147483647; "MONTH"
      ! int MYY  (BUFR_records); -2147483647; "DAY"
      ! int MGG  (BUFR_records); -2147483647; "HOUR"
      ! int NGG  (BUFR_records); -2147483647; "MINUTE"
      ! int DATE (BUFR_records); -2147483647; "Date as YYYYMMDD"
      ! int TIME (BUFR_records); -2147483647; "Time as HHMMSS"
      !---------------------------------------------------------
      if (yname_v == 'DATE') i_met_DATE = j
      if (yname_v == 'TIME') i_met_TIME = j
      !------------------------------------------------------------------------------------------
      ! double MLAH (BUFR_records) ; 9.96920996838687e+36 ; "LATITUDE (HIGH ACCURACY)" ; "DEGREE"
      ! double MLOH (BUFR_records) ; 9.96920996838687e+36 ; "LONGITUDE (HIGH ACCURACY)"; "DEGREE"
      ! int    MHP  (BUFR_records) ; -2147483647          ; "HEIGHT OF STATION" ;"M"
      !------------------------------------------------------------------------------------------
      if (yname_v == 'MLAH') i_met_MLAH = j
      if (yname_v == 'MLOH') i_met_MLOH = j
      if (yname_v == 'MHP' ) i_met_MHP  = j
      !-----------------------------------------------------------------------------------------
      ! int   MTISI (BUFR_records); -2147483647 ; "TIME SIGNIFICANCE"; "CODE_TABLE"     (23) ???
      ! int   NGGTP (BUFR_records); -2147483647 ; "TIME PERIOD OR DISPLACEMENT"; "MINUTE"    ???
      ! float MPPP  (BUFR_records); 9.96921e+36f; "PRESSURE" ; "PA"
      ! float MTN   (BUFR_records); 9.96921e+36f; "TEMPERATURE/DRY BULB TEMPERATURE"
      ! int   MUUU  (BUFR_records); -2147483647 ; "RELATIVE HUMIDITY" ; "%"
      ! int   NQFGD (BUFR_records); -2147483647 ; "QUALITY FLAGS FOR GROUND-BASED GNSS DATA"
      !-----------------------------------------------------------------------------------------
      if (yname_v == 'NQFGD') i_met_NQFGD = j
      !------------------------------------------------------------------------------------------------------
      ! int   MTOTN (BUFR_records)   ; -2147483647 ; "TOTAL NUMBER(ACCUMULATION/AVERAGE)"  (0 or invalid) ???
      ! int   MSACL (BUFR_records, Loop_000_maxlen); "SATELLITE CLASSIFIKATION" ; "CODE_TABLE"
      ! int   MPTID (BUFR_records, Loop_000_maxlen); "PLATFORM TRANSMITTER ID NUMBER"
      ! float MDA   (BUFR_records, Loop_000_maxlen); "BEARING OR AZIMUTH" ; "DEGREE_TRUE"
      ! float MDE   (BUFR_records, Loop_000_maxlen); "ELEVATION" ; "DEGREE"
      ! float NADES (BUFR_records, Loop_000_maxlen); "ATMOSPHERIC PATH DELAY IN SAT. SIGNAL"   ; "M"
      ! float NEERR (BUFR_records, Loop_000_maxlen); "ESTIMATED ERR IN ATMOSPHERIC PATH DELAY" ; "M"
      !------------------------------------------------------------------------------------------------------
      if (yname_v == 'MDA'  ) i_met_MDA   = j
      if (yname_v == 'MDE'  ) i_met_MDE   = j
      if (yname_v == 'NADES') i_met_NADES = j
      if (yname_v == 'NEERR') i_met_NEERR = j
      !-----------------------------------------------------------------------------------------------------------------
      ! int MSSMS(BUFR_records)    ; -2147483647         ; "SAMPLE SCANNING MODE SIGNIFICANCE"              5
      ! float NDPDL(BUFR_records)  ; 9.96921e+36f        ; "DIFF. IN PATH DELAYS F. LIMB VIEWS ...." ; "M"  -0.00999, _,
      ! float NEEPDD(BUFR_records) ; 9.96921e+36f        ; "ESTIMATED ERR IN PATH DELAY DIFFERENCE" ; "M"   0, _,
      ! int MSSMS0(BUFR_records)   ; -2147483647         ; "SAMPLE SCANNING MODE SIGNIFICANCE"              6
      ! float NDPDL0(BUFR_records) ; 9.96921e+36f        ; "DIFF. IN PATH DELAYS F. LIMB VIEWS ...." ; "M"  _, -0.09999
      ! float NEEPDD0(BUFR_records); 9.96921e+36f        ; "ESTIMATED ERR IN PATH DELAY DIFFERENCE" ; "M"   0, _,
      ! float NCZWV(BUFR_records)  ; 9.96921e+36f        ; "COMP OF ZENITH PATH DELAY DUE TO WATER V" ; "M"
      ! float NWLN(BUFR_records)   ; 9.96921e+36f        ; "PRECIPITABLE WATER" ; "KG/M**2"
      ! float MLIED(BUFR_records)  ; 9.96921e+36f        ; "LOG10 OF INTEGRATED ELECTRON DENSITY" ; "LOG(1/M**2)"
      !-----------------------------------------------------------------------------------------------------------------
    enddo
    !--------------------------------------
    ! printout dimension indices to be used
    !--------------------------------------
    if ( lpr_std .and. lpr_extd )                                                  &
    write (6,'(4(a12,i3))') ' i_met_YSOSN ',i_met_YSOSN,' i_met_DATE ',i_met_DATE, &
                            ' i_met_TIME  ',i_met_TIME, ' i_met_MLAH ',i_met_MLAH, &
                            ' i_met_MLOH  ',i_met_MLOH, ' i_met_MHP ' ,i_met_MHP,  &
                            ' i_met_NQFGD ',i_met_NQFGD,' i_met_MDA ' ,i_met_MDA,  &
                            ' i_met_MDE '  ,i_met_MDE,  ' i_met_NADES',i_met_NADES,&
                            ' i_met_NEERR ',i_met_NEERR
    !------------------------
    ! get dimension of fields
    !------------------------
    status = nf90_inq_varid        (ncid, 'NADES' ,  varid_NADES)
    status = nf90_Inquire_Variable (ncid, varid_NADES, ndims=numDims, dimids=dimids, natts=numAtts)

!   nc1 = 0
!   nc2 = 0
    len_slants  = 0
    len_report  = 0

    if(numDims >= 2) then
      dimid_slants = dimids(1)
      dimid_report = dimids(2)
      status = nf90_Inquire_Dimension(ncid, dimid_slants, len=len_slants)
      status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
!     nc1       = len_slants
!     nc2       = len_report
    endif

    !----------------------------------------
    ! define number of reports in source-file
    !----------------------------------------
    i_source = len_report
    !----------------
    ! allocate arrays
    !----------------
    n = min (len_report, rept_use(OT_GPSGB)% max_proc)
!   allocate (ysosn (n))
!   allocate (date  (n))
!   allocate (time  (n))
!   allocate (mlah  (n))
!   allocate (mloh  (n))
    allocate (nqfgd (n))
    allocate (mhp   (n))
    allocate (mda   (n))
    allocate (mde   (n))
    allocate (nades (n))
    allocate (neerr (n))

    !--------------------------------
    ! read variables from NetCDF file
    !--------------------------------
    call get_int  (nqfgd, 'NQFGD', -1      , (/  n/), (/  rec1/))
    call get_real (mhp,   'MHP'  , -999._sp, (/  n/), (/  rec1/))
    call get_real (mda,   'MDA'  , -999._sp, (/1,n/), (/1,rec1/))
    call get_real (mde,   'MDE'  , -999._sp, (/1,n/), (/1,rec1/))
    call get_real (nades, 'NADES', -999._sp, (/1,n/), (/1,rec1/))
    call get_real (neerr, 'NEERR',   -1._sp, (/1,n/), (/1,rec1/))

    !-------------------------------
    ! preset total number of reports
    !-------------------------------
    entry = sum (source(1:ifile-1)% entries)

    !------------------
    ! loop over reports
    !------------------
    do is = 1, n

      entry1  = entry    + 1
      entry   = entry    + 1

      !---------------
      ! initialize use
      !---------------
      use = use_0

      !--------------------
      ! define head section
      !--------------------
      hd            = head
      hd% time      = stime   (is)
      hd% db_time   = db_time (is)
      hd% source    = ifile
      hd% record    = is
      hd% id        = entry1
      hd% center    = s1cent  (is)
      hd% subcenter = s1cents (is)
      if ( lpr_std .and. is < npr_extd ) then
        write (6,'()')
        write (6,'( 8(a16, i8,/),   &
                  & 2(a16, a ,/),   &
                  & 6(a16, i8,/) )' )                      &
          'pe='         ,dace% pe,                         &
          'head is='    ,is,                               &
          'obstype='    , hd% obstype  ,                   &
          'dbkz='       , hd% dbkz     ,                   &
          'modtype='    , hd% modtype  ,                   &
          'buf_type='   , hd% buf_type ,                   &
          'buf_subtype=', hd% buf_subtype,                 &
          'codetype='   , hd% codetype ,                   &
          'time='       , cyyyymmddhhmmss (hd% time)   ,   &
          'db_time='    , cyyyymmddhhmmss (hd% db_time),   &
          'dbk='        , hd% idbk     ,                   &
          'source='     , hd% source   ,                   &
          'record='     , hd% record   ,                   &
          'id='         , hd% id       ,                   &
          'center='     , hd% center   ,                   &
          'subcenter='  , hd% subcenter
      endif

      !--------------------------------------------
      ! perform simple generic check on report type
      !--------------------------------------------
      call check_report_0 (use, hd, 1)
      if (use% state <= STAT_DISMISS) cycle

      !------------------
      ! create new report
      !------------------
      spt0             = empty
      spt0% use        = use
      spt0% hd         = hd
      spti             = spt0

      !----------------------------------
      ! check for valid GNSS station name
      !----------------------------------
      if (ystidn(is) == '') then
         ! Station name empty, reject report
         call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
              comment='read_netcdf_307022: empty station name')
      end if

      if (verify(trim(ystidn(is)), GNSSset) > 0) then
         ! Station name with invalid characters (not in GNSSset)
         ! => reject station
         call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
              comment='read_netcdf_307022: station name with invalid character')
      end if
      
      if (lpr_std) then
         if (spti%use%state <= STAT_DISMISS) then
            write(*,*) '307022: ZTD rejected because of invalid station name : >>>',  &
                        trim(ystidn(is)), '<<<'
         end if
      end if

      !--------------------------
      ! check for sufficient data
      !--------------------------
      if (mhp   (is) == -999._sp .or. &
          mda   (is) == -999._sp .or. &
          mde   (is) == -999._sp .or. &
          nades (is) == -999._sp .or. &
          mlah  (is) == rmissing .or. &
          mloh  (is) == rmissing      ) then
        call decr_rpt_use (spti, CHK_INSDAT, comment='read_netcdf_307022: insufficient data')

        if ( lpr_std ) then
          print *,'pe=',dace% pe,' read_std_netcdf after decr_rpt_use CHK_INSDAT is=',is
        endif
        cycle
      endif

      !------------------------------
      ! process report header entries
      !------------------------------
      spti% corme        = max ( s1updat(is), 0)
      spti% col% c% dlat = mlah     (is)
      spti% col% c% dlon = mloh     (is)
      spti% z            = mhp      (is)
      spti% actual_time  = obs_time (is)
      spti% ident        = istidn   (is)
      spti% col% nlev    = 1
      call set_xuv (spti)
      !------------------------
      ! set center / processing
      !------------------------
!!$      call name2spec (ystidn(is),                              &
!!$                      spti% hd% center,   spti% hd% subcenter, &
!!$                      spti% col% c% dlat, spti% col% c% dlon,  &
!!$                      station,            spti% stret,         &
!!$                      center                                   )

      call name2spec (                            &
             name_center   = ystidn(is),          & !<= E-GVAP name NNNN-PPPP
             wmo_center    = spti% hd% center,    & !<= BUFR WMO center
             wmo_subcenter = spti% hd% subcenter, & !<= BUFR WMO subcenter
             lat           = spti% col% c% dlat,  & !<= station latitude
             lon           = spti% col% c% dlon,  & !<= station longitude
             station       = station,             & !=> station name (no product)
             processing    = spti% stret,         & !=> product ID
             center        = spti% center_id     )  !=> processing center ID

      if (spti% stret < 0) then
         ! Station could not be associated with any product => reject
         call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
              comment='read_netcdf_307022: station without product info')
      end if

      if (lpr_std .and. spti% stret < 0) then
         if (spti%use%state <= STAT_DISMISS) then
            write(*,*) '307022: ZTD rejected because of missing product. ',          &
                       'Station name = >>>', trim(ystidn(is)), '<<<, product ID = ', &
                       spti% stret, ', processing center ID = ',              &
                       spti% center_id
         end if
      end if

      cnt_prd      = ystidn (is) (6:9)
      if (index(ystidn(is), '-') > 0) then
         ! E-GVAP station name with center/product: remove '-'
         ystidn(is) = ystidn(is)(:index(ystidn(is), '-')-1) // &
                      ystidn(is)(index(ystidn(is), '-')+1:)
      end if
      spti% statid = ystidn (is)          ! currently keep old station name
!     spti% statid = station              ! possibly new station name

      if ( lpr_std  .and. is < npr_extd ) then
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
              'spti% col% nlev    = ', spti% col% nlev
      endif

      !-----------------------
      ! check for invalid data
      !-----------------------
      ! Invalid GNSS station coordinates: lon = lat = 0.000 (lon/lat in degrees)
      if ( abs(spti% col% c% dlat) < 1.0E-5_wp .and.       &
           abs(spti% col% c% dlon) < 1.0E-5_wp      ) then
         call decr_rpt_use (spti, CHK_INSDAT, STAT_DISMISS, &
              comment='read_std_netcdf: lon = lat = 0.000')
      end if

      !----------------
      ! standard checks
      !----------------
      call check_report_1 (spti)
      lk = spti% use% state > STAT_DISMISS

      if (lk) then
        !----------------------------------------------
        ! set slants meta data (currently for ZTD only)
        !----------------------------------------------
        slants% Ntot  = 1 ! number of supporting points along the slant path
        slants% Nmod  = 1 ! number of supporting points below top model layer
        slants% Nup   = 0 ! number of supporting points above top model layer
        slants% Nhor  = 1 ! number of supporting points inside the model grid
        slants% Naloc = 1 ! number of points used (either Nmod or 0 if Nhor<Nmod)
        !        slants% assimilate = .true.
        slants% GeoidCorr      = -999._sp         ! geoid undulation not available
        slants% obs% time      = mjd(obs_time(is))! Modified Julian date of data
        slants% obs% site      = 0                ! GPS station ID
        slants% obs% satellite = 0                ! site ID, last 2 digits of STA
        slants% obs% elevation = mde  (is) *d2r   ! elevation in radian
        slants% obs% azimuth   = mda  (is) *d2r   ! azimuth   in radian
        slants% obs% slant     = nades(is) ! slant delay in m
        slants% obs% zslant    = nades(is) ! slant delay mapped into zenith direction

        call check_store_std (slants, spti, obs, nqfgd(is), lk, VN_ZPD, cnt_prd, station)
        if (lk) nkeep = nkeep + 1
      endif
   end do

  end subroutine read_netcdf_307022

! !------------------------------------------
! ! read observations from bufrx2netcdf files
! !------------------------------------------
!     !------------------
!     ! loop over reports
!     !------------------
!     do j = 1, n_report
!       !---------------------------------------------------
!       ! store header data, insert DWD dbkz, perform checks
!       !---------------------------------------------------
!       head% modtype    = GPSGB     ! module type
!       head% obstype    = OT_GPSGB  ! observation type
!       head% center     = WMO0_DWD  ! originating center (DWD)
!       head% subcenter  = 173       !         sub-center (GFZ)
!       head% buf_type   = buf_type
!       head% source     = if        ! source file number
!       head% record     = j         ! record in file
!       head% time       = time_mjd (t_report(j))
!       head% dbkz       = 94
!       use = use_0
!       !-------------
!       ! basic checks
!       !-------------
!       call check_report_0 (use, head, 1)
!       i = f_slant (j)
!       if (.not.associated (STD(i)% station)) then
!         call decr_use (use, STAT_DISMISS, FL_PRACTICE)
!         use% state = min (use% state, STAT_DISMISS)
!         use% check = FL_PRACTICE
!         write (6,*)                                               &
!          'read_std_ascii: no station info present for station # ',&
!           STD(i)% Obs% station
!       endif
!       if (use% state <= STAT_DISMISS) cycle
!
!       !-------------------------------------------
!       ! Final preparation of observation type data
!       ! Storage into components of 'obs'
!       !-------------------------------------------
!       i = f_slant (j)
!       spt = empty
!       spt% use                 = use
!       spt% hd                  = head
!       spt% hd% satid           = STD(i)% obs%     satellite
!       spt% ident               = STD(i)% station% id
!       spt% statid              = STD(i)% station% name
!       spt% col% c% dlat        = STD(i)% station% CoordGeo(2) * r2d
!       spt% col% c% dlon        = STD(i)% station% CoordGeo(1) * r2d
!       spt% col%    nlev        = 1
!       call set_xuv (spt)
!       call check_store_std (STD(i:i+n_slant(j)-1), spt, obs, lkeep)
!     end do
! end subroutine read_std_grib
!------------------------------------------------------------------------------
  subroutine read_std_ascii (obs)
  type (t_obs)     ,intent(inout)         :: obs    ! observations data type
  !----------------------------------------------------
  ! read Slant Totel Delay observations from ASCII file
  ! and store in 3dvar/LETKF observation data types
  !----------------------------------------------------
    integer              :: i, j        ! index variable
    integer              :: if          ! file index variable
    character(len=256)   :: fname       ! file name
    integer              :: n_report
    integer ,allocatable :: n_slant (:)
    integer ,allocatable :: f_slant (:)
    real(dp),allocatable :: t_report(:)
    integer              :: n_total

    type(t_spot)         :: spt       ! report meta data variable
    type(t_spot)         :: empty     ! default initialised
    type(t_head)         :: head      ! data usually stored in BUFR header
    type(t_use)          :: use       ! state of the report
!   type(t_time)         :: time
    logical              :: lkeep

    !------------------------------
    ! loop over list of input files
    !------------------------------
    do if=1, n_source
      !-------------------------------------------------
      ! check for observation type, file type, processor
      !-------------------------------------------------
      if (source(if)% obstype  /= OT_GPSGB)    cycle
      if (source(if)% filetype /= FT_SPECIFIC) cycle
      if (source(if)% pe       /= dace% pe)    cycle
      !--------------
      ! read the data
      !--------------
      fname = path_file (source(if)% path, source(if)% file)
      call read_std_obs (fname)  !  (obs% o)

!-------------------
! check what we have
!-------------------
!print *,&
!'i,stat,mjd,mjd(time),time,slant,ass(station),Name,ID,CoordGeo,ass(Steps,Cart,Ell,Geo,Idx),Ntot,Nmod,Nup,assimilate'
!do i = 1, size (STD)
!if(    i<11.&
!   .or.i>size (STD)-10&
!   .or.(i>32020 .and. i<32050))then
!time = time_mjd (STD(i)% Obs% time)
!
!if (associated (STD(i)% station)) then
!  print *, i,&
!  STD(i)% Obs% station,&
!  STD(i)% Obs% time,&
!  mjd(time),&
!  cyyyymmddhhmmss(time),&
!  STD(i)% Obs% slant,&
!  associated (STD(i)% station),&
!  STD(i)% station% Name,&
!  STD(i)% station% ID,&
!  STD(i)% station% CoordGeo,&
!  associated (STD(i)% SlantSteps),&
!  associated (STD(i)% SlantPointsCart),&
!  associated (STD(i)% SlantPointsEll),&
!  associated (STD(i)% SlantPointsGeo),&
!  associated (STD(i)% SlantPointsIdx),&
!  STD(i)% Ntot,&
!  STD(i)% Nmod,&
!  STD(i)% Nup,&
!  STD(i)% assimilate
!else
!  print *, i,&
!  STD(i)% Obs% station,&
!  STD(i)% Obs% time,&
!  mjd(time),&
!  cyyyymmddhhmmss(time),&
!  STD(i)% Obs% slant,&
!  associated (STD(i)% station)
!endif
!endif
!end do

      if (associated(STD)) then

         ! Was tun, wenn "read_std_obs" keine Beobachtungen findet?
         ! size (STD) = 1 auch wenn der Pointer = Null ist.
         ! Ganz ohne Daten testen ...

         !------------------------------------------------------------------
         ! Store geoid undulation in STD(i)%GeoidCorr
         ! The station coordinates in STD(i)%Station% ... will not be copied
         ! and get lost.
         !------------------------------------------------------------------
         do i = 1, size (STD)
            STD(i)%GeoidCorr = STD(i)%Station%GeoidCorr
         end do

      !------------------------------------------------------------------
      ! regard all slants related to a station and a time as one report
      ! we assume that slants are sorted according to station id and time
      !------------------------------------------------------------------
      n_report = 1
      do i = 2, size (STD)
        if (STD(i)% Obs% station /= STD(i-1)% Obs% station .or.   &
            STD(i)% Obs% time    /= STD(i-1)% Obs% time         ) &
           n_report = n_report + 1
      end do
      allocate (n_slant  (n_report))
      allocate (f_slant  (n_report))
      allocate (t_report (n_report))
      n_slant     = 0
      n_slant (1) = 1
      f_slant (1) = 1
      t_report(1) = STD(1)% Obs% time
      n_total     = 1
      j = 1
      do i = 2, size (STD)
        if (STD(i)% Obs% station /= STD(i-1)% Obs% station .or. &
            STD(i)% Obs% time    /= STD(i-1)% Obs% time         ) then
          j = j + 1
          f_slant (j) = i
          t_report(j) = STD(i)% Obs% time
        endif
        n_slant(j) = n_slant(j) + 1
        n_total    = n_total    + 1
      end do
      !------------------
      ! loop over reports
      !------------------
      do j = 1, n_report
        !---------------------------------------------------
        ! store header data, insert DWD dbkz, perform checks
        !---------------------------------------------------
        head% modtype    = GPSGB     ! module type
        head% obstype    = OT_GPSGB  ! observation type: GNSS ground-based
        head% codetype   = OC_GPSGB  ! GNSS STDs
        head% center     = WMO0_DWD  ! originating center (DWD)
        head% subcenter  = 173       !         sub-center (GFZ)
!       head% buf_type   = buf_type
        head% source     = if        ! source file number
        head% record     = j         ! record in file
        head% time       = time_mjd (t_report(j))
!       head% dbkz       = 94
        use = use_0
        !-------------
        ! basic checks
        !-------------
        call check_report_0 (use, head, 1)
        i = f_slant (j)
        if (.not.associated (STD(i)% station)) then
          call decr_use (use, STAT_DISMISS, FL_PRACTICE)
!         use% state = min (use% state, STAT_DISMISS)
!         use% check = FL_PRACTICE
          write (6,*)                                               &
           'read_std_ascii: no station info present for station # ',&
            STD(i)% Obs% station
        endif
        if (use% state <= STAT_DISMISS) cycle

        !-------------------------------------------
        ! Final preparation of observation type data
        ! Storage into components of 'obs'
        !-------------------------------------------
        i = f_slant (j)
        spt = empty
        spt% use                 = use
        spt% hd                  = head
        spt% hd% satid           = STD(i)% obs%     satellite
        spt% ident               = STD(i)% station% id
        spt% statid              = STD(i)% station% name
        spt% col% c% dlat        = STD(i)% station% CoordGeo(2) * r2d
        spt% col% c% dlon        = STD(i)% station% CoordGeo(1) * r2d
        spt% z                   = STD(i)% station% CoordGeo(3)
        spt% col%    nlev        = 1
        spt% actual_time         = time_mjd (t_report(j))
        call set_xuv (spt)
        !-------------------------------------------------------------------
        !  BUFR table 0 33 038 :
        !  Quality flags for ground-based GNSS* data
        !  Bit No.                                         (fortran ibset)
        !  1 Total zenith delay quality is considered poor  0              1
        !  2 GALILEO satellites used                        1              2
        !  3 GLONASS satellites used                        2              4
        !  4 GPS satellites used                            3              8
        !  5 Meteorological data applied                    4
        !  6 Atmospheric loading correction applied         5
        !  7 Ocean tide loading applied                     6
        !  8 Climate quality data processing                7
        !  9 Near-real time data processing                 8
        !  All 10 Missing value
        !-------------------------------------------------------------------
        ! pcc = 8                                 ! GPS used
        ! pcc = ibset (obs% body (i1:in)% pcc, 3) ! GPS used
        !-------------------------------------------------------------------
        call check_store_std (STD(i:i+n_slant(j)-1), spt, obs, 8, lkeep, VN_SPD, '', '')
      end do
      deallocate (STD)
      deallocate (n_slant)
      deallocate (f_slant)
      deallocate (t_report)
      end if    ! if (allocated(STD)) then
    end do
  end subroutine read_std_ascii
!------------------------------------------------------------------------------
  subroutine check_store_std (slants, spot, obs, pcc, lkeep, lvarno, &
                              cnt_prd, station)
  type(SlantData) ,intent(inout)  :: slants (:)
  type(t_spot)    ,intent(inout)  :: spot
  type(t_obs)     ,intent(inout)  :: obs
  integer         ,intent(in)     :: pcc
  logical         ,intent(out)    :: lkeep   ! keep or reject observation
  integer         ,intent(in)     :: lvarno  ! set varno to VN_ZPD or VN_SPD
  character(len=*),intent(in)     :: cnt_prd ! center/product code
  character(len=*),intent(in)     :: station ! station name (w/o product code)

    type(t_spot),pointer :: spt          ! temporary
    type(t_datum)        :: bod          ! body derived type
    integer              :: n_slants     ! number of slants in report
    integer              :: id           ! observation id
    integer              :: i1, in       ! index range
    integer              :: status       ! active passive etc

    lkeep = .true.

    !----------------------
    ! evaluate quality flag
    !----------------------
!   call decr_rpt_use (spot, CHK_NOTUSED, STAT_PASSIVE, &
!                                         comment='....')

    !------------------------------
    ! Select slants with valid data
    ! Check consistency of metadata
    !------------------------------
    n_slants = size (slants)

    !------------------------------
    ! exit if no valid data present
    !------------------------------
!   if (nray == 0) then
!     lkeep = .false.
!     call decr_rpt_use (spot, CHK_INSDAT, comment='nray <= 0')
!     return
!   endif

    !-------------------------------
    ! sort slants if necessary
    ! (time)
    !-------------------------------

    !---------------------
    ! metadata consistency
    !---------------------

    !--------------------------
    ! modify observation errors
    !--------------------------

    !----------------------
    ! pick up slants
    !----------------------

    !---------
    ! printout
    !---------

    !====================================
    ! store data in observation data type
    !====================================
!   if (occ_int_size == 0) call set_size
    !----------
    ! meta data
    !----------
    spot%         int_type    = ITY_MCOLS
    spot%         cost        = n_slants * 10._wp
    spot%         nr          = n_slants           ! diagonal R so far
    spot%         char        = CHR_NONL+CHR_EXP
    !spot%         actual_time = spot% hd% time

    !--------------------------
    ! report selection (filter)
    !--------------------------
    call check_report_1 (spot)
!   if (sat_ids(1) /= -1 .and. all(sat_ids /= spot% hd% satid)) &
!     call decr_rpt_use (spot, CHK_BLACKLIST)

    call std_status (int (spot% hd% center),    &! <-  data provider
                     int (spot% hd% subcenter), &! <-  data provider
                          lvarno,               &! <-  ZTD/STD
                          spot% stret,          &! <-  product
                          cnt_prd,              &! <-  center/product
                          station,              &! <-  station name
                          status,               &!  -> passive, accepted etc
                          spot% pcc             )!  -> preference
    if (status > -1) then
      call reverse_code (status, stats)          ! feedback -> DACE convention
      call decr_rpt_use (spot, CHK_NONE, status, 'namelist filter')
    else
      call decr_rpt_use (spot, CHK_WHITELIST, STAT_DISMISS, 'not whitelisted')
    endif
    if (spot% use% state > STAT_DISMISS) then
      !--------------------------------------------------
      ! new report header entry in DACE  observation list
      !--------------------------------------------------
      call new_spot (obs,1, set_id=.true.)
      spt => obs% spot (obs% n_spot)
      id      = spt% id
      spt     = spot
      spt% id = id
      !------------------------------------------------
      ! new report body entry in DACE  observation list
      !------------------------------------------------
      if (lvarno == VN_SPD) then
         bod % mn          = 'std'
      else if (lvarno == VN_ZPD) then
         bod % mn          = 'ztd'
      end if
      bod % use % state = STAT_ACTIVE
      bod % use % check = CHK_NONE
      call new_obs (obs, n_slants, spot=spt)
      spt% nr                                                = n_slants
      i1 = spt% o% i+1
      in = spt% o% i + spt% o% n
      obs % varno (i1 : in)             = lvarno  ! VN_ZPD or VN_SPD
      obs %  olev (i1 : in)             = slants% obs% elevation * r2d
      obs %  body (i1 : in)             = bod
      obs %  body (i1 : in)% o          = slants% obs% slant
      obs %  body (i1 : in)% eo         = ZTDerror  ! from STD_OBS namelist
      obs %  body (i1 : in)% obs_par(1) = slants% obs% azimuth   * r2d
      obs %  body (i1 : in)% lev_sig    = slants% obs% satellite
      obs %  body (i1 : in)% lev_typ    = VN_ELEV
      obs %  body (i1 : in)% pcc        = pcc

      !------------------------------------------------
      ! store information read in observation data type
      !------------------------------------------------

!do i=1,spt% o% n
!write(6,*)'check_store_std ',spt% id, spt% statid, slants(i)% obs% time, slants(i)% obs% satellite
!end do

      call store_std (obs, spt, slants)
      call destruct            (slants)
    else
      lkeep = .false.
    endif

  end subroutine check_store_std
!------------------------------------------------------------------------------
  subroutine set_size
  !--------------------------------------------------------------
  ! store sizes of derived data type SlantData and its components
  ! (in termes of size of component OBS% PAR) into
  !-------------------------------------------------------------
    type (SlantData)        :: sd
!   type (SlantTotalDelay)  :: std
    type (t_obs)            :: obs
    if (sd_size == 0) then
        sd_size = size (transfer (sd   ,obs% par))
       int_size = size (transfer (1    ,obs% par))
        wp_size = size (transfer (1._wp,obs% par))
    endif
  end subroutine set_size
!------------------------------------------------------------------------------
  subroutine store_std (obs, spot, slants)
  type(t_obs)     ,intent(inout) :: obs        ! data of all observations
  type(t_spot)    ,intent(inout) :: spot       ! meta data of this observation
  type(SlantData) ,intent(in)    :: slants (:) ! slant meta data
  !-----------------------------------------------------------------------
  ! Store the data from variables SLANTS in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------
    integer ,pointer :: par (:)
    integer          :: n1, nn, m, i
    integer          :: ns
    integer          :: np, np3, npp
    if (sd_size == 0) call set_size
    ns  = size (slants)
    m   = sd_size * ns
    np  = 0
    npp = 0
    np3 = 0
    do i = 1, ns
       if (associated (slants(i)% SlantSteps)) then
        np  = slants(i)% Ntot                     * wp_size
        np3 = slants(i)% Ntot  * 3                * wp_size
        npp = slants(i)% Naloc * slants(i)% Nnghb
        m   = m + (np + np3 * 3 + npp) *  wp_size &
                +                 npp  * int_size
      endif
    end do

    if (spot% p% i < 0) call new_par (obs, m, spot=spot)
    if (m > spot% p% n) call new_par (obs, m, spot=spot)
    if (m < spot% p% n) spot% p% n = m
!   if (m > spot% p% n) call finish('store_std','m > spot% p% n')
    par => obs % par (spot% p% i+1 : spot% p% i + spot% p% n)
    n1 = 1; nn = sd_size * ns
    par (n1:nn) = transfer (slants ,par)
    do i = 1, ns
      if (associated (slants(i)% SlantSteps)) then
        np  = slants(i)% Ntot                     * wp_size
        np3 = slants(i)% Ntot  * 3                * wp_size
        npp = slants(i)% Naloc * slants(i)% Nnghb
        n1 = nn + 1; nn = nn + np
        par (n1:nn) = transfer (slants(i)% SlantSteps      ,par)
        n1 = nn + 1; nn = nn + np3
        par (n1:nn) = transfer (slants(i)% SlantPointsCart ,par)
        n1 = nn + 1; nn = nn + np3
        par (n1:nn) = transfer (slants(i)% SlantPointsEll  ,par)
        n1 = nn + 1; nn = nn + np3
        par (n1:nn) = transfer (slants(i)% SlantPointsGeo  ,par)
        n1 = nn + 1; nn = nn + npp * int_size
        par (n1:nn) = transfer (slants(i)% idxi            ,par)
        n1 = nn + 1; nn = nn + npp *  wp_size
        par (n1:nn) = transfer (slants(i)% wi              ,par)
      endif
    end do

  end subroutine store_std
!------------------------------------------------------------------------------
  subroutine load_std (obs, spot, slants)
  type(t_obs)     ,intent(in) :: obs        ! data of all observations
  type(t_spot)    ,intent(in) :: spot       ! meta data of this observation
  type(SlantData) ,pointer    :: slants (:) ! slant meta data
  !-----------------------------------------------------------------
  ! retrieve GPSGB specific table (type(SlantData)) from 3dvar table
  !-----------------------------------------------------------------
    integer         ,pointer :: par (:)
    integer                  :: n1, nn
    integer                  :: ns
    integer                  :: i
    integer                  :: istat
    integer                  :: np, np3, npp, nt, na, ng
    type(SlantData) ,pointer :: s

    !-------------------------------
    ! retrieve containers for slants
    !-------------------------------
    if (sd_size == 0) call set_size
    ns = spot% o% n
    n1 = 1; nn = sd_size * ns
    allocate (slants (ns))
    if (spot% p% n < nn) then
      write(0,*)'load_std: p%n (',spot% p% n,') smaller than',nn
      call finish ('load_std','GPSGB: p%n too small')
    endif
    par => obs % par (spot% p% i+1 : spot% p% i + spot% p% n)
    slants = transfer(par (n1:nn), slants)
    !-----------------------------------------------------------
    ! retrieve components:
    !  SlantPointsCart, SlantPointsEll, SlantPointsGeo, idxi, wi
    !-----------------------------------------------------------
    do i = 1, ns
      s => slants(i)
      if (associated (s% SlantSteps)) then
        na  = s% Naloc
        ng  = s% Nnghb
        nt  = s% Ntot
        np  = s% Ntot             * wp_size
        np3 = s% Ntot  * 3        * wp_size
        npp = s% Naloc * s% Nnghb
        n1 = nn + 1; nn = nn + np
        allocate (s% SlantSteps (nt))
        s% SlantSteps = transfer (par (n1:nn), (/1._wp/))
        n1 = nn + 1; nn = nn + np3
        allocate (s% SlantPointsCart (nt,3))
        s% SlantPointsCart = reshape (transfer (par(n1:nn),(/1._wp/)),(/nt,3/))
        n1 = nn + 1; nn = nn + np3
        allocate (s% SlantPointsEll  (nt,3))
        s% SlantPointsEll  = reshape (transfer (par(n1:nn),(/1._wp/)),(/nt,3/))
        n1 = nn + 1; nn = nn + np3
        allocate (s% SlantPointsGeo  (nt,3))
        s% SlantPointsGeo  = reshape (transfer (par(n1:nn),(/1._wp/)),(/nt,3/))
        n1 = nn + 1; nn = nn + npp * int_size
        allocate (s% idxi            (ng, na))
        s% idxi           = reshape (transfer (par(n1:nn),(/1    /)),(/ng,na/))
        n1 = nn + 1; nn = nn + npp *  wp_size
        allocate (s% wi              (ng, na))
        s% wi             = reshape (transfer (par(n1:nn),(/1._wp/)),(/ng,na/))
      else
        nullify (s% SlantPointsCart)
        nullify (s% SlantPointsEll)
        nullify (s% SlantPointsGeo)
        nullify (s% idxi)
        nullify (s% wi)
      endif
      nullify (s% SlantPointsIdx)
      if (associated (s% Station)) then
        istat = GNSShash(s% Obs% Station)
        if (istat .gt. 0) then
          s% Station => CoordList(istat)
        else
          call finish('load_std','GPSGB: GNSShash <= 0')
        endif
      endif
    end do
  end subroutine load_std
!------------------------------------------------------------------------------
  subroutine read_fdbk_gpsgb (o)
  type (t_obs) ,intent(inout) :: o          ! observation data type variable
  !-----------------------------------------------
  ! restore SlantData from t_spot, body
  ! after observation was read from feedback file
  !-----------------------------------------------
    type (t_spot) ,pointer        :: s         ! pointer to report
    integer                       :: i         ! report      index
    integer                       :: j         ! observation index
    integer                       :: k         ! slant index
    integer                       :: n         ! number of STD/ZTDs in report
    type (SlantData), allocatable :: slants(:) ! slant meta data and observations

    !------------------
    ! loop over reports
    !------------------
    do i = 1, o% n_spot
      if  (o% spot(i)% hd% obstype /= OT_GPSGB) cycle
      s => o% spot(i)
      n =  s% o% n
      allocate (slants(n))
      !-------------------------
      ! GPSGB specific meta data
      !-------------------------
      s% int_type    = ITY_MCOLS
      s% cost        = n * 10._wp
      s% nr          = n * 1            ! diagonal R so far
      s% char        = CHR_NONL+CHR_EXP
      !--------------------------------------------------------------
      ! restore slants meta data from t_spot (currently for ZTD only)
      !--------------------------------------------------------------
      !if (s% o% n /= 1) call finish ('read_fdbk_gpsgb','more than 1 obs in spot')
      do k = 1, n
        j = s% o% i + k

        slants(k)% Ntot  = 0 ! number of supporting points along the slant path
        slants(k)% Nmod  = 0 ! number of supporting points below top model layer
        slants(k)% Nup   = 0 ! number of supporting points above top model layer
        slants(k)% Nhor  = 0 ! number of supporting points inside the model grid
        slants(k)% Naloc = 0 ! number of points used (either Nmod or 0 if Nhor<Nmod)

        if (o% varno (j) == VN_ZPD) then
          ! ???
        else if (o% varno (j) == VN_SPD) then
          ! ???
        else
          call finish ('read_fdbk_gpsgb','GPSGB: no ZTD or STD')
        end if

        slants(k)% obs% time      = mjd(s% actual_time)! Modified Julian date of data
        slants(k)% obs% site      = 0                  ! GPS station ID
        slants(k)% obs% satellite = o% body (j)% lev_sig
        slants(k)% obs% elevation = o% olev (j)             *d2r ! elevation in radian
        slants(k)% obs% azimuth   = o% body (j)% obs_par(1) *d2r ! azimuth   in radian
        slants(k)% obs% slant     = o% body (j)% o               ! slant delay in m
        slants(k)% obs% zslant    = o% body (j)% o  ! STD mapped into zenith direction
      end do
      !------------------------------
      ! store in component t_obs% par
      !------------------------------
      call store_std (o, s, slants(1:n))
      call destruct (slants)
      deallocate (slants)
    end do

  end subroutine read_fdbk_gpsgb
!------------------------------------------------------------------------------
!==============================================================================
! Routines to handle the data structure used by the STD
! operators to store atmospheric fields and adjoints:
!
! Variables used:
!
!   cols  type(t_cols)        atmospheric data, gathered at observation PEs
!   xi    type(t_vector_segm) atmospheric data, interpolated
!   mc    type(p_column)      atmospheric data, conditioned for STD operator
!   spot  type(t_spot)        observation meta data, individual report
!   obs   type(t_obs_block)   observation meta data
!   s     type(SlantData)     slant observation meta data and data
!
! Auxiliary subroutines defined below:
!
!   std_col2xi                            convert t,q ,gp  to  t,rh,gp
!   std_2xi2col                           convert t,rh,gp  to  t,q ,gp
!==============================================================================
  subroutine set_std_data (mc, s, grid, spot, cols, xi, err)
    type(p_column)      ,intent(inout)        :: mc(:,:) ! model columns
    type(SlantData)     ,intent(in)           :: s       ! actual slant
    type (t_grid)       ,intent(in)           :: grid    ! atmospheric grid
    type(t_spot)        ,intent(in)           :: spot
    type(t_cols)        ,intent(in)           :: cols
    type(t_vector_segm) ,intent(in) ,optional :: xi
    integer             ,intent(out)          :: err

    integer                     :: i, j, k, l, ii
    integer                     :: ic, ncol
    integer                     :: npar, ke
    type(t_col)        ,pointer :: c
    real(wp)                    :: tmp
    real(wp),       allocatable :: t    (:) ! temperature
    real(wp),       allocatable :: q    (:) ! spec. humidity
    real(wp),       allocatable :: lnp  (:) ! (log) pressure
    real(wp),       allocatable :: geo  (:) ! geopotential
    real(wp),       allocatable :: d2t  (:) ! spline, 2nd derivative
    real(wp),       allocatable :: d2q  (:) ! spline, 2nd derivative
    real(wp),       allocatable :: d2lnp(:) ! spline, 2nd derivative
    type(p_column), allocatable :: tmpmc(:) ! temporary model columns
    real(wp),       parameter   :: eps = .01_wp ! lat/lon increment

    err  = 0
    ke   = spot% mke
    npar = 3 * ke
    ncol = size (spot% imcol)
    allocate (  t(ke),  q(ke),  lnp(ke), geo(ke))
    allocate (d2t(ke),d2q(ke),d2lnp(ke))
    allocate (tmpmc(ncol))

    if (present (xi)) then
      if (spot% n_spt * npar /= spot% i% n) &
        call finish ('set_std_data','GPSGB: n_spt*(3*ke)/=n')
    end if

    !-----------------------------------------------
    ! Check for fill values in ICON-LAM outer border
    !-----------------------------------------------
    if (grid% gridtype == DWD6_ICON .and. .not. grid% global) then
       do l = 1, ncol
          ic              =  spot% imcol(l)% imc(1)
          c               => cols% col(ic)
          tmp = c% t(1)
          !if (c% t(1) == 9999) then
          if (tmp == 9999) then
             ! Column with fill values => skip observation
             err = 1
             return
          end if
       end do
    end if

    do l = 1, ncol
       ic              =  spot% imcol(l)% imc(1)
       c               => cols% col(ic)
       allocate (tmpmc(l)% gpm (size(c% geo)))
       allocate (tmpmc(l)% p   (size(c% p)))
       allocate (tmpmc(l)% q   (size(c% q)))
       allocate (tmpmc(l)% t   (size(c% t)))
       tmpmc(l)% geoid =       c% s% geoid
       tmpmc(l)% dlon  =       c% c% dlon
       tmpmc(l)% dlat  =       c% c% dlat
       tmpmc(l)% gpm   =       c% geo / gacc     ! geometric height (fixed!)
       tmpmc(l)% p     =  exp (c% p)
       tmpmc(l)% q     =       c% q
       tmpmc(l)% t     =       c% t
       if (present (xi)) then
          ii = spot% i% i + (l-1)*npar
          call std_xi2col (t          (:),    & ! -> t
                           q          (:),    & ! -> q
                           geo        (:),    & ! -> geo
                           tmpmc(l)% p(:),    & ! <- p
                           xi% x(ii+1:ii+npar)) ! <- t,rh,geo
          tmpmc(l)% t(:) = t
          tmpmc(l)% q(:) = q
          tmpmc(l)% gpm(:) = geo(:) / gacc
!         if (l == 1) print *, "###gpm0:", &
!              c% geo(ke) / gacc,  tmpmc(l)% gpm(ke), &
!              c% geo(ke) / gacc - tmpmc(l)% gpm(ke)
          lnp = c% p !log (tmpmc(l)% p(:))
          call init_spline (geo, t, d2t, err)
          if (err /= 0) then
             ! An error in init_spline suggests that there are problems
             ! with the columns and that the code cannot proceed. In order to
             ! avoid divide by zero errors the program is stopped.
             call finish('set_std_data',                            &
                         'GPSGB: Geopotential not monotonously increasing')

             ! Debug columns: Pint columns with errors to ASCII files
             !call print_col(c, s)
             !cycle
          end if
          call init_spline (geo, q,   d2q)
          call init_spline (geo, lnp, d2lnp)
          do i = 1, ke
             call spline (geo, t,   d2t,   c% geo(i), tmpmc(l)% t(i))
             call spline (geo, q,   d2q,   c% geo(i), tmpmc(l)% q(i))
             call spline (geo, lnp, d2lnp, c% geo(i), tmp)
             tmpmc(l)% p(i) = exp (tmp)
          end do
          tmpmc(l)% gpm(:) = c% geo / gacc
       endif
    end do

    do k = 1, s% Naloc
       do j = 1, s% Nnghb
          l              =  s% idxi (j,k)
          ic             =  spot% imcol(l)% imc(1)
          c              => cols% col(ic)
          allocate (mc(j,k)% gpm (size(c% geo)))
          allocate (mc(j,k)% p   (size(c% p)))
          allocate (mc(j,k)% q   (size(c% q)))
          allocate (mc(j,k)% t   (size(c% t)))
          mc(j,k)% geoid = tmpmc(l)% geoid
          mc(j,k)% dlon  = tmpmc(l)% dlon
          mc(j,k)% dlat  = tmpmc(l)% dlat
          mc(j,k)% gpm   = tmpmc(l)% gpm
          mc(j,k)% p     = tmpmc(l)% p
          mc(j,k)% q     = tmpmc(l)% q
          mc(j,k)% t     = tmpmc(l)% t
          !------------------------------------------------------
          ! for ZTD nearest neighbour
          ! provide dummy coordinates surrounding the observation
          !------------------------------------------------------
          if (ncol == 1) then
            if (s% Nnghb >= 3) then
              mc(j,k)% dlon  = spot% col%c% dlon + (-1)** j    * eps
              mc(j,k)% dlat  = spot% col%c% dlat + (-1)**(j/2) * eps
            else
              mc(j,k)% dlon  = spot% col%c% dlon
              mc(j,k)% dlat  = spot% col%c% dlat
            endif
          endif
       end do
    end do

    do l = 1, ncol
       deallocate (tmpmc(l)% gpm)
       deallocate (tmpmc(l)% p  )
       deallocate (tmpmc(l)% q  )
       deallocate (tmpmc(l)% t  )
    end do

  end subroutine set_std_data
!------------------------------------------------------------------------------
  subroutine std_col2xi (t, q, geo, p, xi, fg)
    !------------------------------------------------------------------------
    ! converts columns of t,q  and geopotential (from different PEs)
    !                  to t,rh and geopotential (interpolation space)
    ! called by subroutine interpolate, module psasutil
    !------------------------------------------------------------------------
    real(wp) ,intent(in)  :: t  (:)  ! temperature         [K]
    real(wp) ,intent(in)  :: q  (:)  ! specific humidity   [kg/kg]
    real(wp) ,intent(in)  :: geo(:)  ! geopotential        [m2/s2]
    real(wp) ,intent(in)  :: p  (:)  ! pressure            [Pa]
    real(wp) ,intent(out) :: xi (:)  ! virtual temperature [K]
                                     ! relative humidity   [1]
                                     ! geopotential        [m2/s2]
    logical  ,intent(in)  :: fg      ! if true: account for rh0fg

    integer :: ke
    ke = size (t)
    xi(1:3*ke  :3) = t               ! temperature         [K]
    xi(2:3*ke  :3) = gh_rh(         &! generalised humidity[1]
                     rh_q (q, t, p) &! relative humidity
                   , fg, 200._wp, p )
    xi(3:3*ke  :3) = geo             ! geopotential height [gpm]
  end subroutine std_col2xi
!------------------------------------------------------------------------------
  subroutine std_xi2col (t, q, geo, p, xi)
  !------------------------------------------------------------------------
    ! converts columns of t,rh and geopotential  (interpolation space)
    !                  to t,rh and geopotential
    ! called by set_std_data (to set GNSS specific data type)
    !------------------------------------------------------------------------
    real(wp) ,intent(out) :: t  (:)  ! temperature         [K]
    real(wp) ,intent(out) :: q  (:)  ! specific humidity   [kg/kg]
    real(wp) ,intent(out) :: geo(:)  ! geopotential        [m2/s2]
    real(wp) ,intent(in)  :: p  (:)  ! pressure            [Pa]
    real(wp) ,intent(in)  :: xi (:)  ! virtual temperature [K]
                                     ! relative humidity   [1]
                                     ! geopotential        [m2/s2]
    integer  :: ke
    real(wp) :: rh (size(q))

    ke    = size (t)
    call rh_gh (rh, xi(2:3*ke  :3))
    t     =         xi(1:3*ke  :3)         ! temperature         [K]
    q     = q_rh   (rh, t, p )             ! specific humidity   [kg/kg]
    geo   =         xi(3:3*ke  :3)         ! geopotential        [m2/s2]

  end subroutine std_xi2col
!==============================================================================
  subroutine set_std (time)
    type(t_time) ,intent(in) :: time
    !-----------------------------------------------
    ! initialize auxiliary data used by STD operator
    !-----------------------------------------------
    integer                     :: month
    type(error_status) ,pointer :: errstat ! error status variable

    if (rept_use(OT_GPSGB)% use(CHK_NONE) <= STAT_DISMISS) return
    nullify (errstat)
    call enter_callee ('set_std',errstat)
    month = imm (time)
    call msis_init (month, errstat)
    if (error_callee (errstat)) then
      call display_status (errstat)
      call finish ('set_std','GPSGB: error in msis_init')
    endif
    call clear_status (errstat)
    if (associated (errstat)) deallocate (errstat)
  end subroutine set_std
!==============================================================================

  !---------------------------------------------------------------------
  ! subroutine print_col
  !---------------------------------------------------------------------
  !
  !> @brief Print model columns for debugging
  !>
  !> <b> call print_col (col, slant)  </b>
  !>
  !> @param[in]            col    One model column
  !> @param[in, optional]  slant  Slant observation
  !>
  !
  !---------------------------------------------------------------------
  subroutine print_col (col, slant)

    type(t_col), intent(in)                :: col     ! model column
    type (SlantData) ,intent(in), optional :: slant   ! slant meta data

    integer, save :: Nfile = 20
    integer :: k, j, i, Nlev, Ncol
    integer :: ounit
    character (len=12), allocatable :: cprof(:,:)
    character (len=30) :: fmt
    character (len=2)  :: ch2
    !-------------------------------------------------------------------
    Nfile = Nfile + 1

    write(*,*) 'Print profile to file fort ', dace% pe*100+Nfile

    ! Output unit, unit not open, write output files fort.unit, e.g. fort.1321
    ounit = dace% pe*100+Nfile

    write(ounit,*) '# Model column'
    write(ounit,*) '#'
    write(ounit,*) '# Written by print_col, mo_std'
    write(ounit,*) '#'

    if (present(slant)) then
       write(ounit,*) 'Delay , STD = ', slant%Obs%slant, ' m (obs)'
       write(ounit,*) '        at an elevation of ', slant%Obs%elevation*r2d,' degrees'
       write(ounit,*)
       write(ounit,*)                                                               &
            'Station ', slant% Station%SName,                              char(10),&
            'Station ID = ', slant% Station%ID,                            char(10),&
            'Geoid/Ellipsoid lon = ', slant% Station%CoordGeo(1)*r2d,           &
            slant% Station%CoordEll(1)*r2d,  char(10),&
            'Geoid/Ellipsoid lat = ', slant% Station%CoordGeo(2)*r2d,           &
            slant% Station%CoordEll(2)*r2d,  char(10),&
            'Geoid/Ellipsoid hei = ', slant% Station%CoordGeo(3) ,                  &
            slant% Station%CoordEll(3), char(10),         &
            'Cartesian X/Y/Z m   = ', slant% Station%CoordCart
       write(ounit,*)

       write(ounit,*)                                               &
        'Slant info: ',                              char(10),&
        ' Station ',    slant%Station%SName,         char(10),& ! station short name
        ' ID ',         slant%Station%ID,            char(10),& ! station ID
        ' time ',       slant%Obs%time,              char(10),& ! observation time
        ' Sat ',        slant%Obs%satellite,         char(10),& ! observed satellite
        ' Elev ',       slant%Obs%elevation*r2d,     char(10),& ! elevation
        ' Azi ',        slant%Obs%azimuth*r2d,       char(10),& ! azimuth
        ' ZTD_obs ',    slant%Obs%zslant,            char(10),& ! observed STD
        ' STD_obs ',    slant%Obs%slant,             char(10),& ! STD mapped to zenith
        ' STD_gerade ', slant%STDline,               char(10),& ! model STD, line
        ' STD_ray ',    slant%STD                               ! model STD, raytracer
       write(ounit,*)
    end if

    ! Print one column, several variables
    write(ounit,*)
    write(ounit,'(a)') 'Model column:'
    write(ounit,'(a,2(f8.3,tr2))') 'Coordinates, lon/lat, degrees: ', &
                                    col%c%dlon, col%c%dlat
    write(ounit,'(a,3(i8,tr2))')   'Grid indices i,j,l           : ', &
                                    col%i, col%j, col%l
    write(ounit,*)

    if (associated(col%t)) then
       Nlev = size(col%t)
    else if (associated(col%p)) then
       Nlev = size(col%p)
    else if (associated(col%geo)) then
       Nlev = size(col%geo)
    else
       ! Assume 65 level in ICON-D2
       Nlev = 65
    end if
    allocate( cprof(Nlev+2,10) )
    cprof = ''

    Ncol = 0

    ! First column with level index
    Ncol = Ncol + 1
    cprof(1,Ncol) = ' level'
    cprof(2,Ncol) = ' -----'
    do i=1, Nlev
       write(cprof(i+2,Ncol),'(i5)') i
    end do

    ! Pressure profile => is ln(p)!!!
    if (associated(col%p)) then
       Ncol = Ncol + 1
       cprof(1,Ncol) = '  ln(p)'
       cprof(2,Ncol) = '  -----'
       do i=1, Nlev
          write(cprof(i+2,Ncol),'(f10.3)') col%p(i)
       end do
    else
       write(ounit,'(a)') 'Pointer col%p not associated'
    end if

    ! Temperature profile
    if (associated(col%t)) then
       Ncol = Ncol + 1
       cprof(1,Ncol) = '  T [K]'
       cprof(2,Ncol) = '  -----'
       do i=1, Nlev
          write(cprof(i+2,Ncol),'(f10.3)') col%t(i)
       end do
    else
       write(ounit,'(a)') 'Pointer col%t not associated'
    end if

    ! Profile of specific humidity
    if (associated(col%q)) then
       Ncol = Ncol + 1
       cprof(1,Ncol) = '  q    '
       cprof(2,Ncol) = '  -----'
       do i=1, Nlev
          write(cprof(i+2,Ncol),'(f10.3)') col%q(i)
       end do
    else
       write(ounit,'(a)') 'Pointer col%q not associated'
    end if

    ! Profile of geopotential
    if (associated(col%geo)) then
       Ncol = Ncol + 1
       cprof(1,Ncol) = '  geo  '
       cprof(2,Ncol) = '  -----'
       do i=1, Nlev
          write(cprof(i+2,Ncol),'(f10.3)') col%geo(i)
       end do
    else
       write(ounit,'(a)') 'Pointer col%geo not associated'
    end if

    write(ounit,*)

    ! Print profiles in column
    write(ch2,'(i2)') Ncol
    fmt = '(' // ch2 // '(a12,tr2))'
    !write(ounit,*) fmt

    do i=1, Nlev+2
       write(ounit,fmt) cprof(i,1:Ncol)
    end do

  end subroutine print_col


 !---------------------------------------------------------------------
 ! subroutine print_field
 !---------------------------------------------------------------------
 !
 !> @brief Print internal Dace field to ASCII file (for debugging)
 !>
 !> <b> call print_field (atm, file, field, klev, msg) </b>
 !>
 !>
 !> @param[in]  atm    Field in t_atm
 !> @param[in]  file   File <=> not yet used
 !> @param[in]  field  Field to write (name from t_atm)
 !> @param[in]  klev   Level to write,
 !>                    klev=-1 - write all levels, one per file
 !> @param[in]  msg    Message, written to file header
 !>
 !> The model field is gathered to the I/O processor and written to
 !> an ASCII file. For each model level one separate ASCII file is
 !> written. The ASCII file contains 3 colums:
 !>   lon  lat  field
 !>
 !> The ASCII file names are
 !> Dace-field-FF-KKK.dat
 !> with field name FF and level KKK.
 !>
 !> The ASCII files might be plotted with the Python script
 !> PlotIconD2.py
 !>
 !> Requires
 !>
 !>  USE mo_atm_transp, ONLY: scatter_level  ! gather_multi     ! gather  multi level field
 !>
 !>  use mo_atm_state,   only: t_atm           ! atmospheric state derived type
 !>
 !> Example
 !> call print_field (atm(1), 'dummy', 't', -1, 'apply_H_m - atm(1)')
 !
 !---------------------------------------------------------------------
 subroutine print_field (atm, file, field, klev, msg)

    type(t_atm),   intent(in)    :: atm     ! atmospheric state
    CHARACTER (*), intent(in)    :: file
    CHARACTER (*), intent(in)    :: field
    integer,       intent(in)    :: klev
    CHARACTER (*), intent(in)    :: msg

    character (len=64) :: outf1, outfile
    character (len=3)  :: ch3
    INTEGER                      :: iu, ios
    INTEGER                      :: i, j, k, d, k1, k2
    REAL (kind=wp) _POINTER  :: ff(:,:,:,:)   => Null()  ! Local copy of full field

    !write(*,*) 'print_field start ...'
    !write(*,*) 'associated(atm%t) = ', associated(atm%t)

    if (dace% lpio) then
       allocate( ff(atm%grid%lbg(1):atm%grid%ubg(1), atm%grid%lbg(2):atm%grid%ubg(2), &
                    atm%grid%lbg(3):atm%grid%ubg(3), atm%grid%lbg(4):atm%grid%ubg(4)) )

!!$       write(*,*) 'size(ff), shape(ff) = ', size(ff), shape(ff)
!!$       write(*,*) 'size(atm%t), shape(atm%t) = ', size(atm%t), shape(atm%t)
!!$       write(*,*) 'size(atm%t(:,:,60,:)), shape(atm%t(:,:,60,:)) = ', &
!!$            size(atm%t(:,:,60,:)), shape(atm%t(:,:,60,:))
!!$       write(*,*) 'size(atm%grid%rlon), shape(atm%grid%rlon) = ', size(atm%grid%rlon), shape(atm%grid%rlon)
!!$       write(*,*) 'size(atm%grid%rlat), shape(atm%grid%rlat) = ', size(atm%grid%rlat), shape(atm%grid%rlat)

       !write(*,*) 'size(), shape() = ', size(), shape()


    end if

    if (field == 't') then
       outf1 = 'Dace-field-t'
       call gather_multi (ff, atm%t, atm%grid%dc, dace%pio)
    else if (field == 'pp') then
       outf1 = 'Dace-field-pp'
       call gather_multi (ff, atm%pp, atm%grid%dc, dace%pio)
    else if (field == 'q') then
       outf1 = 'Dace-field-q'
       call gather_multi (ff, atm%q, atm%grid%dc, dace%pio)
    else
       write(*,*) 'print_field: Unsupported field: ', trim(field)
       return
    end if

    !write(*,*) 'print_field nach gather '

  if (dace% lpio) then

     k = 60

     iu = get_unit_number ()

     !write(*,*) 'Open file field.dat, unit = ', iu

     ! Write one file per level
     if (klev > 0) then
        ! Write one file with this level
        k1 = klev
        k2 = klev
     else
        ! Write files for all levels
        k1 = atm%grid%lbg(3)
        k2 = atm%grid%ubg(3)
     end if

     do k=k1, k2
        write(ch3,'(i3.3)') k
        outfile = trim(outf1) // '-' // ch3 // '.dat'

        OPEN (iu, file=outfile,    &
             status='unknown', action='write', iostat=ios)
        if (ios == 0) then
           write(iu,'(a,a)') '# ', outfile
           write(iu,'(a)') '# ' // trim(msg)
           write(iu,'(a,4(i10,tr2))') '# shape(field)         = ', shape(ff)
           write(iu,'(a,4(i10,tr2))') '# shape(atm%grid%rlon) = ', shape(atm%grid%rlon)
           write(iu,'(a,l4)')         '# associated(atm%grid) = ', associated(atm%grid)
           write(iu,'(a,l4)')         '# associated(atm%grid%rlon) = ', associated(atm%grid%rlon)
           write(iu,'(a,l4)')         '# associated(atm%grid%rlat) = ', associated(atm%grid%rlat)
           write(iu,'(a)') '# '
           write(iu,'(a,i8,tr2,i6)') '# Verification time, JD, Sec. of day: ', &
                                        atm%time%days, atm%time%secs
           write(iu,'(a,i8,tr2,i6)') '# Reference time,    JD, Sec. of day: ', &
                                        atm%ref_time%days, atm%ref_time%secs
           write(iu,'(a)') '# '
           write(iu,'(a,2(f12.4,tr2))') '# Field min/max = ', minval(ff), maxval(ff)
           write(iu,'(a)') '# '
           write(iu,'(a,i3)') '# Level k = ', k
           write(iu,'(a)') '# '
           write(iu,'(a)') '# lon (deg.)  lat (deg.) field    '

           DO d     = LBOUND(ff,4),UBOUND(ff,4)
              DO j   = LBOUND(ff,2),UBOUND(ff,2)
                 DO i = LBOUND(ff,1),UBOUND(ff,1)

                    ! Write lon/lat/field to ASCII file
                    write(iu,'(3(f9.3,tr2))') &
                         atm%grid%rlon(i,j,1,d)*r2d, atm%grid%rlat(i,j,1,d)*r2d, &
                         ff(i,j,k,d)

                 END DO
              END DO
           END DO

           CLOSE (iu)

        else
           write(*,*) 'Error opening file field.dat, ios = ', ios
        end if

     end do   ! Loop levels / files

     CALL return_unit_number (iu)

  end if    ! if (dace% lpio) then

  end subroutine print_field

!==============================================================================
!
! specific MPI-bcast for derived type GNSSStation
!
#define VECTOR
#define DERIVED type(GNSSStation),dimension(:)
#define p_bcast_DERIVED p_bcast_GNSSStation
#undef MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
end module mo_std
