!
!+ GRIB inventory derived type definition and low level operators
!
MODULE mo_grib_invt
!
! Description:
!   GRIB inventory derived type definition.
!   Low level operators acting on this derived type.
!   Higher level routines are located in mo_grib_handling.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Andreas Rhodin
!  changes for GRIB2 API
! V1_20        2012-06-18 Andreas Rhodin
!  some cleanup
! V1_22        2013-02-13 Harald Anlauf
!  changes for GRIB API / CDI
! V1_23        2013-03-26 Harald Anlauf
!  GRIB_API implementation
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON
!  Add WMO code for MeteoSwiss
!  disable uuid for __ibm__
!  handle WMO4_3_HOURS for GRIB1
! V1_27        2013-11-08 Harald Anlauf
!  Fixes for generalized vertical coordinate; distinguish init_ana,init_fc
! V1_28        2014/02/26 Harald Anlauf
!  Fix GRIB2 for COSMO; analysis increment; clean up for ICON, IAU scheme
! V1_29        2014/04/02 Harald Anlauf
!  GRIB2: bugfix for (initialized) ensemble forecast
! V1_31        2014-08-21 Harald Anlauf
!  get_inventory: productDefinitionTemplateNumber table bugfix
! V1_35        2014-11-07 Harald Anlauf
!  grib2 inventory, typeOfStatisticalProcessing: handle min/max via 'range'
! V1_42        2015-06-08 Harald Anlauf
!  ICON local patch; adaptions to GRIB2
! V1_43        2015-08-19 Andreas Rhodin
!  Improve GRIB encoding of single-level/surface fields for flake
! V1_44        2015-09-30 Harald Anlauf
!  GRIB2 inventory: correctly detect analysis increments
! V1_45        2015-12-15 Harald Anlauf
!  Skip isobaric layers when deriving vertical grid; Changes for IFS
! V1_46        2016-02-05 Harald Anlauf
!  Handle DWD local definiton 252, GRIB2 encoding of mean and spread
! V1_47        2016-06-06 Harald Anlauf
!  inventory_2: bugfix for uninitialized p2; handle COSMO sub-hourly forecasts
!   handle initialized analysis for ECMWF ensembles (pf/cf)
! V1_49        2016-10-25 Harald Anlauf
!  Generalized vertical coordinate: use nlev from GRIB2
! V1_50        2017-01-09 Harald Anlauf
!  Scaling factor for ECMWF grib edition 1 fields from mars
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD    2010  initial version
! Harald Anlauf   DWD    2012  extensions for GRIB2, ECMWF
!------------------------------------------------------------------------------

  !=============
  ! Modules used
  !=============
  use mo_kind,       only: i8, dp             ! 8 byte integer kind parameter
  use mo_exception,  only: finish, message    ! abort on error condition
  use mo_time,       only: t_time,           &! date+time data type
                           init_time,        &! time initialisation routine
                           operator(+),      &! add      times
                           operator(==),     &! compare  times
                           ZERO_TIME          ! zero time value
  use mo_emos_grib1, only: t_grib1,          &! GRIB record data type
                           local_ecmwf,      &! ECMWF local extension present ?
                           local_dwd,        &! DWD   local extension present ?
                           INVALID_HANDLE     ! invalid GRIB2-API handle value
  use mo_gribtables, only: ar_des,           &! GRIB entry descriptor type
                           search             ! search for grib entries
  use mo_wmo_tables, only: WMO0_ECMWF,       &!
                           WMO0_DWD,         &!
                           WMO0_COSMO,       &!
                           WMO0_MSWISS,      &!
                           WMO0_COMET,       &!
                           DWD6_ICOSAHEDRON, &!
                           DWD6_ICON,        &!
                           WMO6_LATLON,      &!
                           WMO6_ROTLL,       &!
                           WMO6_GAUSSIAN,    &!
                           WMO6_HARMONIC,    &!
                           WMO4_SECOND,      &!
                           WMO4_MINUTE,      &!
                           DWD4_15_MINUTES,  &!
                           WMO4_HOUR,        &!
                           WMO4_3_HOURS,     &!
                           WMO4_12_HOURS,    &!
                           WMO4_DAY,         &!
                           WMO4_MONTH,       &!
                           WMO3_SURFACE,     &!
                           WMO3_CLD_BASE,    &!
                           WMO3_CLD_TOPS,    &!
                           WMO3_ZERO_ISO,    &!
                           WMO3_NOM_TOA,     &!
                           WMO3_ISOTHERM,    &!
                           WMO3_SEALEVEL,    &!
                           WMO3_ISOBARIC,    &!
                           WMO3_ABOVESUR,    &!
                           WMO3_BELOWSUR,    &!
                           WMO3_HEIGHT,      &! Altitude above MSL [m]
                           WMO3_LAYER,       &!
                           WMO3_HYBRID,      &!
!                          WMO3_HHYBRID,     &!
                           WMO3_HYBRIDB,     &!
                           WMO3_GENV,        &!
                           WMO3_LAKE_BOT,    &!
                           WMO3_SEDIM_BOT,   &!
                           WMO3_MIX_LAYER,   &!
                           WMO3_SURF_HORZ,   &!
                           WMO3_SURF_TANG,   &!
                           DWD5_ENKF_ANA,    &!
                           DWD5_IFS_FC,      &!
                           DWD5_NUDGING,     &!
                           WMO5_ACCU,        &!
                           WMO5_AVERAGE,     &!
                           WMO5_DIFF,        &!
                           WMO5_FORECAST,    &!
                           WMO5_RANGE,       &!
                           WMO5_INIT_ANA,    &!
                           DWD5_INIT_FC,     &!
                           WMO5_VALID         !
! use mo_run_params, only: grib_library       ! GRIB API to use
#ifdef GRIB_API
  use grib_api,       only: grib_get,        &! get item from GRIB message
                            grib_set,        &! set item in   GRIB message
                            GRIB_SUCCESS      ! return status value
  use mo_grib12_dwd,  only: t_par_grib1,     &! Derived type: GRIB1 params.
                            t_par_grib2,     &! Derived type: GRIB2 params.
                            search_grib1,    &! Search param. in GRIB1 tab.
                            search_grib2      ! Search param. in GRIB2 tab.
#endif
  implicit none

  !================
  ! Public entities
  !================
  private
  PUBLIC :: grib_library     ! GRIB API to use
  !-------------------------------------------------
  ! data type definition: t_inventory and components
  !-------------------------------------------------
  public :: t_inventory      ! Inventory Data Type
  public ::   t_ds           !   characteristics of the data set
  public ::   t_ct           !   identification of the center
  public ::   t_pa           !   characteristics of the parameter
  public ::   t_ti           !   time information
  public ::   t_lv           !   level information
  public ::   t_gr           !   horizontal grid
  public ::   t_en           !   ensemble data
  !--------------------------------------------------
  ! routines to derive the inventory from a GRIB file
  !--------------------------------------------------
  public :: inventory_entry  ! derive inventory table entry from GRIB
  PUBLIC :: operator (==)    ! compare components of an inventory entry
  !---------------------------------------
  ! derive quantities from a GRIB 1 record
  !---------------------------------------
  PUBLIC :: ver_time         ! get verification time from grib1 record
  PUBLIC :: ref_time         ! get reference time from grib1 record
  PUBLIC :: db_time          ! get data-base time from grib1 record
  !----------------------------
  ! ensemble derived quantities
  !----------------------------
  PUBLIC :: ENS_MEAN         ! ensemble mean
  PUBLIC :: ENS_SPREAD       ! ensemble spread

  !====================
  ! Inventory Data Type
  !====================
                                     ! characteristics of the data set
  TYPE t_ds                          !
    INTEGER           :: rec         !   record number
#ifdef USE_PBSEEK64
    INTEGER(i8)       :: ptr         !   address of record in file
#else
    INTEGER           :: ptr         !   address of record in file
#endif
    INTEGER           :: blen        !   size of GRIB record in bytes
    TYPE (t_time)     :: db_time     !   data-base time
    INTEGER           :: edition     !   GRIB edition number
    INTEGER           :: packing     !   Data representation (packing)
  END TYPE t_ds                      !
                                     ! identification of the center
  TYPE t_ct                          !
    INTEGER           :: center      !   PDS octet 5
    INTEGER           :: subcenter   !   PDS octet 26
    INTEGER           :: process     !   PDS octet 6
    INTEGER           :: local_ident !   local extension identification
  END TYPE t_ct                      !
                                     ! identification of the parameter
  TYPE t_pa                          !
    CHARACTER(len=24) :: shortname   !   parameter name (postprocessing)
    CHARACTER(len=10) :: iname       !   internally used name
!   CHARACTER(len=8)  :: units       !   units
                                     ! GRIB1 table and code:
    INTEGER           :: table       !   parameter table version (PDS octet 4)
    INTEGER           :: code        !   kpds5, PDS octet 9
                                     ! GRIB2 triple:
    INTEGER           :: discipline  !   discipline
    INTEGER           :: category    !   parameterCategory
    INTEGER           :: number      !   parameterNumber
                                     !
    CHARACTER(len=8)  :: runtype     !   analysis, forecast, etc.
    INTEGER           :: runclass    !   0..3 haupt,vor,ass,test
    INTEGER           :: expid       !   experiment id
    INTEGER           :: zen         !   additional element number
    REAL(dp)          :: factor      !   scaling factor for units conversion
  END TYPE t_pa                      !
                                     ! time information
  TYPE t_ti                          !
    TYPE (t_time)     :: ver_time    !   verification time
    TYPE (t_time)     :: ref_time    !   reference time
    TYPE (t_time)     :: rng_time    !   time range
    INTEGER           :: p1          !   PDS octet 19
    INTEGER           :: p2          !   PDS octet 20
    INTEGER           :: range       !   PDS octet 21
    INTEGER           :: unit        !   PDS octet 14
    INTEGER           :: range_unit  !   GRIB2: indicatorOfUnitForTimeRange
  END TYPE t_ti                      !
                                     ! level information
  TYPE t_lv                          !
    INTEGER           :: leveltype   !   kdps6, PDS octet 10
    INTEGER           :: levelvalue  !   PDS octets 11-12
    INTEGER           :: levels(2)   !   PDS octet 11, PDS octet 12
    INTEGER           :: nlev        !   no. levels for gen.vert.coordinate
    INTEGER           :: grid_num    !   number of vertical grid used
    CHARACTER(len=1 ) :: uuid(16)    !   UUID of vertical grid
  END TYPE t_lv                      !
                                     ! horizontal grid
  TYPE t_gr                          !
    INTEGER           :: gridtype    !   GDS octet 6
    INTEGER           :: ni          !   characteristic for resolution
    INTEGER           :: nxny        !   number of values in BDS
    INTEGER           :: grid_num    !   number of grid used
    INTEGER           :: grid_ref    !   number of grid in reference
    CHARACTER(len=1 ) :: uuid(16)    !   UUID of unstructured grids
  END TYPE t_gr                      !
                                     ! ensemble
  TYPE t_en                          !
    INTEGER           :: id          ! Ensemble identification by table
    INTEGER           :: size        ! Number of ensemble members
    INTEGER           :: no          ! Actual number of ensemble member
  END TYPE t_en                      !
                                     !---------------------------------
                                     ! container
  TYPE t_inventory                   !
    TYPE (t_ds)       :: ds          !   characteristics of the data set
    TYPE (t_ct)       :: ct          !   identification of the center
    TYPE (t_pa)       :: pa          !   identification of the parameter
    TYPE (t_ti)       :: ti          !   time information
    TYPE (t_lv)       :: lv          !   level information
    TYPE (t_gr)       :: gr          !   horizontal grid
    TYPE (t_en)       :: en
    !--------------------------
    ! relations to other fields
    !--------------------------
    INTEGER           :: ifile       ! file index
    INTEGER           :: nvar        ! number of variables
    INTEGER           :: nlev        ! number of levels
    INTEGER           :: ntime       ! number of time slots
    INTEGER           :: nmem        ! number of ensemble members
    INTEGER           :: jlev        ! level number
    INTEGER           :: jtime       ! time slot
    INTEGER           :: ivar        ! variable count
    INTEGER           :: jmem        ! ensemble member count
  END TYPE t_inventory

  !----------------------------
  ! ensemble derived quantities
  !----------------------------
  INTEGER, PARAMETER :: ENS_MEAN    =  -2     ! ensemble mean
  INTEGER, PARAMETER :: ENS_SPREAD  =  -3     ! ensemble spread

  !===========
  ! Interfaces
  !===========

  INTERFACE local_dwd
    MODULE PROCEDURE local_dwd_invt
  END INTERFACE local_dwd

  INTERFACE local_ecmwf
    MODULE PROCEDURE local_ecmwf_invt
  END INTERFACE local_ecmwf

  !-------------------
  ! compare data types
  !-------------------
  INTERFACE OPERATOR (==)
!   MODULE PROCEDURE equal_ct
    MODULE PROCEDURE equal_pa
    MODULE PROCEDURE equal_ti
    MODULE PROCEDURE equal_lv
!   MODULE PROCEDURE equal_gr
    MODULE PROCEDURE equal_en
    MODULE PROCEDURE equal_invt
  END INTERFACE

  !=================
  ! module variables
  !=================
#ifdef GRIB_API
  integer :: grib_library = 2   ! API: 1=(C)GRIBEX  2=GRIB2-API (default)
#else
  integer :: grib_library = 1   ! API: 1=(C)GRIBEX
#endif

contains
!==============================================================================

  FUNCTION inventory_entry (grib) result (gi)
  TYPE (t_grib1), intent(in) :: grib  ! GRIB 1 record
  TYPE (t_inventory)         :: gi    ! table entry
    !-------------------------------------------------------
    ! derives a GRIB file inventory entry from a GRIB record
    !-------------------------------------------------------
    select case (grib_library)
    case (1)
      gi = inventory_1 (grib)
    case (2)
      gi = inventory_2 (grib)
    end select
  END FUNCTION inventory_entry

!------------------------------------------------------------------------------

  FUNCTION inventory_1 (grib) result (gi)
  TYPE (t_grib1), intent(in) :: grib  ! GRIB 1 record
  TYPE (t_inventory)         :: gi    ! table entry
  !-------------------------------------------------------
  ! derives a GRIB file inventory entry from a
  ! GRIB1 record decoded by (C)GRIBEX
  !-------------------------------------------------------

    TYPE (ar_des) :: entry  ! GRIB table entry

    gi%ds% rec         = grib% krec
    gi%ds% ptr         = grib% kpos
    gi%ds% blen        = grib% blen
    gi%ds% edition     = grib% isec0% edition
    gi%ds% db_time     = db_time (grib)
    gi%ds% packing     = 0              ! currently not handled

    gi%ct% center      = grib% isec1% center
    gi%ct% subcenter   = grib% isec1% sub_center
    gi%ct% process     = grib% isec1% process
    gi%ct% local_ident = grib% isec1% local_ident

    gi%ti% ver_time    = ver_time (grib)
    gi%ti% ref_time    = ref_time (grib)
    gi%ti% unit        = grib% isec1% time_unit
    gi%ti% range_unit  = grib% isec1% time_unit
    gi%ti% p1          = grib% isec1% p1
    gi%ti% p2          = grib% isec1% p2
    gi%ti% range       = grib% isec1% time_range
    call init_time (gi% ti% rng_time, hh = 0)

    gi%pa% shortname   = ''
    gi%pa% table       = grib% isec1% table
    gi%pa% code        = grib% isec1% code
    gi%pa% discipline  = 255
    gi%pa% category    = 255
    gi%pa% number      = 255
    gi%pa% runtype     = 'unknown '

    IF (gi%ti% range == WMO5_FORECAST .AND. gi%ti% p1 == 0) gi%pa% runtype = 'analysis'
    IF (gi%ti% range == WMO5_FORECAST .AND. gi%ti% p1 /= 0) gi%pa% runtype = 'forecast'
    IF (gi%ti% range == DWD5_INIT_FC)                       gi%pa% runtype = 'init_fc'
    IF (gi%ti% range == WMO5_INIT_ANA)                      gi%pa% runtype = 'init_ana'
    IF (gi%ti% range == WMO5_RANGE)                         gi%pa% runtype = 'range   '
    IF (gi%ti% range == WMO5_AVERAGE)                       gi%pa% runtype = 'average '
    IF (gi%ti% range == WMO5_ACCU)                          gi%pa% runtype = 'accum   '
    IF (gi%ti% range == WMO5_VALID)                         gi%pa% runtype = 'valid   '
    IF (gi%ti% range == DWD5_NUDGING)                       gi%pa% runtype = 'nudging '
    IF (gi%ti% range == DWD5_IFS_FC)                        gi%pa% runtype = 'ifs_fc  '
    IF (gi%ti% range == DWD5_ENKF_ANA)                      gi%pa% runtype = 'enkf_ana'

    gi%lv% leveltype  = grib% isec1% level_type
    gi%lv% levelvalue = grib% isec1% level_st
    gi%lv% levels     = (/grib% isec1% level_st, grib% isec1% level_b/)
    !--------------------------------------------------------------------
    ! For soil variables, derive levelvalue from mean of layer boundaries
    ! (rounded "up").
    !--------------------------------------------------------------------
    if (gi%lv% leveltype == WMO3_BELOWSUR .and. grib% isec1% level_b > 0) then
       gi%lv% levelvalue = (grib% isec1% level_st + grib% isec1% level_b + 1) / 2
    end if
!if (gi%lv% leveltype == WMO3_BELOWSUR) then
!print *, "Belowsurf:"
!print *, "Levelvalue, Levels(:) =", gi%lv% levelvalue, gi%lv% levels(:)
!end if

    entry             = search (  code= gi%pa% code,    &
                                 table= gi%pa% table,   &
                                levtyp= gi%lv% leveltype)
    if (entry% name == '') &
      entry           = search (  code= gi%pa% code,    &
                                 table= gi%pa% table    )
    gi%pa% shortname  = entry% name
    gi%pa% iname      = entry% iname

    gi%pa% factor     = 1._dp
    !-----------------------------------------------------------
    ! Save conversion factor for selected ECMWF fields from mars
    !-----------------------------------------------------------
    if (gi%ct% center == WMO0_ECMWF .and. gi%ds% edition == 1     &
                                    .and. entry% fak     /= 1._dp &
                                    .and. entry% units   /= ""    ) then
       gi%pa% factor  = entry% fak
    end if

    gi%gr% nxny       = grib% isec4% n_data
    gi%gr% gridtype   = grib% isec2 (1)
    gi%gr% ni         = 0

    SELECT CASE (gi%gr% gridtype)
    CASE (DWD6_ICOSAHEDRON)
      gi%gr% ni = grib% tri%    ni
    CASE (WMO6_LATLON, WMO6_ROTLL)
      gi%gr% ni = grib% latlon% ni
    CASE (WMO6_GAUSSIAN)
      gi%gr% ni = grib% gauss%  nj * 2
    CASE DEFAULT
      write(0,*) 'get_inventory: invalid gridtype =',gi%gr% gridtype
      call finish('get_inventory','invalid gridtype')
    END SELECT

    gi%pa% runclass =   0
    gi%pa% expid    =   0 ! experiment id
    gi%pa% zen      = 255 ! additional element number
    gi%en% id       =   0 ! Ensemble identification by table
    gi%en% size     =   0 ! Number of ensemble members
    gi%en% no       =  -1 ! Actual number of ensemble member
    if (local_dwd (grib)) then
      if (gi% ct% local_ident == 253) then
        gi%en% id     = grib% s1_dwd% ensemble_id
        gi%en% size   = grib% s1_dwd% ensemble_size
        gi%en% no     = grib% s1_dwd% ensemble_no
      endif
      gi%pa% runclass = grib% s1_dwd% run_type
      gi%pa% expid    = grib% s1_dwd% exp
      gi%pa% zen      = grib% s1_dwd% element_no
    else if (local_ecmwf (grib)) then
      if (grib% s1_ecmwf% n_en > 1) then
        gi%en% id   = grib% s1_ecmwf% type    ! Grib type (11=perturbed fc)
        gi%en% size = grib% s1_ecmwf% n_en    ! Ensemble size incl. control
        gi%en% no   = grib% s1_ecmwf% i_en    ! Member id
      end if
      if (grib% s1_ecmwf% type == 2)   gi%pa% runclass = 2

      if (verify (grib% s1_ecmwf% expid,'0123456789') == 0) then
        read (grib% s1_ecmwf% expid,*) gi%pa% expid   ! ASCII->int
!     else
!       write (0,*) 'inventory_1: bad expid=',grib% s1_ecmwf% expid
      end if
    else
    endif

    !-------------------------------------------------------------------
    ! Fix wrong coding of time range for selected ECMWF fields from mars
    !-------------------------------------------------------------------
    if (gi%ct% center == WMO0_ECMWF .and. gi%ds% edition == 1) then
       if (gi%pa% iname == 'tot_prec') then
          gi%ti% range = WMO5_ACCU
          gi%ti% p2 = max (gi%ti% p1, gi%ti% p2)
          gi%ti% p1 = 0
       end if
    end if
    !---------------
    ! set time range
    !---------------
    select case (gi% ti% range)
    case (WMO5_RANGE,WMO5_AVERAGE,WMO5_ACCU)
      gi% ti% rng_time = time_p12 (gi% ti% p2 - gi% ti% p1, gi% ti% range_unit)
    end select

  END FUNCTION inventory_1

!------------------------------------------------------------------------------

  FUNCTION inventory_2 (grib) result (gi)
  TYPE (t_grib1), intent(in) :: grib  ! GRIB record
  TYPE (t_inventory)         :: gi    ! table entry
  !-------------------------------------------------------
  ! derives a GRIB file inventory entry from a
  ! GRIB1 record using the GRIB2 API
  !-------------------------------------------------------

#ifdef GRIB_API

    integer           :: stat, ffs, sfs, scf, scv
    integer           :: topd, togp, pdtn, tosp
    integer           :: relscale       ! Relative scale
    integer           :: i, ms, mt, tr
    integer           :: status
    logical           :: getlev
    character(len=8)  :: c
    TYPE (ar_des)     :: entry          ! GRIB table entry
    type(t_par_grib1) :: code1          ! GRIB1 parameters
    type(t_par_grib2) :: code2          ! GRIB2 parameters
    character(len=1)  :: uuid(16)       ! Grid UUID
    logical           :: lrange         ! Time range?
    integer           :: runit          ! Range units

    !--------------------------------
    ! characteristics of the data set
    !--------------------------------
    gi%ds% rec        = grib% krec
    gi%ds% ptr        = grib% kpos
    gi%ds% blen       = grib% blen
    call grib_get (grib% handle, 'editionNumber',gi%ds% edition)
    gi%ds% db_time    = db_time (grib)
    if (gi%ds% edition == 2) then
       call grib_get (grib% handle, 'dataRepresentationTemplateNumber', &
                                    gi%ds% packing)
    else
       gi%ds% packing = 0           ! currently not handled
    end if

    !-----------------------------
    ! identification of the center
    !-----------------------------
    call grib_get (grib% handle,'centre',                     gi%ct% center)
    call grib_get (grib% handle,'subCentre',                  gi%ct% subcenter)
    call grib_get (grib% handle,'generatingProcessIdentifier',gi%ct% process)
    call grib_get (grib% handle,'localDefinitionNumber'      ,gi%ct% local_ident,&
                   status)
    if (status /= GRIB_SUCCESS) gi%ct% local_ident = 0
    !-----------------
    ! time information
    !-----------------
!   gi%ti% ver_time    = ver_time (grib)
    gi%ti% ref_time    = ref_time (grib)
    if (gi%ds% edition == 2) then
      ! Code table 1.2
      call grib_get (grib% handle,'significanceOfReferenceTime',i)
      if (i < 0 .or. i > 1) then
         write(*,*) 'significanceOfReferenceTime =', i
         call finish ('get_inventory','bad significanceOfReferenceTime')
      end if
    end if

!   call grib_get (grib% handle,'stepType',i)
!   call grib_get (grib% handle,'stepType',c)
    call init_time (gi% ti% rng_time, hh = 0)

    if (gi%ds% edition == 1) then
!     call grib_get (grib% handle,'stepUnits'         ,gi%ti% unit)
      call grib_get (grib% handle,'unitOfTimeRange'   ,gi%ti% unit)
      gi%ti% range_unit = gi%ti% unit
      if (gi%ti% unit == DWD4_15_MINUTES) then
         !--------------------------------------------------------------
         ! Enforce conversion to minutes to work around GRIB_API problem
         !--------------------------------------------------------------
         call grib_set (grib% handle,'stepUnits'      ,WMO4_MINUTE)
         gi%ti% unit = WMO4_MINUTE
      end if
      call grib_get (grib% handle,'startStep'         ,gi%ti% p1)
      call grib_get (grib% handle,'endStep'           ,gi%ti% p2)
      call grib_get (grib% handle,'timeRangeIndicator',gi%ti% range)
    end if
    gi%ti% ver_time  = ver_time (grib)
    !-----------------
    ! derive 'runtype'
    !-----------------
    if (gi%ds% edition > 1) then
      ! Code table 1.4
      call grib_get (grib% handle,'typeOfProcessedData'            ,topd)
      ! Code table 4.0
      call grib_get (grib% handle,'productDefinitionTemplateNumber',pdtn)

      lrange = .false.
      select case (pdtn)
      case (0,1,7,32,33,40,41,44,45,48) ! Analysis or forecast
      case (57)                         ! Analysis or forecast (aerosols etc.)
      case (2)                          ! Derived forecast
      case (8,11,12,15,34,42,43,46,47)  ! Average or accumulated values
         lrange = .true.
      case default
         write(*,*) 'productDefinitionTemplateNumber=', pdtn
         call finish ('get_inventory','unsupported template')
      end select

      ! Code table 4.3
      call grib_get (grib% handle,'typeOfGeneratingProcess',togp)
      select case (togp)
      case (0)   ! Analysis
      case (1)   ! Initialization
      case (2)   ! Forecast
      case (3)   ! Bias corrected forecast
      case (4)   ! Ensemble forecast
      case (5)   ! Probability forecast
      case (6)   ! Forecast error
      case (7)   ! Analysis error
      case (8)   ! Observation
      case (9)   ! Climatological
      case (10)  ! Probability-weighted forecast
      case (11)  ! Bias-corrected ensemble forecast
      case (12)  ! Post-processed analysis
      case (13)  ! Post-processed forecast
      case (14)  ! Nowcast
      case (15)  ! Hindcast
      case (16)  ! Physical retrieval
      case (17)  ! Regression analysis
      case (18)  ! Difference between two forecasts
      case (192) ! bias corrected ensemble forecast (DWD)
      case (193) ! calibrated forecast (DWD)
      case (194) ! calibrated ensemble forecast (DWD)
      case (195) ! Interpolated analysis/forecast (e.g. ECMWF background)
      case (196) ! Invariant data (DWD)
      case (197) ! smoothed forecast (DWD)
      case (198) ! smoothed and calibrated forecast (DWD)
      case (199) ! probability forecast derived by neighbourhood method (DWD)
      case (200) ! Difference init. analysis - analysis
      case (201) ! Difference analysis - first guess
      case (202) ! Nudging
      case (203) ! Nudgcast
      case (204) ! product derived by statistical model (DWD)
      case (205) ! Deterministic forecast from ensemble mean (DWD)
      case (206) ! Time-filtered bias (DWD)
      case (207) ! Time-filtered assimilation increment (DWD)
      case (208) ! Time-filtered increment, weighted with cos(2pi loc.time/day)
      case (220) ! postprocessing (DWD)
      case (209:219) ! DWD reserved
      case (221:254) ! DWD reserved
      case default
         write(*,*) 'typeOfGeneratingProcess =', togp
         call finish ('get_inventory','unsupported Process')
      end select

      ! Code table 4.4
      call grib_get (grib% handle,'indicatorOfUnitOfTimeRange',gi%ti% unit)
      if (gi%ti% unit == 13) gi%ti% unit = WMO4_SECOND
      call grib_get (grib% handle,'forecastTime',gi%ti% p1)
!     call grib_get (grib% handle,'stepRange'   ,gi%ti% p2)

      if (lrange) then
         !---------------
         ! get time range
         !---------------
         ! Code table 4.10: Type of statistical processing
         call grib_get (grib% handle,'typeOfStatisticalProcessing',tosp)

         ! Code table 4.4
         call grib_get (grib% handle,'indicatorOfUnitForTimeRange',runit)
         if (runit == 13) runit = WMO4_SECOND
         gi%ti% range_unit = runit
         call grib_get (grib% handle,'lengthOfTimeRange'          ,tr)
         gi%ti% rng_time = time_p12 (tr, gi% ti% range_unit)
         gi%ti% p2       = gi%ti% p1 + tr       ! Derive endStep
!        if (gi%ti% p1 == 0) then
!           gi%ti% unit = runit
!        else if (runit /= gi%ti% unit) then
!           write (*,*) 'Inconsistent units for time range:', gi%ti% unit, runit
!           call finish ('get_inventory','unsupported: inconsistent time units')
!        end if
      else
         tosp      = -1
         gi%ti% p2 = gi%ti% p1                  ! Emulate endStep=startStep
      end if

      select case (tosp)
      case (-1)
         gi%ti% range = WMO5_FORECAST
      case (0)
         gi%ti% range = WMO5_AVERAGE
      case (1,11)       ! Treat summation as accumulation
         gi%ti% range = WMO5_ACCU
      case (2:3)        ! Maximum, minimum
         gi%ti% range = WMO5_RANGE
      case (4)
         gi%ti% range = WMO5_DIFF
      case default
         write(*,*) 'typeOfStatisticalProcessing =', tosp
         call finish ('get_inventory','unsupported statistical processing')
      end select

      if (.not. lrange) then
         select case (topd)
         case (0)       ! Analysis
            if (togp == 202) gi%ti% range = DWD5_NUDGING
         case (1,3:5)   ! Forecast
            if (gi%ti% p1 == 0 .and. togp /= 0  &
                               .and. togp /= 201) then
                             gi%ti% range = WMO5_INIT_ANA
            else
                             gi%ti% range = WMO5_FORECAST
              if (togp == 1) gi%ti% range = DWD5_INIT_FC
            end if
            if (togp == 195) gi%ti% range = DWD5_IFS_FC
         case (2)       ! Analysis and forecast
!        case (3)       ! Control forecast
!        case (4)       ! Perturbed forecast
!        case (5)       ! Control and perturbed forecast
         case (192)     ! Perturbed analysis
         case (255)     ! Missing value
         case default
            write(*,*) 'typeOfProcessedData =', topd
            call finish ('get_inventory','unsupported typeOfProcessedData')
         end select
      end if
    end if ! edition == 2

    gi%pa% runtype    = 'unknown '

    IF (gi%ti% range == WMO5_FORECAST .AND. gi%ti% p1 == 0) then
       if (togp == 201) then
                                             gi%pa% runtype = 'ana_inc'
       else
                                             gi%pa% runtype = 'analysis'
       end if
    end IF
    IF (gi%ti% range == WMO5_FORECAST .AND. gi%ti% p1 /= 0) &
                                             gi%pa% runtype = 'forecast'
    IF (gi%ti% range == DWD5_INIT_FC)        gi%pa% runtype = 'init_fc'
    IF (gi%ti% range == WMO5_INIT_ANA)       gi%pa% runtype = 'init_ana'
    IF (gi%ti% range == WMO5_RANGE)          gi%pa% runtype = 'range   '
    IF (gi%ti% range == WMO5_DIFF)           gi%pa% runtype = 'diff    '
    IF (gi%ti% range == WMO5_AVERAGE)        gi%pa% runtype = 'average '
    IF (gi%ti% range == WMO5_ACCU)           gi%pa% runtype = 'accum   '
    IF (gi%ti% range == DWD5_NUDGING)        gi%pa% runtype = 'nudging '
    IF (gi%ti% range == DWD5_IFS_FC)         gi%pa% runtype = 'ifs_fc  '
    IF (gi%ti% range == DWD5_ENKF_ANA)       gi%pa% runtype = 'enkf_ana'

    !--------------------------------
    ! identification of the parameter
    !--------------------------------
    call grib_get (grib% handle,'shortName'             ,gi%pa% shortname)
    gi%pa% discipline  = 255
    gi%pa% category    = 255
    gi%pa% number      = 255
    gi%pa% table       = 255
    gi%pa% code        = 255
    if (gi%ds% edition == 1) then
      call grib_get (grib% handle,'table2Version'       ,gi%pa% table)
      call grib_get (grib% handle,'indicatorOfParameter',gi%pa% code)
    else
      call grib_get (grib% handle,'discipline'          ,gi%pa% discipline)
      call grib_get (grib% handle,'parameterCategory'   ,gi%pa% category)
      call grib_get (grib% handle,'parameterNumber'     ,gi%pa% number)
    end if

    !--------------------------------------------------
    ! Fixup of effective runtype for selected variables
    !--------------------------------------------------
    if (gi%ti% range == WMO5_AVERAGE) then
      select case (gi%pa% shortname)
      case ("U_10M_AV", "V_10M_AV")     ! 10-minute averaged 10m wind
        gi%pa% runtype = 'forecast'
      end select
    end if

    !------
    ! level
    !------
    ffs           = -1
    sfs           = -1
    getlev        = .false.
    gi%lv% levels = -1

    if (gi%ds% edition == 1) then
      call grib_get (grib% handle,'indicatorOfTypeOfLevel' ,gi%lv% leveltype)
!     call grib_get (grib% handle,'typeOfLevel', i)
!     call grib_get (grib% handle,'typeOfLevel', c)
      select case (gi%lv% leveltype)
      case (WMO3_ISOBARIC,WMO3_ABOVESUR,WMO3_BELOWSUR,WMO3_HYBRID,WMO3_HEIGHT)
         call grib_get (grib% handle,'level'      ,gi%lv% levels(1))
      case (WMO3_HYBRIDB)
         call grib_get (grib% handle,'topLevel'   ,gi%lv% levels(1))
         call grib_get (grib% handle,'bottomLevel',gi%lv% levels(2))
      end select
    else
      call grib_get (grib% handle,'typeOfFirstFixedSurface' ,ffs, status=status)
      if (status /= GRIB_SUCCESS) ffs = -1
      call grib_get (grib% handle,'typeOfSecondFixedSurface',sfs, status=status)
      if (status /= GRIB_SUCCESS) sfs = -1
      relscale = 0
      select case (ffs)
      case (1)
         gi%lv% leveltype    = WMO3_SURFACE  ! 1
         select case (sfs)
         case (WMO3_LAKE_BOT,WMO3_MIX_LAYER)
            gi%lv% leveltype = sfs
         end select
      case (100)
         gi%lv% leveltype    = WMO3_ISOBARIC ! 100
         getlev              = .true.
         relscale            = -2            ! Pa -> hPa
      case (101)
         gi%lv% leveltype    = WMO3_SEALEVEL ! 102
      case (102)
         gi%lv% leveltype    = WMO3_HEIGHT   ! 103
         getlev              = .true.
      case (103)
         select case (sfs)
         case default
            gi%lv% leveltype = WMO3_ABOVESUR ! 105
         case (1:180,185:254)
            gi%lv% leveltype = WMO3_LAYER    ! 106
         end select
         getlev              = .true.
      case (106)
         gi%lv% leveltype    = WMO3_BELOWSUR ! 111
         getlev              = .true.
         relscale            = 2             ! m -> cm (or mm?)
!        select case (sfs)
!        case (106)          ! Soil water/ice content
!        end select
      case (105)
         select case (sfs)
!        case (101,255)
         case (1,101,255)
            gi%lv% leveltype = WMO3_HYBRID   ! 109
         case (105)
            gi%lv% leveltype = WMO3_HYBRIDB  ! 110
         case default
            write(0,*) 'param = ', trim (gi%pa% shortname), &
                 gi%pa% discipline, gi%pa% category, gi%pa% number,ffs
            write(0,*) 'unknown typeOfSecondFixedSurface =', sfs
            call finish('get_inventory','unknown typeOfSecondFixedSurface')
         end select
      case (WMO3_GENV) !150
         ! COSMO/ICON generalized vertical coordinate
         gi%lv% leveltype    = WMO3_GENV     ! 150
         select case (sfs)
         case (101,255,WMO3_GENV)
         case default
            write(0,*) 'param = ', trim (gi%pa% shortname), &
                 gi%pa% discipline, gi%pa% category, gi%pa% number,ffs
            write(0,*) 'invalid typeOfSecondFixedSurface =', sfs
            call finish('get_inventory','invalid typeOfSecondFixedSurface')
         end select
      case (WMO3_LAKE_BOT,WMO3_SEDIM_BOT,WMO3_MIX_LAYER)
         gi%lv% leveltype    = ffs
      case (WMO3_CLD_BASE,WMO3_CLD_TOPS,WMO3_NOM_TOA, &
            WMO3_ZERO_ISO,WMO3_ISOTHERM               )
         write(0,*) 'not implemented: typeOfFirstFixedSurface =', ffs
         gi%lv% leveltype = ffs
      case (WMO3_SURF_HORZ,WMO3_SURF_TANG)  ! SYNOP verification of radiation
         gi%lv% leveltype = ffs
      case default
         write(0,*) 'not implemented: typeOfFirstFixedSurface =', ffs
         call message('get_inventory','unknown typeOfFirstFixedSurface')
         gi%lv% leveltype = ffs
      end select
      if (getlev) then
         call grib_get (grib% handle,'scaleFactorOfFirstFixedSurface' ,scf)
         call grib_get (grib% handle,'scaledValueOfFirstFixedSurface' ,scv)
         gi%lv% levels(1)    = nint (scv * 10._dp**(relscale - scf))
         select case (sfs)
         case (100,102,103,105,106,WMO3_GENV)
            call grib_get (grib% handle,'scaleFactorOfSecondFixedSurface',scf)
            call grib_get (grib% handle,'scaledValueOfSecondFixedSurface',scv)
            gi%lv% levels(2) = nint (scv * 10._dp**(relscale - scf))
         case default
            gi%lv% levels(2) = -1       ! Second surface missing/not handled!
         end select
      else
         select case (gi%lv% leveltype)
         case (WMO3_SURFACE, WMO3_SEALEVEL, WMO3_SURF_HORZ, -1)
            ! No meaningful level values expected.
         case default
           call grib_get (grib% handle,'topLevel'      ,gi%lv% levels(1))
           call grib_get (grib% handle,'bottomLevel'   ,gi%lv% levels(2))
         end select
      end if
    end if
    gi%lv% levelvalue = gi%lv% levels(1)
    !--------------------------------------------------------------------
    ! For soil variables, derive levelvalue from mean of layer boundaries
    ! (rounded "up").
    !--------------------------------------------------------------------
    if (gi%lv% leveltype == WMO3_BELOWSUR .and. gi%lv% levels(2) > 0) then
       gi%lv% levelvalue = (gi%lv% levels(1) + gi%lv% levels(2) + 1) / 2
    end if
!if (gi%lv% leveltype == WMO3_BELOWSUR) then
!print *, "Belowsurf:"
!print *, "Levelvalue, Levels(:) =", gi%lv% levelvalue, gi%lv% levels(:)
!end if

    if (gi%lv% leveltype == WMO3_GENV) then
      call grib_get (grib% handle,'nlev',              gi%lv% nlev    )
      call grib_get (grib% handle,'numberOfVGridUsed', gi%lv% grid_num)
!#if !defined(__ibm__)
      call grib_get (grib% handle,'uuidOfVGrid',       gi%lv% uuid    )
!#else /* GRIB_API version 1.11.0 not available */
!     gi%lv% uuid     = achar (0)
!#endif
    else
      gi%lv% nlev     = -1
      gi%lv% grid_num = -1
      gi%lv% uuid     = achar (0)
    end if

    if (gi%ds% edition == 2) then
      !------------------------------------------------
      ! Look up entry in DWD table of "local concepts",
      ! derive corresponding GRIB1 code, table.
      !------------------------------------------------
      code2 = search_grib2 (name   =gi%pa% shortname,  &
                            dis    =gi%pa% discipline, &
                            cat    =gi%pa% category,   &
                            num    =gi%pa% number,     &
                            levtyp1=ffs)
!print *, "dis,num,cat,ffs =", gi%pa% discipline, gi%pa% category, gi%pa% number,ffs
!if (code2% paramID > 0) print *, "OK1", code2% paramID
      !------------------------
      ! Last resort: any levtyp
      !------------------------
      if (code2% paramID < 0) then
        code2 = search_grib2 (name   =gi%pa% shortname,  &
                              dis    =gi%pa% discipline, &
                              cat    =gi%pa% category,   &
                              num    =gi%pa% number      )
!if (code2% paramID > 0) print *, "OK2", code2% paramID
      end if
      !--------------------------
      ! Despair: ignore shortname
      !--------------------------
      if (code2% paramID < 0) then
        code2 = search_grib2 (dis    =gi%pa% discipline, &
                              cat    =gi%pa% category,   &
                              num    =gi%pa% number      )
      end if
      !------------------------------
      ! Resignation: ignore shortname
      !------------------------------
      if (code2% paramID < 0) then
        code2 = search_grib2 (dis    =gi%pa% discipline, &
                              cat    =gi%pa% category,   &
                              num    =gi%pa% number      )
      end if
      if (code2% paramID < 0) then
         call message ("inventory_2",&
           "parameter not in GRIB2 table:"//trim (gi%pa% shortname))
      else
         code1 = search_grib1 (paramID = code2% paramID)
         if (code1% paramID < 0) then
!           call message ("inventory_2",&
!                "parameter not in GRIB1 table:"//trim (gi%pa% shortname))
         else
            gi%pa% code  = code1% code
            gi%pa% table = code1% table
!           print *, "code, table =", gi%pa% code , gi%pa% table
         end if
      end if
    end if
    !-----------------------------
    ! derive name from table, code
    !-----------------------------
    entry             = search (  code= gi%pa% code,    &
                                 table= gi%pa% table,   &
                                levtyp= gi%lv% leveltype)
    if (entry% name == '') &
      entry           = search (  code= gi%pa% code,    &
                                 table= gi%pa% table    )
    if ( gi%pa% shortname == ""                                .or.  &
        (gi%pa% shortname == "unknown" .and. entry% name /= "")    ) &
      gi%pa% shortname = entry% name
    gi%pa%   iname     = entry% iname

    gi%pa% factor      = 1._dp
    !-----------------------------------------------------------
    ! Save conversion factor (e.g. some ECMWF fields from mars)
    !-----------------------------------------------------------
    if (entry% fak /= 1._dp .and. entry% units /= "") then
       gi%pa% factor   = entry% fak
    end if
    !-----
    ! grid
    !-----

! Only for GME grid?
!   call grib_get (grib% handle,'totalNumberOfGridPoints',gi%gr% nxny, stat)
    call grib_get (grib% handle,'numberOfPoints'         ,gi%gr% nxny)

    gi%gr% gridtype = -1
    if (gi%ds% edition == 1) then
       call grib_get (grib% handle,'dataRepresentationType' , i, stat)
!      print *, "dataRepresentationType", i, stat
       if (stat == 0) then
          gi%gr% gridtype = i
       end if
    else
       call grib_get (grib% handle,'gridDefinitionTemplateNumber', i, stat)
!      print *,'gridDefinitionTemplateNumber', i, stat
       if (stat == 0) then
          select case (i)
          case (0)
             gi%gr% gridtype = WMO6_LATLON
          case (1)
             gi%gr% gridtype = WMO6_ROTLL
          case (40)
             gi%gr% gridtype = WMO6_GAUSSIAN
          case (50)
             gi%gr% gridtype = WMO6_HARMONIC
          case (100)
             gi%gr% gridtype = DWD6_ICOSAHEDRON
          case (101)
             gi%gr% gridtype = DWD6_ICON
          case default
             write(0,*) 'unknown gridDefinitionTemplateNumber =',i
             call finish('get_inventory','invalid gridtype')
          end select
       end if
    end if

    gi%gr% ni       = 0
    gi%gr% grid_num = 0
    gi%gr% grid_ref = 0

    SELECT CASE (gi%gr% gridtype)
    CASE (DWD6_ICOSAHEDRON)
      call grib_get (grib% handle,'Ni', gi%gr% ni)
    CASE (WMO6_LATLON, WMO6_ROTLL)
      call grib_get (grib% handle,'Ni', gi%gr% ni)
    CASE (WMO6_GAUSSIAN)
      call grib_get (grib% handle,'Nj', gi%gr% ni)
      gi%gr% ni = gi%gr% ni * 2
    CASE (WMO6_HARMONIC)
      call grib_get (grib% handle,'J',  gi%gr% ni)
    CASE (DWD6_ICON)
      call grib_get (grib% handle,'numberOfGridUsed'       , gi%gr% grid_num)
      call grib_get (grib% handle,'numberOfGridInReference', gi%gr% grid_ref)
!     print *, "ICON: number of points:", gi%gr% nxny
      select case (gi%gr% grid_ref)
      case (1)  ! Triangle centers
         gi%gr% ni = nint (sqrt ( gi%gr% nxny    / 20._dp))
      case (2)  ! Triangle vertices
         gi%gr% ni = nint (sqrt ((gi%gr% nxny-2) / 10._dp))
      case (3)  ! Midpoints of triangle sides (edges)
         gi%gr% ni = nint (sqrt ( gi%gr% nxny    / 30._dp))
      case default
         gi%gr% ni = 0
         print *, "Ignoring numberOfGridInReference =", gi%gr% grid_ref
      end select
      uuid = achar (0)
#if !defined(__ibm__)
! Retrieving the UUID requires the proper template "template.3.101.def"
! and GRIB_API version >= 1.11.0
      call grib_get (grib% handle,'uuidOfHGrid', uuid)
#endif
      gi%gr% uuid = uuid
    CASE DEFAULT
      write(0,*) 'get_inventory: invalid gridtype =',gi%gr% gridtype
      call finish('get_inventory','invalid gridtype')
    END SELECT

    !---------------------------------------------
    ! local extensions: expid, runclass, ensembles
    !---------------------------------------------
    gi%pa% runclass =   0
    gi%pa% expid    =   0 ! experiment id
    gi%pa% zen      = 255 ! additional element number
    gi%en% id       =   0 ! Ensemble identification by table
    gi%en% size     =   0 ! Number of ensemble members
    gi%en% no       =  -1 ! Actual number of ensemble member
    if (local_dwd (gi)) then
      if (any (gi% ct% local_ident == [152,153,252,253])) then
        if (gi%ds% edition == 1) then
          call grib_get (grib% handle,'localEnsembleIdentification'      ,gi%en% id)
          call grib_get (grib% handle,'localActualNumberOfEnsembleNumber',gi%en% no)
          call grib_get (grib% handle,'localNumberOfEnsembleMembers'     ,gi%en% size)
        else
          call grib_get (grib% handle,'localTypeOfEnsembleForecast',gi%en% id)
          select case (pdtn)
          case default
           if (any (gi% ct% local_ident == [152,153])) then
            call grib_get(grib% handle,'localPerturbationNumber'   ,gi%en% no)
           else
            call grib_get(grib% handle,'perturbationNumber'        ,gi%en% no)
           end if
          case (2,12)
           call grib_get(grib% handle,'derivedForecast'            ,i)
           select case (i)
           case (0)
              gi%en% no = ENS_MEAN
           case (4)
              gi%en% no = ENS_SPREAD
           end select
          end select
          if (any (gi% ct% local_ident == [152,153])) then
           call grib_get(grib% handle,'localNumberOfForecastsInEnsemble',gi%en% size)
          else
           call grib_get(grib% handle,'numberOfForecastsInEnsemble'     ,gi%en% size)
          end if
        end if
      endif
      if (gi%ds% edition == 1) then
        call grib_get (grib% handle,'localElementNumber',gi%pa% zen)
        call grib_get (grib% handle,'localVersionNumber',gi%pa% expid)
        gi%pa% runclass =      gi%pa% expid / (2**14)
        gi%pa% expid    = mod (gi%pa% expid ,  2**14)
      else
        call grib_get (grib% handle,'productionStatusOfProcessedData',i)
        select case (i)
        case (0) ! Operational      (routine)
        case (1) ! Operational test (parallel suite)
        case (2) ! Research         (experiments)
        case default
        end select
        call grib_get (grib% handle,'localInformationNumber' ,gi%pa% zen)
        call grib_get (grib% handle,'localNumberOfExperiment',gi%pa% expid)
!       call grib_get (grib% handle,'backgroundGeneratingProcessIdentifier',i,stat)
        call grib_get (grib% handle,'backgroundProcess'      ,i,stat)
        if (stat == 0) then
          gi%pa% runclass = i
        else
          call grib_get (grib% handle,'localVersionNumber',i)
          gi%pa% runclass = i / (2**14)
        end if
      end if
    else if (local_ecmwf (gi)) then
      call grib_get (grib% handle,'marsType'                   ,mt)
      call grib_get (grib% handle,'marsStream'                 ,ms)
      call grib_get (grib% handle,'experimentVersionNumber'    ,c)
      if (verify (c,'0123456789 ') == 0) then
         read (c,*) gi%pa% expid   ! ASCII->int
!        print *, "ECMWF: experimentVersionNumber='",trim(c),"'",gi%pa% expid
      end if
      !--------------------------------------
      ! Cosmetic fixes for inventory printout
      !--------------------------------------
      select case (gi%ti% range)
      case (WMO5_FORECAST, DWD5_INIT_FC, WMO5_INIT_ANA, DWD5_IFS_FC)
        select case (mt)
        case (1)  ! First guess
          gi%pa% runtype = 'forecast'
        case (2)  ! Analysis
          gi%pa% runtype = 'analysis'
        case (3)  ! Intialized analysis
          gi%pa% runtype = 'init_ana'
        case (9,10,11)  ! 9=Forecast, 10=Control Forecast, 11=Perturbed Forecast
          if (gi%ti% p1 == 0) then
            gi%pa% runtype = 'init_ana'
          else
            gi%pa% runtype = 'forecast'
          end if
        end select
      end select

      select case (ms)
      case (1025)       ! Deterministic
      case (1030:1035)  ! Ensembles
!       call grib_get (grib% handle,'etyp???           '         ,gi%en% id)
        call grib_get (grib% handle,'perturbationNumber'         ,gi%en% no)
        call grib_get (grib% handle,'numberOfForecastsInEnsemble',gi%en% size)
      end select
    else if (gi%ds% edition == 2) then
      call grib_get (grib% handle,'backgroundProcess',i,stat)
      if (stat == 0) then
        gi%pa% runclass = i
      end if
    endif

    !-------------------------------------------------------------------
    ! Fix wrong coding of time range for selected ECMWF fields from mars
    !-------------------------------------------------------------------
    if (gi%ct% center == WMO0_ECMWF .and. gi%ds% edition == 1) then
       if (gi%pa% iname == 'tot_prec') then
          gi%ti% range = WMO5_ACCU
          gi%ti% p2 = max (gi%ti% p1, gi%ti% p2)
          gi%ti% p1 = 0
       end if
    end if
    !-----------------------------
    ! set time range for edition 1
    !-----------------------------
    if (gi%ds% edition == 1) then
      select case (gi%ti% range)
      case (WMO5_RANGE,WMO5_AVERAGE,WMO5_ACCU)
        gi%ti% rng_time = time_p12 (gi%ti% p2 - gi%ti% p1, gi%ti% range_unit)
      end select
    end if
#else
    call finish ('inventory_2','not linked with GRIB API')
#endif

  END FUNCTION inventory_2

!==============================================================================
  function local_dwd_invt (invt) result (local_dwd)
  type(t_inventory) ,intent(in) :: invt       ! GRIB inventory table entry
  logical                       :: local_dwd  ! local DWD extension present

    local_dwd = (invt% ct% center      == WMO0_DWD     .or.  &
                 invt% ct% center      == WMO0_COSMO   .or.  &
                 invt% ct% center      == WMO0_MSWISS  .or.  &
                 invt% ct% center      == WMO0_COMET ) .and. &
                (invt% ct% subcenter   == 0   .or.           &
                 invt% ct% subcenter   == 255     )    .and. &
                (invt% ct% local_ident == 152 .or.           &
                 invt% ct% local_ident == 153 .or.           &
                 invt% ct% local_ident >= 252 .and.          &
                 invt% ct% local_ident <= 254)

  end function local_dwd_invt
!------------------------------------------------------------------------------
  function local_ecmwf_invt (invt) result (local_ecmwf)
  type(t_inventory) ,intent(in) :: invt         ! GRIB inventory table entry
  logical                       :: local_ecmwf  ! local ECMWF extension present

    local_ecmwf = invt% ct% center      == WMO0_ECMWF .and. &
!                 invt% ct% sub_center  == xxx .and. &
                  invt% ct% local_ident >    0

  end function local_ecmwf_invt
!==============================================================================

  FUNCTION ref_time (grib)
  !-------------------------------------------------------------
  ! derive reference (start) time from GRIB record or API handle
  !-------------------------------------------------------------
  TYPE (t_time)  :: ref_time  ! reference time
  TYPE (t_grib1) :: grib      ! decoded GRIBEX record or GRIB API handle

#ifdef GRIB_API
    integer :: yyyymmdd, hhmm
#endif
    !--------------------------
    ! decode from GRIBEX record
    !--------------------------
    if (grib% handle == INVALID_HANDLE) then
      CALL init_time (ref_time, yyyy = grib% isec1% year              &
                                     +(grib% isec1% century-1) * 100, &
                                mo   = grib% isec1% month,            &
                                dd   = grib% isec1% day,              &
                                hh   = grib% isec1% hour,             &
                                mi   = grib% isec1% minute            )
    !----------------------------
    ! alternatively use GRIB2-API
    !----------------------------
    else
#ifdef GRIB_API
      call grib_get (grib% handle, 'dataDate',yyyymmdd)
      call grib_get (grib% handle, 'dataTime',hhmm    )
      call init_time (ref_time, yyyymmdd = yyyymmdd, &
                                hhmmss   = hhmm*100  )
#else
      call finish ('ref_time','not linked with GRIB2 API')
#endif
    endif

  END FUNCTION ref_time

!------------------------------------------------------------------------------

  FUNCTION db_time (grib)
  !-----------------------------------------------------
  ! derive data base time from GRIB record or API handle
  !-----------------------------------------------------
  TYPE (t_time)  :: db_time  ! data base time
  TYPE (t_grib1) :: grib     ! decoded GRIBEX record or GRIB API handle

#ifdef GRIB_API
    integer :: edition, yyyy, mo, dd, hh, mi, ss, localdef, status
#endif
    !--------------------------
    ! set default invalid value
    !--------------------------
    db_time = ZERO_TIME
    !--------------------------
    ! decode from GRIBEX record
    !--------------------------
    if (grib% handle == INVALID_HANDLE) then
      if (local_dwd (grib)) then
        CALL init_time (db_time, yyyy = grib% s1_dwd% year + 1900, &
                                 mo   = grib% s1_dwd% month,       &
                                 dd   = grib% s1_dwd% day,         &
                                 hh   = grib% s1_dwd% hour,        &
                                 mi   = grib% s1_dwd% minute       )
      endif
    !----------------------------
    ! alternatively use GRIB2-API
    !----------------------------
    else
#ifdef GRIB_API
      call grib_get (grib% handle, 'editionNumber'        ,edition)
      call grib_get (grib% handle, 'localDefinitionNumber',localdef, status)
      if (status /= GRIB_SUCCESS) localdef = 0
      select case (edition)
      case (1)
        select case (localdef)
        case (252,253,254) ! DWD
          call grib_get (grib% handle, 'localDecodeDateYear'  ,yyyy)
          call grib_get (grib% handle, 'localDecodeDateMonth' ,mo  )
          call grib_get (grib% handle, 'localDecodeDateDay'   ,dd  )
          call grib_get (grib% handle, 'localDecodeDateHour'  ,hh  )
          call grib_get (grib% handle, 'localDecodeDateMinute',mi  )
          yyyy = yyyy + 1900
          CALL init_time (db_time, yyyy = yyyy, &
                                   mo   = mo,   &
                                   dd   = dd,   &
                                   hh   = hh,   &
                                   mi   = mi    )
!      case (1,30)         ! ECWMF
!      case default
       end select
      case (2)
        select case (localdef)
        case (152,153,252:254) ! DWD
          call grib_get (grib% handle, 'localCreationDateYear'  ,yyyy)
          call grib_get (grib% handle, 'localCreationDateMonth' ,mo  )
          call grib_get (grib% handle, 'localCreationDateDay'   ,dd  )
          call grib_get (grib% handle, 'localCreationDateHour'  ,hh  )
          call grib_get (grib% handle, 'localCreationDateMinute',mi  )
          call grib_get (grib% handle, 'localCreationDateSecond',ss  )
          CALL init_time (db_time, yyyy = yyyy, &
                                   mo   = mo,   &
                                   dd   = dd,   &
                                   hh   = hh,   &
                                   mi   = mi,   &
                                   ss   = ss    )
!       case (1)           ! ECWMF
!       case default
       end select
      case default
        call finish ('inventory_2','invalid GRIB edition')
      end select
#else
      call finish ('db_time','not linked with GRIB2 API')
#endif
    endif

  END FUNCTION db_time

!------------------------------------------------------------------------------

  FUNCTION time_p12 (p12, unit)
  !-------------------------------------------
  ! derive time range from GRIB record entries
  !-------------------------------------------
  TYPE (t_time)       :: time_p12  ! time
  integer ,intent(in) :: p12       ! time range
  integer ,intent(in) :: unit      ! time range unit

    !----------------------------------
    ! time range, average, accumulation
    !----------------------------------
    SELECT CASE (unit)
    CASE (WMO4_SECOND)
                          CALL init_time (time_p12, ss = p12)
    CASE (WMO4_MINUTE)
                          CALL init_time (time_p12, mi = p12)
    CASE (DWD4_15_MINUTES)
                          CALL init_time (time_p12, mi = p12 * 15)
    CASE (WMO4_HOUR)
                          CALL init_time (time_p12, hh = p12)
    CASE (WMO4_12_HOURS)
                          CALL init_time (time_p12, hh = p12 * 12)
    CASE (WMO4_DAY)
                          CALL init_time (time_p12, dd = p12)
    CASE (WMO4_MONTH)
                          CALL init_time (time_p12, mo = p12)
!   case (WMO4_YEAR)
!   case (WMO4_DECADE)
!   case (WMO4_NORMAL)
!   case (WMO4_CENTURY)
    CASE (WMO4_3_HOURS)
                          CALL init_time (time_p12, hh = p12 * 3)
!   case (WMO4_6_HOURS)
    CASE default
      write(0,*) 'time_p12: time unit',unit,'not yet implemented'
      CALL finish ('time_p12','time unit not yet implemented')
    END SELECT

  end FUNCTION time_p12
!------------------------------------------------------------------------------

  FUNCTION ver_time (grib)
  !--------------------------------------------------------
  ! derive verification time from GRIB record or API handle
  !--------------------------------------------------------
  TYPE (t_time)  :: ver_time  ! verification time
  TYPE (t_grib1) :: grib      ! GRIB record or API handle

    TYPE (t_time)    :: p1
    CHARACTER(len=3) :: c_range
#ifdef GRIB_API
    integer          :: yyyymmdd, hhmm
#endif
    !--------------------
    ! decode GRIB1 record
    !--------------------
    if (grib% handle == INVALID_HANDLE) then
      CALL init_time (ver_time, yyyy = grib% isec1% year              &
                                     +(grib% isec1% century-1) * 100, &
                                mo   = grib% isec1% month,            &
                                dd   = grib% isec1% day,              &
                                hh   = grib% isec1% hour,             &
                                mi   = grib% isec1% minute            )

      SELECT CASE (grib% isec1% time_range)
      CASE (WMO5_FORECAST, WMO5_INIT_ANA, DWD5_INIT_FC, DWD5_NUDGING, &
            DWD5_IFS_FC, DWD5_ENKF_ANA                                )
        !--------------------------------
        ! forecast or analysis or nudging
        !--------------------------------
        SELECT CASE (grib% isec1% time_unit)
        CASE (WMO4_SECOND)
                              CALL init_time (p1, ss = grib% isec1% p1)
        CASE (WMO4_MINUTE)
                              CALL init_time (p1, mi = grib% isec1% p1)
        CASE (DWD4_15_MINUTES)
                              CALL init_time (p1, mi = grib% isec1% p1 * 15)
        CASE (WMO4_HOUR)
                              CALL init_time (p1, hh = grib% isec1% p1)
        CASE (WMO4_3_HOURS)
                              CALL init_time (p1, hh = grib% isec1% p1 * 3)
        CASE (WMO4_12_HOURS)
                              CALL init_time (p1, hh = grib% isec1% p1 * 12)
        CASE (WMO4_DAY)
                              CALL init_time (p1, dd = grib% isec1% p1)
        CASE (WMO4_MONTH)
                              CALL init_time (p1, mo = grib% isec1% p1)
!       case (WMO4_YEAR)
!       case (WMO4_DECADE)
!       case (WMO4_NORMAL)
!       case (WMO4_CENTURY)
!       case (WMO4_6_HOURS)
        CASE default
          write(*,*) "time_range,time_unit,p1 =", grib% isec1% time_range, &
               grib% isec1% time_unit, grib% isec1% p1
          CALL finish ('ver_time','time unit not yet implemented')
        END SELECT
        ver_time = ver_time + p1
      CASE (WMO5_RANGE,WMO5_AVERAGE,WMO5_ACCU)
        !----------------------------------
        ! time range, average, accumulation
        !----------------------------------
        ver_time = ver_time + time_p12 (grib% isec1% p2, grib% isec1% time_unit)
      CASE (WMO5_VALID)
        !--------------------------------------------
        ! Valid at ref.time+P1 (P1 occupies 2 octets)
        !--------------------------------------------
        SELECT CASE (grib% isec1% time_unit)
        CASE (WMO4_MINUTE)
                              CALL init_time (p1, mi = grib% isec1% p2 &
                                               + 256 * grib% isec1% p1 )
        CASE (DWD4_15_MINUTES)
                              CALL init_time (p1, mi =(grib% isec1% p2 &
                                               + 256 * grib% isec1% p1 ) * 15)
        CASE (WMO4_HOUR)
                              CALL init_time (p1, hh = grib% isec1% p2 &
                                               + 256 * grib% isec1% p1 )
        CASE (WMO4_DAY)
                              CALL init_time (p1, dd = grib% isec1% p2 &
                                               + 256 * grib% isec1% p1 )
        CASE (WMO4_MONTH)
                              CALL init_time (p1, mo = grib% isec1% p2 &
                                               + 256 * grib% isec1% p1 )
        CASE default
          CALL finish ('ver_time','time unit not yet implemented')
        END SELECT
        ver_time = ver_time + p1
      CASE default
        !--------------------
        ! not yet implemented
        !--------------------
        WRITE (c_range,'(i3)') grib% isec1% time_range
        CALL init_time (ver_time, yyyy=1,mo=1,dd=1)
      END SELECT
    !----------------------------
    ! alternatively use GRIB2-API
    !----------------------------
    else
#ifdef GRIB_API
      call grib_get (grib% handle, 'validityDate',yyyymmdd)
      call grib_get (grib% handle, 'validityTime',hhmm)
      call init_time (ver_time, yyyymmdd = yyyymmdd, &
                                hhmmss   = hhmm*100  )
#else
      call finish ('ver_time','not linked with GRIB2 API')
#endif
    endif

  END FUNCTION ver_time

!==============================================================================

  elemental FUNCTION equal_pa (pa1, pa2)
  LOGICAL                 :: equal_pa
  TYPE (t_pa) ,INTENT(in) :: pa1, pa2

    LOGICAL :: eq1, eq2

    ! GRIB1: compare table, code (if valid)
    eq1 = pa1% table /= 255        .AND. &
          pa2% table /= 255        .AND. &
          pa1% table == pa2% table .AND. &
          pa1% code  == pa2% code

    ! GRIB2: compare discipline, parameterCategory, parameterNumber
    eq2 = pa1% discipline /= 255             .AND. &
          pa2% discipline /= 255             .AND. &
          pa1% discipline == pa2% discipline .AND. &
          pa1% category   == pa2% category   .AND. &
          pa1% number     == pa2% number

    equal_pa = (eq1 .OR. eq2)                   &
         .AND. pa1% shortname == pa2% shortname &
         .AND. pa1% iname     == pa2% iname     &
         .AND. pa1% runtype   == pa2% runtype

  END FUNCTION equal_pa

!------------------------------------------------------------------------------

  elemental FUNCTION equal_ti (ti1, ti2)
  LOGICAL                 :: equal_ti
  TYPE (t_ti) ,INTENT(in) :: ti1, ti2

    equal_ti = ti1% ver_time == ti2% ver_time &
         .AND. ti1% p1       == ti2% p1       &
         .AND. ti1% p2       == ti2% p2       &
         .AND. ti1% range    == ti2% range

  END FUNCTION equal_ti

!------------------------------------------------------------------------------

  elemental FUNCTION equal_lv (lv1, lv2)
  LOGICAL                 :: equal_lv
  TYPE (t_lv) ,INTENT(in) :: lv1, lv2

    equal_lv =     lv1% leveltype  == lv2% leveltype  &
         .AND.     lv1% levelvalue == lv2% levelvalue &
         .AND. ALL(lv1% levels     == lv2% levels)

  END FUNCTION equal_lv

!------------------------------------------------------------------------------

  elemental FUNCTION equal_en (en1, en2)
  LOGICAL                 :: equal_en
  TYPE (t_en) ,INTENT(in) :: en1, en2

    equal_en =     en1% id   == en2% id   &
         .AND.     en1% size == en2% size &
         .AND.     en1% no   == en2% no

  END FUNCTION equal_en

!------------------------------------------------------------------------------

  elemental function equal_invt (i1, i2)
  type(t_inventory) ,intent(in) :: i1, i2    ! inventories to compare
  logical                       :: equal_invt ! result

    equal_invt =                         &!
!                (i1% ct == i2% ct)  &! identification of the center
                 (i1% pa == i2% pa)  &! characteristics of the parameter
           .and. (i1% ti == i2% ti)  &! time information
           .and. (i1% lv == i2% lv)  &! level information
!          .and. (i1% gr == i2% gr)  &! horizontal grid
           .and. (i1% en == i2% en)   ! ensemble data

  end function equal_invt

!==============================================================================
end module mo_grib_invt
