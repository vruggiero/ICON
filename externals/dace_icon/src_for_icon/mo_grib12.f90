!
!+ conversion of GRIB1 / GRIB2 messages using the GRIB API
!
MODULE mo_grib12
!
! Description:
!   Conversion of GRIB1 / GRIB2 messages using the GRIP API
!   Method:
!    A) flag 'grib_library' == 2 :
!       use the GRIB API to transparently encode/decode GRIB 1/2
!    B) flag 'grib_library' == 1 :
!       1) Use the GRIB API to convert the message
!       2) Afterwards fix explicitely entries not converted correctly,
!          i.e.: production status, experiment number,
!                icosahedral grid
!                dwd local extensions
!                data base time
!                time range
!                .....
!       3) Use the (C)GRIBEX library do decode/encode the message
!
!   Implementation in the 3dvar code:
!     Provide generic routines to read+decode or encode+write
!     GRIB 1 or GRIB 2 messages/files
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
!  precautions for GRIB2 API
! V1_22        2013-02-13 Harald Anlauf
!  Move grib_library,grib_edition handling to mo_grib*: set_grib
! V1_23        2013-03-26 Harald Anlauf
!  GRIB_API implementation
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON; FTRACE instrumentation.
! V1_27        2013-11-08 Harald Anlauf
!  Fixes for generalized vertical coordinate
!  GRIB2: properly distinguish init_ana,init_fc
! V1_45        2015-12-15 Harald Anlauf
!  get_inventory: don't crash on IFS spherical harmonics
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!------------------------------------------------------------------------------
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!------------------------------------------------------------------------------

  !-------------
  ! Modules used
  !-------------
  use mo_kind,        only: dp, wp                 ! double precision kind
  use mo_exception,   only: finish                 ! abort on error condition
  use mo_emos_grib1,  only: t_grib1,              &! GRIB record derived type
                            local_dwd,            &! check local DWD extension
                            pbgrib,               &! read message from file
                            gribex,               &! en/decode message
                            reallocate_buffer,    &! (re)allocate % kgrib
                            reallocate_data,      &! (re)allocate % rsec4
                            INVALID_HANDLE         ! invalid GRIB_API handle
  use mo_grib_invt,   only: t_inventory,          &! file inventory entry
                            inventory_entry,      &! derive entry from record
                            grib_library           ! GRIB API to use
  use mo_t_grib_api                                ! GRIB 2 Tables
  use mo_wmo_tables                                ! GRIB 1 Tables
#ifdef GRIB_API
  use mo_time,        only: i_time                 ! derive integers from time
  use grib_api,       only: grib_new_from_message,&! GRIB handle from message
                            grib_clone,           &! clone         GRIB message
                            grib_get_message_size,&! get size of   GRIB message
                            grib_copy_message,    &! get message from handle
                            grib_get,             &! get item from GRIB message
                            grib_set,             &! set item in   GRIB message
                            grib_release,         &! release       GRIB message
                            grib_get_size,        &! get size of an array item
                            kindOfSize,           &! kind parameter
                            GRIB_SUCCESS           ! return status value
#endif
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: grib1togrib2  ! convert GRIB 1 message to GRIB 2
  public :: grib2togrib1  ! convert GRIB 2 message to GRIB 1
  public :: read_gribex   ! read GRIB 1/2 message from file and decode


!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  subroutine grib1togrib2 (grib1, igrib2)
  !-----------------------------------
  ! convert a GRIB record to edition 2
  !-----------------------------------
  type(t_grib1) ,intent(in)  ::  grib1 ! GRIB1 derived type
  integer       ,intent(out) :: igrib2 ! GRIB2 handle

#ifndef GRIB_API
    call finish('grib1togrib2','not linked with GRIB API !')
#else
    !----------------
    ! Local variables
    !----------------
    integer           :: expid                     ! experiment Id
    integer           :: yyyy, mm, dd, hh, mi, ss  ! time as integers
    type(t_inventory) :: gi                        ! GRIB file inventory

    !========================================================
    ! get GRIB 2 handle, convert to GRIB 1 using the GRIB API
    !========================================================
    call grib_new_from_message (igrib2, grib1% kgrib)
    call grib_set (igrib2, 'editionNumber', 2)

    !=============================
    ! manually set/correct entries
    !=============================

    !------------------------------------------
    ! derive inventory entry from GRIB 1 record
    !------------------------------------------
    gi = inventory_entry (grib1)

    !-------------------------------------
    ! production status, experiment number
    !-------------------------------------
    expid = gi%pa% expid
    select case (expid)
    case (0:50)
      !----------------------
      ! Routine (operational)
      !----------------------
      call grib_set (igrib2, 'productionStatusOfProcessedData', 0)
      call grib_set (igrib2, 'localNumberOfExperiment',         0)
    case (51:99)
      !---------------------------------
      ! Parallel Suite (pre-operational)
      !---------------------------------
      call grib_set (igrib2, 'productionStatusOfProcessedData', 1)
      call grib_set (igrib2, 'localNumberOfExperiment',  50-expid)
    case (100:)
      !-----------
      ! Experiment
      !-----------
      call grib_set (igrib2, 'productionStatusOfProcessedData', 2)  ! changes short name !!!!!
      call grib_set (igrib2, 'localNumberOfExperiment',     expid)
    case default
      write(0,*) 'invalid experiment ID =',expid
      call finish('grib1togrib2','invalid experiment ID')
    end select

    !--------------------
    ! set time range, etc
    !--------------------
    select case (gi%pa% runtype)
    case ('analysis')
      call grib_set (igrib2, 'typeOfProcessedData',         GRIB_1_4_ANALYSIS)
      call grib_set (igrib2, 'significanceOfReferenceTime', GRIB_1_2_ANALYSIS)
    case ('forecast')
!   case ('init_ana')
!   case ('range   ')
!   case ('average ')
!   case ('accum   ')
!   case ('nudging ')
!   case ('ifs_fc  ')
!   case ('enkf_ana')
    case default
      call finish ('grib1togrib2', 'invalid runtype: '//gi%pa% runtype)
    end select

    !-----
    ! Grid
    !-----
    select case (grib1% isec2(1))
!   case (WMO6_LATLON)
!   case (WMO6_GAUSSIAN)
!   case (WMO6_ROTLL)
    case (DWD6_ICOSAHEDRON)
      call grib_set (igrib2,                                             &
!                   'latitudeOfThePolePointOfTheIcosahedronOnTheSphere', &
                    'latitudeOfThePolePoint',                            &
                     lat2 (grib1% tri% lat_pole)                         )
      call grib_set (igrib2, &
!                   'longitudeOfThePolePointOfTheIcosahedronOnTheSphere',&
                    'longitudeOfThePolePoint',                           &
                     lon2 (grib1% tri% lon_pole)                         )
      call grib_set (igrib2,                                             &
!                   'longitudeOfTheCenterLineOfTheFirstDiamondOfTheIcosahedronOnTheSphere',&
                    'longitudeOfFirstDiamondCenterLine',                 &
                     lon2 (grib1% tri% lon_dia1                         ))
    case default
    end select

! not handlesd so far:
! exponentOf2ForTheNumberOfIntervalsOnMainTriangleSides
! exponentOf3ForTheNumberOfIntervalsOnMainTriangleSides
! numberOfIntervalsOnMainTriangleSidesOfTheIcosahedron
! numberOfDiamonds
! gridPointPosition ('3.8.table',masterDir,localDir);
! numberingOrderOfDiamonds 'grib2/tables/[tablesVersion]/3.9.table';
! scanningModeForOneDiamond 'grib2/tables/[tablesVersion]/3.10.table';
! totalNumberOfGridPoints  : dump ;

!
!        repr        !  1 | 192  Data representation type (table 6).
!        ni2         !  2 | Number of factor 2 in factorisation of NI.
!        ni3         !  3 | Number of factor 3 in factorisation of NI.
!        nd          !  4 | Number of diamonds.
!        ni          !  5 | Number of triangular subdivisions.
!        orient      !  6 | Flag for orientation of diamonds.
!        lat_pole    !  7 | Latitude of pole point.
!        lon_pole    !  8 | Longitude of pole point.
!        lon_dia1    !  9 | Longitude of the first diamond.
!        scan_mode   ! 10 | Flag for storage sequence.
!        reserved_11 ! 11 |
!        nvcp        ! 12 ! Number of vertical coordinate parameters.

    !----------
    ! Parameter
    !----------
!   call grib_set (igrib2, 'shortName', gi% pa% name)

    !---------------------------
    ! handle DWD local extension
    !---------------------------
    if (local_dwd (grib1)) then
      yyyy =   grib1% s1_dwd% year + 1900
      if (yyyy < 1950) yyyy = yyyy + 100
      call grib_set (igrib2, 'localCreationDateYear',                  yyyy       )
      call grib_set (igrib2, 'localCreationDateMonth',  grib1% s1_dwd% month      )
      call grib_set (igrib2, 'localCreationDateDay',    grib1% s1_dwd% day        )
      call grib_set (igrib2, 'localCreationDateHour',   grib1% s1_dwd% hour       )
      call grib_set (igrib2, 'localCreationDateMinute', grib1% s1_dwd% minute     )
      call grib_set (igrib2, 'localCreationDateSecond', 0                         )
      call grib_set (igrib2, 'backgroundGeneratingProcessIdentifier', grib1% s1_dwd% run_type)
      call i_time   (gi% ti% ver_time, yyyy, mm, dd, hh, mi, ss)
      call grib_set (igrib2, 'localValidityDateYear',   yyyy)
      call grib_set (igrib2, 'localValidityDateMonth',  mm  )
      call grib_set (igrib2, 'localValidityDateDay',    dd  )
      call grib_set (igrib2, 'localValidityDateHour',   hh  )
      call grib_set (igrib2, 'localValidityDateMinute', mi  )
      call grib_set (igrib2, 'localValidityDateSecond', ss  )

!  not handled so far:
!  unsigned[1] localHostIdentifier: dump;
!  unsigned[2] localTypeOfEnsembleForecast: dump;
!  unsigned[2] localInformationNumber: dump;
!  unsigned[2] reserved;

!    unused values in GRIB 1:
!    INTEGER :: local_ident     ! 41    | local GRIB use identifier.       (254)
!    INTEGER :: day_number      ! 42    | not used any more                (255)
!    INTEGER :: record_number   ! 43-45 | not used any more                (200)
!    INTEGER :: decoding        ! 46    | not used any more                (255)
!    INTEGER :: element_no      ! 47    | not used any more                (255)

    endif

#endif

  end subroutine grib1togrib2
!==============================================================================
  function lon2 (lon1)
  !------------------------------------
  ! convert longitude: GRIB 1 -> GRIB 2
  !------------------------------------
  integer ,intent(in) :: lon1
  integer             :: lon2
    lon2 = lon1 * 1000
    if (lon2 < 0) lon2 = lon2 + 360000000
  end function lon2
!------------------------------------------------------------------------------
  function lat2 (lat1)
  !-----------------------------------
  ! convert latitude: GRIB 1 -> GRIB 2
  !-----------------------------------
  integer ,intent(in) :: lat1
  integer             :: lat2
    lat2 = lat1 * 1000
  end function lat2
!------------------------------------------------------------------------------
  function lon1 (lon2)
  !------------------------------------
  ! convert longitude: GRIB 2 -> GRIB 1
  !------------------------------------
  integer ,intent(in) :: lon2
  integer             :: lon1
    lon1 = (lon2+500) / 1000
  end function lon1
!------------------------------------------------------------------------------
  function lat1 (lat2)
  !-----------------------------------
  ! convert latitude: GRIB 2 -> GRIB 1
  !-----------------------------------
  integer ,intent(in) :: lat2
  integer             :: lat1
    lat1 = (lat2+500) / 1000
  end function lat1
!==============================================================================
  subroutine grib2togrib1 (grib1, hoper, x, kret, scanmode)
  !-----------------------------------
  ! convert a GRIB record to edition 1
  !-----------------------------------
  type(t_grib1) ,intent(inout)           :: grib1     ! GRIB1 derived type
  character     ,intent(in)              :: hoper     ! decoding parameter
  real(wp)      ,intent(inout) ,optional :: x (:,:,:) ! field to extract
  integer       ,intent(out)   ,optional :: kret      ! error return value
  integer       ,intent(in)    ,optional :: scanmode  ! scanning mode for x

#ifndef GRIB_API
    call finish('grib2togrib1','not linked with GRIB API !')
#else

    integer                :: igrib1      ! GRIB 1 clone to convert
    integer                :: igrib2      ! GRIB 2 handle
    integer(kindOfSize)    :: byte_size   ! size of message in byte
    integer                ::  int_size   ! size of message in int
    character ,allocatable :: message (:) ! message (temporary)
    integer                :: runtype     ! ana, fc, ..
    integer                :: status      ! return status
    integer                :: ffs, sfs    ! temporaries for fixup
    integer                :: scf, scv    !     of soil variables
    integer                :: topd        ! typeOfProcessedData
    integer                :: togp        ! typeOfGeneratingProcess

    !-----------------------------------------
    ! get GRIB API handle, clone GRIB 2 record
    !-----------------------------------------
    call grib_new_from_message (igrib2, grib1% kgrib, status=status)
    if (status/=GRIB_SUCCESS)                                    &
      call finish('grib2togrib1','ERROR in grib_new_from_message')

    call grib_clone (igrib2, igrib1, status=status)
    if (status/=GRIB_SUCCESS)                         &
      call finish('grib2togrib1','ERROR in grib_clone')

    !========================================
    ! 2) convert GRIB 2 to GRIB 1 by GRIB API
    !========================================
    call grib_set (igrib1, 'editionNumber', 1, status=status)
    if (status/=GRIB_SUCCESS)                       &
      call finish('grib2togrib1','ERROR in grib_set')

    !---------------------------------
    ! convert from GRIB API to CGRIBEX
    !---------------------------------
    call grib_get_message_size (igrib1, byte_size, status=status)
    if (status/=GRIB_SUCCESS)                                    &
      call finish('grib2togrib1','ERROR in grib_get_message_size')

    int_size = (byte_size + 3) / 4
    allocate (message (int_size*4))
    call reallocate_buffer     ( grib1, int_size )
    call grib_copy_message     (igrib1, message  )
    call grib_release          (igrib1)
    grib1% kgrib = transfer (message, grib1% kgrib)
    deallocate (message)

    !-----------------------------
    ! 3) decode  GRIB 1 by CGRIBEX
    !-----------------------------
!   call gribex (grib1, 'J')        ! Decode sections
    call gribex (grib1, 'D')        ! Decode, but keep reduced grid
!   call gribex (grib1, 'R')        ! Decode, expanding reduced grid

    !==========================================
    ! fix entries not converted properly so far
    !==========================================

    !----
    ! PDS
    !----
    select case (grib1% isec1% level_type)
    case (WMO3_BELOWSUR)
!     print *, "Below surface:", grib1% isec1% level_st, grib1% isec1% level_b
      call grib_get   (igrib2,'typeOfFirstFixedSurface'        ,ffs)
      call grib_get   (igrib2,'scaleFactorOfFirstFixedSurface' ,scf)
      call grib_get   (igrib2,'scaledValueOfFirstFixedSurface' ,scv)
      grib1% isec1%   level_st = nint (scv * 10._dp**(2-scf))   ! m -> cm
      call grib_get   (igrib2,'typeOfSecondFixedSurface'       ,sfs)
      if (sfs /= 255) then
        call grib_get (igrib2,'scaleFactorOfSecondFixedSurface',scf)
        call grib_get (igrib2,'scaledValueOfSecondFixedSurface',scv)
        grib1% isec1% level_b  = nint (scv * 10._dp**(2-scf))   ! m -> cm
      else
        grib1% isec1% level_b  = 0
      end if
    end select

    !-----
    ! Grid
    !-----
    select case (grib1% isec2(1))
    case (WMO6_LATLON,WMO6_ROTLL)
       ! In edition 2, all longitudes are between 0 and 360.
       if (grib1% latlon% lon_first > grib1% latlon% lon_last) then
          grib1% latlon% lon_first = grib1% latlon% lon_first - 360000
       end if
    case (WMO6_GAUSSIAN)
       ! In edition 2, all longitudes are between 0 and 360.
       if (grib1% gauss% lon_first > grib1% gauss% lon_last) then
          grib1% gauss% lon_first = grib1% gauss% lon_first - 360000
       end if
    case (DWD6_ICOSAHEDRON)
      call grib_get (igrib2,                                             &
!                   'latitudeOfThePolePointOfTheIcosahedronOnTheSphere', &
                    'latitudeOfThePolePoint',                            &
                     grib1% tri% lat_pole                                )
      call grib_get (igrib2, &
!                   'longitudeOfThePolePointOfTheIcosahedronOnTheSphere',&
                    'longitudeOfThePolePoint',                           &
                     grib1% tri% lon_pole                                )
      call grib_get (igrib2,                                             &
!                   'longitudeOfTheCenterLineOfTheFirstDiamondOfTheIcosahedronOnTheSphere',&
                    'longitudeOfFirstDiamondCenterLine',                 &
                     grib1% tri% lon_dia1                                )
      grib1% tri% lat_pole = lat1 (grib1% tri% lat_pole)
      grib1% tri% lon_pole = lon1 (grib1% tri% lon_pole)
      grib1% tri% lon_dia1 = lon1 (grib1% tri% lon_dia1)
    case default
    end select

    !--------------------
    ! set time range, etc
    !--------------------
    call grib_get (igrib2, 'significanceOfReferenceTime', runtype)
    select case (runtype)
    case (GRIB_1_4_ANALYSIS)
      grib1% isec1% time_range = WMO5_FORECAST
    case (GRIB_1_4_FORECAST)
      grib1% isec1% time_range = WMO5_FORECAST
    end select

    ! Code table 1.4
    call grib_get (igrib2,'typeOfProcessedData'    ,topd)
    call grib_get (igrib2,'typeOfGeneratingProcess',togp)
    select case (topd)
    case (0)   ! Analysis
       if (togp == 202) grib1% isec1% time_range = DWD5_NUDGING
    case (1)   ! Forecast
       if (grib1% isec1% p1 == 0) then
                        grib1% isec1% time_range = WMO5_INIT_ANA
       else
                        grib1% isec1% time_range = WMO5_FORECAST
         if (togp == 1) grib1% isec1% time_range = DWD5_INIT_FC
       end if
       if (togp == 195) grib1% isec1% time_range = DWD5_IFS_FC
    case (2)   ! Analysis and forecast
    case (3)   ! Control forecast
    case (4)   ! Perturbed forecast
    case (5)   ! Control and perturbed forecast
    case (192) ! Perturbed analysis
    case (255) ! Missing value
    case default
    end select

!  GRIB_1_4_ANALYSIS      =   0 ,&! Analysis products
!  GRIB_1_4_FORECAST      =   1 ,&! Forecast products
!  GRIB_1_4_ANA_FC        =   2 ,&! Analysis and forecast products
!  GRIB_1_4_CNTR_FC       =   3 ,&! Control forecast products
!  GRIB_1_4_PERT_FC       =   4 ,&! Perturbed forecast products
!  GRIB_1_4_CNTR_PERT_FC  =   5 ,&! Control and perturbed forecast products
!  GRIB_1_4_SAT_OBS       =   6 ,&! Processed satellite observations
!  GRIB_1_4_RADAR_OBS     =   7 ,&! Processed radar observations
!  GRIB_1_4_EVENT_PROB    =   8 ,&! Event Probability
!  GRIB_1_4_MISSING       = 255 ,&! Missing
!  EDZW1_4_PERT_ANA       = 192   ! Perturbed analysis

!   WMO5_FORECAST  =   0, &! Forecast valid for reference time + P1
!   WMO5_INIT_ANA  =   1, &! Initialized analysis (P1=0).
!   WMO5_RANGE     =   2, &! Range between ref.time+P1 and ref.time+P2
!   WMO5_AVERAGE   =   3, &! Average (ref.time+P1 to ref.time+P2)
!   WMO5_ACCU      =   4, &! Accumulation (ref.time+P1 to ref.time+P2)
!   WMO5_DIFF      =   5, &! Difference (ref.time+P1 minus ref.time+P2)
!   WMO5_VALID     =  10, &! Valid at ref.time+P1 (P1 occupies 2 octets)
!   DWD5_NUDGING   =  13, &! nudging analysis (DWD special)
!   DWD5_IFS_FC    =  14, &! IFS forecast
!   DWD5_ENKF_ANA  =  15, &! EnKF    analysis (DWD special)
!   WMO5_FC_AV_1   = 113, &! Average of N forecasts.
!   WMO5_FC_AC_1   = 114, &! Accumulation of N forecasts.
!   WMO5_FC_AV_R   = 115, &! Average of N forecasts, with the same ref.time.
!   WMO5_FC_AC_R   = 116, &! Accumulation of N forecasts, with same ref.time.
!   WMO5_FC_AV_2   = 117, &! Average of N forecasts.
!   WMO5_VAR_IANA  = 118, &! Temporal variance/covariance, of N init.analyses.
!   WMO5_STDEV_FC  = 119, &! Standard deviation of N forecasts.
!   WMO5_AV_IANA   = 123, &! Average of N uninitialized analyses.
!   WMO5_AC_IANA   = 124   ! Accumulation of N uninitialized analyses.

    !----------------
    ! local extension
    !----------------
    call grib_get (igrib2,'localDefinitionNumber',grib1% isec1% local_ident,&
                   status)
    if (status /= GRIB_SUCCESS) grib1% isec1% local_ident = 0
    select case (grib1% isec1% local_ident)
    case (254,253)
      !------------------------------------
      ! local extension flag and identifier
      !------------------------------------
      grib1% isec1%  local_flag    = 1
      grib1% s1_dwd% local_ident   = grib1% isec1% local_ident
      !---------------
      ! unused entries
      !---------------
      select case (runtype)
      case (GRIB_1_4_ANALYSIS)
        grib1% s1_dwd% day_number    = 0
        grib1% s1_dwd% record_number = 0
        grib1% s1_dwd% decoding      = 0
!       grib1% s1_dwd% element_no    = 0
      case default
        grib1% s1_dwd% day_number    = 255
        grib1% s1_dwd% record_number =  90
        grib1% s1_dwd% decoding      = 255
!       grib1% s1_dwd% element_no    = 255
      end select
      !---------------------------
      ! ZEN (Zusatz-Elementnummer)
      !---------------------------
      call grib_get (igrib2,'localInformationNumber',grib1% s1_dwd% element_no)
      !---------------
      ! data base time
      !---------------
      call grib_get (igrib2,'localCreationDateYear'                ,grib1% s1_dwd% year)
      call grib_get (igrib2,'localCreationDateMonth'               ,grib1% s1_dwd% month)
      call grib_get (igrib2,'localCreationDateDay'                 ,grib1% s1_dwd% day)
      call grib_get (igrib2,'localCreationDateHour'                ,grib1% s1_dwd% hour)
      call grib_get (igrib2,'localCreationDateMinute'              ,grib1% s1_dwd% minute)
      grib1% s1_dwd% year = grib1% s1_dwd% year - 1900
      !-----------------
      ! run type, exp_id
      !-----------------
      call grib_get (igrib2,'localNumberOfExperiment'              ,grib1% s1_dwd% exp)
      call grib_get (igrib2,'backgroundGeneratingProcessIdentifier',grib1% s1_dwd% run_type)
      call grib_get (igrib2,'productionStatusOfProcessedData'      ,               runtype)

      select case (runtype)
      case (0) ! Routine (operational)
        grib1% s1_dwd% exp = 1
      case (1) ! Parallel Suite (pre-operational)
        grib1% s1_dwd% exp = grib1% s1_dwd% exp + 50
      case (2) ! Experiment
      end select

      if (grib1% isec1% local_ident == 253) then
         call grib_get (igrib2,'localTypeOfEnsembleForecast',grib1% s1_dwd% ensemble_id)
         call grib_get (igrib2,'perturbationNumber'         ,grib1% s1_dwd% ensemble_no)
         call grib_get (igrib2,'numberOfForecastsInEnsemble',grib1% s1_dwd% ensemble_size)
      end if

!    not yet handled:
!    INTEGER :: user_id         ! 55    ! User id, specified by table
!    INTEGER :: experiment_id   ! 56-57 ! Experiment identifier
!    INTEGER :: major_version   ! 64    ! Model major version number
!    INTEGER :: minor_version   ! 65    ! Model minor version number

    end select

    !------------------------------------------------
    ! Representation: force bits per value to 16 if 0
    !------------------------------------------------
    if (grib1% isec4% bits == 0) then
      grib1% isec4% bits = 16
      call grib_get (igrib2,'getNumberOfValues',grib1% isec4% n_data)
      call grib_get (igrib2,'numberOfMissing'  ,grib1% isec4% non_miss)
      grib1% isec4% non_miss = grib1% isec4% n_data - grib1% isec4% non_miss
      call reallocate_data (grib1, grib1% isec4% non_miss)
      call grib_get (igrib2,'values'           ,grib1% rsec4)
    endif

    call grib_release (igrib2)

    if (present (x)) then
      call gribex (grib1, 'C')
      call gribex (grib1, hoper, x, kret, scanmode)
    endif

#endif

  end subroutine grib2togrib1
!==============================================================================
!     Provide generic routines to read+decode
!     GRIB 1 or GRIB 2 messages/files
!==============================================================================

  subroutine read_gribex (grib, hoper, x, kret, scanmode)
  type(t_grib1) ,intent(inout)           :: grib      ! GRIB derived type
  character     ,intent(in)              :: hoper     ! decoding parameter
  real(wp)      ,intent(inout) ,optional :: x (:,:,:) ! field to extract
  integer       ,intent(out)   ,optional :: kret      ! error return value
  integer       ,intent(in)    ,optional :: scanmode  ! scanning mode for x
  !----------------------------------------------
  ! read GRIB message and decode
  ! if edition number is 2
  ! the message will be converted by the GRIB API
  !----------------------------------------------
    character(len=4) :: edition
    integer          :: status, ed
    if (present (kret)) kret = 0
    !-----------------------------------------
    ! read message using the CGRIBEX interface
    !-----------------------------------------
    call pbgrib (grib, kret)
    if (present (kret)) then
      if (kret /= 0) return
    endif
    select case (grib_library)
    case (1)
      grib% handle = INVALID_HANDLE
      !-------------------------
      ! check for edition number
      !-------------------------
      edition = transfer (grib% kgrib(2), edition)
      ed      = ichar (edition(4:4))
      select case (ed)
      case (2)
        !------------------------------------
        ! GRIB 2 :convert GRIB 2 -> 1, decode
        !------------------------------------
        call grib2togrib1 (grib, hoper, x, kret, scanmode)
      case(1)
        !--------------------
        ! GRIB 1: just decode
        !--------------------
        if (present (x)) then
          call gribex (grib, hoper, x, kret, scanmode)
        else
          call gribex (grib, hoper, kret)
        endif
      case default
        write(6,*) 'read_gribex: invalid edition number =',ed
        call finish ('read_gribex','invalid edition number')
      end select
    case (2)
#ifdef GRIB_API
      if (grib% handle /= INVALID_HANDLE) call grib_release (grib% handle)
FTRACE_BEGIN("read_gribex:new_from_message")
      call grib_new_from_message (grib% handle, grib% kgrib, status)
      if (status /= GRIB_SUCCESS) &
        call finish('read_gribex','error from grib_new_from_message')
FTRACE_END  ("read_gribex:new_from_message")
      select case (hoper)
      case ('J') ! (GRIB1: Decode sections 0,1,2 and integer part of section 4)
        if (present (x)) then
          call finish('read_gribex',&
                      'no decoding of data for hoper = '//trim (hoper))
        end if
FTRACE_BEGIN("read_gribex:decode_head")
        call grib2_decode_head (grib)
FTRACE_END  ("read_gribex:decode_head")
!     case ('D') ! (GRIB1: Decode data from GRIB code)
      case ('R') ! (GRIB1: Decode data from GRIB code, incl. reduced grids)
        if (present (x)) then
FTRACE_BEGIN("read_gribex:decode_data")
          call grib2_decode_data (grib, x, scanmode)
FTRACE_END  ("read_gribex:decode_data")
        else
          call finish('read_gribex',"gribex('R') requires output field")
        end if
      case default
        call finish('read_gribex','not implemented: hoper = '//trim (hoper))
      end select
#else
      call finish('read_gribex','not linked with GRIB API')
#endif
    end select
  end subroutine read_gribex
!==============================================================================
#ifdef GRIB_API
  !-----------------------------------------------------------------
  ! Emulate decoding of essential GRIB1 sections for a GRIB2 message
  !-----------------------------------------------------------------
  subroutine grib2_decode_head (g)
    type(t_grib1) ,intent(inout) :: g      ! GRIB derived type

    integer           :: ed, gr, gridtype, h, j, nv, nt, ng, nr, ni, stat, pvp
    integer           :: sh, cp
    real(wp)          :: r(9)
!   character(len=1)  :: uuid(16)       ! Grid UUID
!   character(len=32) :: c

    h = g% handle
    !----------
    ! Section 0
    !----------
    call grib_get (h,'editionNumber', g% isec0% edition)
    call grib_get (h,'totalLength'  , g% isec0% n_octets)
    ed = g% isec0% edition
    !----------
    ! Section 1
    !----------
    call grib_get (h,'centre',                     g% isec1% center)
    call grib_get (h,'subCentre',                  g% isec1% sub_center)
    call grib_get (h,'generatingProcessIdentifier',g% isec1% process)
    call grib_get (h,'localDefinitionNumber'      ,g% isec1% local_ident,stat)
    if (stat /= GRIB_SUCCESS) g% isec1% local_ident = 0
    !----------
    ! Section 2
    !----------
    if (ed == 1) then
       call grib_get (h,'dataRepresentationType', gridtype)
!      call grib_get (h,'gridType', c)
!      select case (c)
!      case ('regular_ll')
!         gridtype = WMO6_LATLON
!      case ('rotated_ll')
!         gridtype = WMO6_ROTLL
!      case ('regular_gg')
!         gridtype = WMO6_GAUSSIAN
!      case ('triangular_grid')
!         gridtype = DWD6_ICOSAHEDRON
!      case default
!         call finish('grib2_decode_head','unknown gridtype: '//trim (c))
!      end select
    else
       call grib_get (h,'gridDefinitionTemplateNumber', gr)
       select case (gr)
       case (0)
          gridtype = WMO6_LATLON
       case (1)
          gridtype = WMO6_ROTLL
       case (40)
          gridtype = WMO6_GAUSSIAN
       case (100)
          gridtype = DWD6_ICOSAHEDRON
       case (101)
          gridtype = DWD6_ICON
       case default
          write(0,*) 'unknown gridDefinitionTemplateNumber =',gr
          call finish('grib2_decode_head','invalid gridtype')
       end select
    end if
    if (ed == 1) then
       call grib_get (h,'numberOfVerticalCoordinateValues', nv)
    else
       call grib_get (h,'NV', nv)
    end if

    select case (gridtype)
    case (WMO6_LATLON,WMO6_ROTLL)
       g% latlon% repr = gridtype
       g% latlon% nvcp = nv
       call grib_get (h,'Ni', g% latlon% ni)
       call grib_get (h,'Nj', g% latlon% nj)
       call grib_get (h,'latitudeOfFirstGridPointInDegrees' ,r(1))
       call grib_get (h,'longitudeOfFirstGridPointInDegrees',r(2))
       call grib_get (h,'latitudeOfLastGridPointInDegrees'  ,r(3))
       call grib_get (h,'longitudeOfLastGridPointInDegrees' ,r(4))
       call grib_get (h,'iDirectionIncrementInDegrees'      ,r(5))
       call grib_get (h,'jDirectionIncrementInDegrees'      ,r(6))
       g% latlon% lat_first = nint (r(1) * 1000)
       g% latlon% lon_first = nint (r(2) * 1000)
       g% latlon% lat_last  = nint (r(3) * 1000)
       g% latlon% lon_last  = nint (r(4) * 1000)
       g% latlon% di        = nint (r(5) * 1000)
       g% latlon% dj        = nint (r(6) * 1000)
       if (ed == 2) then
          ! In edition 2, all longitudes are between 0 and 360.
          if (g% latlon% lon_first > g% latlon% lon_last) then
             g% latlon% lon_first = g% latlon% lon_first - 360000
          end if
       end if
       call grib_get (h,'jScansPositively',j)
       g% latlon% scan_mode = j
       if (j == 1) g% latlon% scan_mode = WMO8_J_POSITIVE
       if (gridtype == WMO6_ROTLL) then
          call grib_get (h,'latitudeOfSouthernPoleInDegrees' ,r(7))
          call grib_get (h,'longitudeOfSouthernPoleInDegrees',r(8))
          call grib_get (h,'angleOfRotationInDegrees'        ,r(9))
          g% latlon% lat_rot   = nint (r(7) * 1000)
          g% latlon% lon_rot   = nint (r(8) * 1000)
          g% rsec2%  rot_angle =       r(9)
       end if
    case (WMO6_GAUSSIAN)
       g% gauss% repr = gridtype
       g% gauss% nvcp = nv
       call grib_get (h,'Ni', g% gauss% ni)
       call grib_get (h,'Nj', g% gauss% nj)
       call grib_get (h,'latitudeOfFirstGridPointInDegrees' ,r(1))
       call grib_get (h,'longitudeOfFirstGridPointInDegrees',r(2))
       call grib_get (h,'latitudeOfLastGridPointInDegrees'  ,r(3))
       call grib_get (h,'longitudeOfLastGridPointInDegrees' ,r(4))
       call grib_get (h,'iDirectionIncrementInDegrees'      ,r(5))
       g% gauss% lat_first = nint (r(1) * 1000)
       g% gauss% lon_first = nint (r(2) * 1000)
       g% gauss% lat_last  = nint (r(3) * 1000)
       g% gauss% lon_last  = nint (r(4) * 1000)
       if (abs (r(5)) <= 360._wp) then
         g% gauss% di      = nint (r(5) * 1000)
       else
         g% gauss% di      = 0
       end if
       call grib_get (h,'jScansPositively',j)
       g% gauss% scan_mode = j
       if (j == 1) g% gauss% scan_mode = WMO8_J_POSITIVE
    case (DWD6_ICOSAHEDRON)
       g% tri% repr = gridtype
       g% tri% nvcp = nv
       call grib_get (h,'n2', g% tri% ni2)
       call grib_get (h,'n3', g% tri% ni3)
       call grib_get (h,'nd', g% tri% nd)
       call grib_get (h,'Ni', g% tri% ni)
       if (ed == 1) then
         call grib_get (h,'latitudeOfIcosahedronPole'        ,g% tri% lat_pole)
         call grib_get (h,'longitudeOfIcosahedronPole'       ,g% tri% lon_pole)
         call grib_get (h,'longitudeOfFirstDiamondCentreLine',g% tri% lon_dia1,&
                                                              stat             )
         if (stat /= 0) &
         call grib_get (h,'longitudeOfFirstDiamondCenterLine',g% tri% lon_dia1)
!        call grib_get (h,'latitudeOfIcosahedronPoleInDegrees'        ,r(1))
!        call grib_get (h,'longitudeOfIcosahedronPoleInDegrees'       ,r(2))
!        call grib_get (h,'longitudeOfFirstDiamondCenterLineInDegrees',r(3))
       else
         call grib_get (h,'latitudeOfThePolePoint'           ,g% tri% lat_pole)
         call grib_get (h,'longitudeOfThePolePoint'          ,g% tri% lon_pole)
         call grib_get (h,'longitudeOfFirstDiamondCentreLine',g% tri% lon_dia1,&
                                                              stat             )
         if (stat /= 0) &
         call grib_get (h,'longitudeOfFirstDiamondCenterLine',g% tri% lon_dia1)
         g% tri% lat_pole = g% tri% lat_pole / 1000
         g% tri% lon_pole = g% tri% lon_pole / 1000
         g% tri% lon_dia1 = g% tri% lon_dia1 / 1000
         call grib_get (h,'totalNumberOfGridPoints', nt)
         if (nt /= 10*(g% tri% ni + 1)**2) then
            write(0,*) "Unexpected number of gridpoints:", &
                 nt,"/=", 10*(g%tri% ni + 1)**2
            call finish('grib2_decode_head','gridtype GME')
         end if
       end if
    case (DWD6_ICON)
       g% tri% repr = gridtype
       g% tri% nvcp = nv
       call grib_get (h,'numberOfDataPoints',      nt)
       call grib_get (h,'numberOfGridUsed',        ng)
       call grib_get (h,'numberOfGridInReference', nr)
       ni = 0
       select case (nr)
       case (1)  ! Triangle centers
          ni = nint (sqrt ( nt    / 20._dp))
       case (2)  ! Triangle vertices
          ni = nint (sqrt ((nt-2) / 10._dp))
       case (3)  ! Midpoints of triangle sides (edges)
          ni = nint (sqrt ( nt    / 30._dp))
       case default
          write(0,*) "unsupported numberOfGridInReference =", nr
          call finish('grib2_decode_head','gridtype ICON')
       end select
       g% tri% ni = ni
!      print *, "numberOfGridUsed =", ng
    case (WMO6_HARMONIC)
       g% sph% repr = gridtype
       g% sph% nvcp = nv
       call grib_get (h,'J', g% sph% j)
       call grib_get (h,'K', g% sph% k)
       call grib_get (h,'M', g% sph% m)
       call grib_get (h,'sphericalHarmonics',sh)
       call grib_get (h,'complexPacking'    ,cp)
       g% sph% repr_type = sh                     ! See code table 9
       g% sph% repr_mode = merge (1, 2, cp == 0)  ! See code table 10
!!! Unsupported:
!!     lat_rot     ! 13 | Latitude of the southern pole of rotation.
!!     lon_rot     ! 14 | Longitude of the southern pole of rotation.
!!     lat_strech  ! 15 | Latitude of the pole of stretching.
!!     lon_strech  ! 16 | Longitude of the pole of stretching.
    case default
       if (ed == 1) then
          write(0,*) 'unsupported dataRepresentationType =',gridtype
       else
          write(0,*) 'unknown gridDefinitionTemplateNumber =',gr
       end if
       call finish('grib2_decode_head','invalid gridtype')
    end select

    g% rsec2% vcp = -HUGE (0._dp)
    if (nv > 0) then
      call   grib_get (h,'PVPresent',pvp)       ! Check if PV really present
      if (pvp == 1) then
        call grib_get (h,'pv'       ,g% rsec2% vcp(1:nv))
      else
        nv = 0          ! PV not present (generalized vertical coord.?)
      end if
    end if

    g% isec2(1)  = gridtype
    g% isec2(12) = nv
!   print *, gridtype, nv

  end subroutine grib2_decode_head
!==============================================================================
  subroutine grib2_decode_data (grib, x, scanmode)
    type(t_grib1)     ,intent(inout) :: grib      ! GRIB derived type
    real(wp)          ,intent(inout) :: x (:,:,:) ! field to extract
    integer ,optional ,intent(in)    :: scanmode  ! scanning mode for x

    integer(kindOfSize) :: np
    integer             :: iscan, jscan, jcons, lb2, ub2, stat
    real(wp)            :: vals(size (x))
!   real(wp)            :: missingValue
    logical             :: yflip

!   call grib_get (grib% handle,'numberOfPoints',np)
    call grib_get_size (grib% handle,'values',np)

    if (np /= size (x)) then
       write(0,*) "Size mismatch while decoding data:", np, size (x)
       call finish ("grib2_decode_data","size mismatch")
    end if

    !-------------------------
    ! Check data storage order
    !-------------------------
    call grib_get (grib% handle,'iScansNegatively'     ,iscan, stat)
    if (stat /= 0) iscan = 0
    call grib_get (grib% handle,'jScansPositively'     ,jscan, stat)
    if (stat /= 0) jscan = -1
    call grib_get (grib% handle,'jPointsAreConsecutive',jcons, stat)
    if (stat /= 0) jcons = 0

    if (iscan /= 0) &
         call finish ("grib2_decode_data","iScansNegatively not implemented")
    if (jcons /= 0) &
         call finish ("grib2_decode_data","jPointsAreConsecutive not supported")
    !------------------------------------
    ! Check whether we need to flip N<->S
    ! (Beware of unstructured grids!)
    !------------------------------------
    yflip = .false.
    if (present (scanmode)) then
       if (iand (scanmode, WMO8_J_CONSECUTIVE) /= 0) &
         call finish ('grib2_decode_data','scanmode==WMO8_J_CONSECUTIVE')
       select case (scanmode)
       case (WMO8_J_POSITIVE)
          if (jScan == 0) yflip = .true.
       case (0)
          if (jScan >  0) yflip = .true.
       case default
         call finish ('grib2_decode_data','invalid scanmode')
       end select
    end if

!   missingValue = -9.e15_wp
!   call grib_set (grib% handle,'missingValue'  ,missingValue)

    call grib_get (grib% handle, 'values', vals)

    if (yflip) then
!      print *, "grib2_decode_data: flip!"
       lb2 = lbound (x,2)
       ub2 = ubound (x,2)
       x(:,ub2:lb2:-1,:) = reshape (vals, shape (x))
    else
       x(:,   :      ,:) = reshape (vals, shape (x))
    end if

!   print *, "grib_get: values =", minval (vals), "...", maxval (vals)

  end subroutine grib2_decode_data
#endif
!==============================================================================
end module mo_grib12
