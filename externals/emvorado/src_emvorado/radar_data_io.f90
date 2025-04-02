! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2005-2024, DWD, KIT
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

#ifdef _DACE_
#define __COSMO__
#endif

MODULE radar_data_io

!------------------------------------------------------------------------------
!
! Description: This module provides basic constants, variables and derived types for the
!              the various output methods, data and formats.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp
  USE radar_data, ONLY : nradsta_max, cmaxlen

!================================================================================
!================================================================================

  IMPLICIT NONE

!================================================================================
!================================================================================

  PRIVATE

  ! Unit numbers for control output to check contents of volume scan output files:
  !-------------------------------------------------------------------------------

  ! unit number of radar beam geometry control output file
  INTEGER, PUBLIC  :: radgeomoutputunit (nradsta_max)
  ! unit number of simul. radar radial wind control output file
  INTEGER, PUBLIC  :: radwindoutputunit (nradsta_max)
  ! unit number of observed radar radial wind control output file
  INTEGER, PUBLIC  :: radwindobsoutputunit (nradsta_max)
  ! unit number of simul. radar reflectivity control output file
  INTEGER, PUBLIC  :: radrefloutputunit (nradsta_max)
  ! unit number of observed radar reflectivity control output file
  INTEGER, PUBLIC  :: radreflobsoutputunit (nradsta_max)
  ! unit number of simul. ZDR control output file
  INTEGER, PUBLIC  :: zdroutputunit (nradsta_max)
  ! unit number of observed ZDR control output file
  INTEGER, PUBLIC  :: zdrobsoutputunit (nradsta_max)
  ! unit number of simul. RhoHV control output file
  INTEGER, PUBLIC  :: rhvoutputunit (nradsta_max)
  ! unit number of observed RhoHV control output file
  INTEGER, PUBLIC  :: rhvobsoutputunit (nradsta_max)
  ! unit number of simul. KDP control output file
  INTEGER, PUBLIC  :: kdpoutputunit (nradsta_max)
  ! unit number of observed KDP control output file
  INTEGER, PUBLIC  :: kdpobsoutputunit (nradsta_max)
  ! unit number of simul. radar extinction control output file
  INTEGER, PUBLIC  :: ahoutputunit (nradsta_max)
  ! unit number of simul. ADP control output file
  INTEGER, PUBLIC  :: adpoutputunit (nradsta_max)
  ! unit number of simul. LDR control output file
  INTEGER, PUBLIC  :: ldroutputunit (nradsta_max)
  ! unit number of observed LDR control output file
  INTEGER, PUBLIC  :: ldrobsoutputunit (nradsta_max)


  INTEGER, PARAMETER, PUBLIC  :: MAX_LEN_FILENAME  = cmaxlen
  INTEGER, PARAMETER, PUBLIC  :: LEN_DATETIME      = 14
  INTEGER, PARAMETER, PUBLIC  :: MAX_LEN_SNIPPET   = 20
  INTEGER, PARAMETER, PUBLIC  :: MAX_LEN_TAB_ENTRY = 200
  INTEGER, PARAMETER, PUBLIC  :: MAX_LEN_ATTR_NAME = 64
  INTEGER, PARAMETER, PUBLIC  :: MAX_LEN_ATTR      = cmaxlen
  INTEGER, PARAMETER, PUBLIC  :: MAX_LEN_DESC      = 800
  INTEGER, PARAMETER, PUBLIC  :: NATTR_MAX         = 25

  INTEGER, PARAMETER, PUBLIC  :: GRIB_MISSING   = 255
  INTEGER, PARAMETER, PUBLIC  :: GRIB_TRUE      = 1
  INTEGER, PARAMETER, PUBLIC  :: GRIB_FALSE     = 0
  INTEGER, PARAMETER, PUBLIC  :: GRIB_EDITION_2 = 2

  INTEGER, PARAMETER, PUBLIC  :: PRODUCT_TYPE_OBS = 1
  INTEGER, PARAMETER, PUBLIC  :: PRODUCT_TYPE_SIM = 2

  ! Type of generating process for grib2:
  INTEGER,  PARAMETER, PUBLIC :: TGP_ANALYSIS = 0
  INTEGER,  PARAMETER, PUBLIC :: TGP_INITRUN  = 1
  INTEGER,  PARAMETER, PUBLIC :: TGP_FORECAST = 2
  INTEGER,  PARAMETER, PUBLIC :: TGP_ENSEMBLE = 4

  ! Which value should be used in a grib2 bitmap:
  INTEGER,  PARAMETER, PUBLIC :: BITMAP_PRI_NOBITMAP        = 0    ! no bitmap at all
  INTEGER,  PARAMETER, PUBLIC :: BITMAP_PRI_MISSVAL         = 1    ! miss_value
  INTEGER,  PARAMETER, PUBLIC :: BITMAP_PRI_UNDEVAL         = 2    ! zero_value
  INTEGER,  PARAMETER, PUBLIC :: BITMAP_PRI_MISSVAL_UNDEVAL = 12   ! miss_value if miss_values present, otherwise zero_value

  ! Derived types for organizing output of volume scans in radar_output_methods.f90:
  !---------------------------------------------------------------------------------

  !
  ! ecCodes-Größen:
  !
  TYPE t_scaledMetadataList
    INTEGER                          :: codeFigure
    CHARACTER(LEN=MAX_LEN_TAB_ENTRY) :: meaning
    REAL(dp)                         :: value
  END TYPE t_scaledMetadataList
  !
  TYPE t_grib_sec_0
    INTEGER                     :: discipline
    INTEGER                     :: editionNumber
  END type t_grib_sec_0
  !
  TYPE t_grib_sec_1
    INTEGER                     :: tablesVersion
    INTEGER                     :: productionStatusOfProcessedData
    INTEGER                     :: typeOfProcessedData
    INTEGER                     :: centre
    INTEGER                     :: subCentre
    INTEGER                     :: localTablesVersion
    INTEGER                     :: significanceOfReferenceTime
    INTEGER                     :: year
    INTEGER                     :: month
    INTEGER                     :: day
    INTEGER                     :: hour
    INTEGER                     :: minute
    INTEGER                     :: second
  END type t_grib_sec_1
  !
  TYPE t_grib_sec_2
    INTEGER                     :: setLocalDefinition
    !
    INTEGER                     :: localDefinitionNumber
    INTEGER                     :: localHostIdentifier
    INTEGER                     :: localCreationDateYear
    INTEGER                     :: localCreationDateMonth
    INTEGER                     :: localCreationDateDay
    INTEGER                     :: localCreationDateHour
    INTEGER                     :: localCreationDateMinute
    INTEGER                     :: localCreationDateSecond
    INTEGER                     :: localValidityDateYear
    INTEGER                     :: localValidityDateMonth
    INTEGER                     :: localValidityDateDay
    INTEGER                     :: localValidityDateHour
    INTEGER                     :: localValidityDateMinute
    INTEGER                     :: localValidityDateSecond
    INTEGER                     :: localInformationNumber
    INTEGER                     :: localVersionNumber
    INTEGER                     :: localNumberOfExperiment
    INTEGER                     :: localTypeOfEnsembleForecast
    !
    INTEGER                     :: localNumberOfRadarStations
    INTEGER                     :: localIndexOfRadarStation
    INTEGER                     :: localCountryId
    INTEGER                     :: localNationalStationId
    CHARACTER(LEN=4)            :: localStationName
    INTEGER                     :: localStationLatitude
    INTEGER                     :: localStationLongitude
    INTEGER                     :: localScaleFactorOfStationHeightAboveMSL
    INTEGER                     :: localScaledValueOfStationHeightAboveMSL
    INTEGER                     :: localScaleFactorOfAntennaElevationAngle
    INTEGER                     :: localScaledValueOfAntennaElevationAngle
    INTEGER                     :: localIndicatorOfUnitOfTimeRange
    INTEGER                     :: localAccumulationInterval
    INTEGER                     :: localTypeOfRadar
    INTEGER                     :: localOperatingMode
    INTEGER                     :: localQualityControlIndicator
    INTEGER                     :: localClutterFilterIndicator
    INTEGER                     :: localTypeOfDataMaskedByBitmap
    INTEGER                     :: localScaleFactorOfBitmapValue
    INTEGER                     :: localScaledValueOfBitmapValue
    INTEGER                     :: localNumberOfScaledMetadata
    INTEGER, ALLOCATABLE        :: localIndexOfScaledMetadata(:)
    INTEGER, ALLOCATABLE        :: localScaledMetadataTableNumber(:)
    INTEGER, ALLOCATABLE        :: localScaledMetadataId(:)
    INTEGER, ALLOCATABLE        :: localScaleFactorOfMetadata(:)
    INTEGER, ALLOCATABLE        :: localScaledValueOfMetadata(:)
    TYPE(t_scaledMetadataList), ALLOCATABLE :: scaledMetadataList(:)
    ! Unscaled values for prospecitve reader routines:
    REAL(dp)                    :: stationLatitudeInDegrees
    REAL(dp)                    :: stationLongitudeInDegrees
    REAL(dp)                    :: stationHeight
    REAL(dp)                    :: elevationAngle
  END type t_grib_sec_2
  !
  TYPE t_grib_sec_3
    INTEGER                     :: gridDefinitionTemplateNumber
    INTEGER                     :: numberOfDataBinsAlongRadials
    INTEGER                     :: numberOfRadials
    INTEGER                     :: latitudeOfCenterPoint
    INTEGER                     :: longitudeOfCenterPoint
    INTEGER                     :: spacingOfBinsAlongRadials
    INTEGER                     :: offsetFromOriginToInnerBound
    INTEGER                     :: iScansNegatively
    INTEGER                     :: jScansPositively
    INTEGER                     :: jPointsAreConsecutive
    INTEGER, ALLOCATABLE        :: startingAzimuth(:)
    INTEGER, ALLOCATABLE        :: azimuthalWidth(:)
  END type t_grib_sec_3
  !
  TYPE t_grib_sec_4
    INTEGER                     :: productDefinitionTemplateNumber
    INTEGER                     :: parameterCategory
    INTEGER                     :: parameterNumber
    INTEGER                     :: typeOfGeneratingProcess
    INTEGER                     :: backgroundProcess
    INTEGER                     :: generatingProcessIdentifier
    INTEGER                     :: hoursAfterDataCutoff
    INTEGER                     :: minutesAfterDataCutoff
    INTEGER                     :: indicatorOfUnitOfTimeRange
    INTEGER                     :: forecastTime
    INTEGER                     :: typeOfFirstFixedSurface
    INTEGER                     :: scaleFactorOfFirstFixedSurface
    INTEGER                     :: scaledValueOfFirstFixedSurface
    INTEGER                     :: typeOfSecondFixedSurface
    INTEGER                     :: scaleFactorOfSecondFixedSurface
    INTEGER                     :: scaledValueOfSecondFixedSurface
    INTEGER                     :: typeOfEnsembleForecast
    INTEGER                     :: perturbationNumber
    INTEGER                     :: numberOfForecastsInEnsemble
    INTEGER                     :: numberOfRadarSitesUsed
    INTEGER                     :: siteLatitude
    INTEGER                     :: siteLongitude
    INTEGER                     :: siteElevation
    INTEGER                     :: siteId
    INTEGER                     :: operatingMode
    INTEGER                     :: reflectivityCalibrationConstant
    INTEGER                     :: qualityControlIndicator
    INTEGER                     :: clutterFilterIndicator
    INTEGER                     :: constantAntennaElevationAngle
    INTEGER                     :: accumulationInterval
    INTEGER                     :: referenceReflectivityForEchoTop
    INTEGER                     :: rangeBinSpacing
    INTEGER                     :: radialAngularSpacing
  END type t_grib_sec_4
  !
  TYPE t_grib_sec_5
    INTEGER                     :: dataRepresentationTemplateNumber
    INTEGER                     :: bitsPerValue
    INTEGER                     :: typeOfOriginalFieldValues
    CHARACTER(LEN=MAX_LEN_TAB_ENTRY):: packingType                 ! Art der Kompression
  END type t_grib_sec_5
  !
  TYPE t_grib_sec_6
    INTEGER                     :: bitmapPresent
    REAL(dp)                    :: bitmapValue
  END type t_grib_sec_6
  !
  TYPE t_grib_sec_7
    REAL(dp), ALLOCATABLE       :: values(:)
    REAL(dp), ALLOCATABLE       :: radials(:)
    REAL(dp), ALLOCATABLE       :: azimuths(:)
  END type t_grib_sec_7
  !
  TYPE t_grib_loc_sec_tab
    INTEGER,                          ALLOCATABLE :: codeFigure(:)
    CHARACTER(LEN=MAX_LEN_TAB_ENTRY), ALLOCATABLE :: meaning(:)
    INTEGER                                       :: tableNumber
    INTEGER                                       :: nentry
  END type t_grib_loc_sec_tab
  !
  TYPE t_grib
    ! Sektionen:
    TYPE(t_grib_sec_0) :: sec0
    TYPE(t_grib_sec_1) :: sec1
    TYPE(t_grib_sec_2) :: sec2
    TYPE(t_grib_sec_3) :: sec3
    TYPE(t_grib_sec_4) :: sec4
    TYPE(t_grib_sec_5) :: sec5
    TYPE(t_grib_sec_6) :: sec6
    TYPE(t_grib_sec_7) :: sec7
    ! Tabellen:
    TYPE(t_grib_loc_sec_tab) :: tab2
    ! Sample file:
    CHARACTER(LEN=MAX_LEN_FILENAME) :: sample_file
    ! Other general information:
    INTEGER :: product_type
    INTEGER :: bitmap_priority
    INTEGER :: localOperatingMode_in
    CHARACTER(len=LEN_DATETIME) :: creationDate
    CHARACTER(len=LEN_DATETIME) :: dateTime_ini
    CHARACTER(len=LEN_DATETIME) :: dateTime_act
    INTEGER :: year, month, day, hour, minute, second
    INTEGER :: forecastseconds
    INTEGER :: forecastminutes
    CHARACTER(len=12) :: varname
    CHARACTER(len=12) :: packingtype      ! 'grid_simple', 'grid_ccsds', 'png'
  END type t_grib
  !
  TYPE t_grib2_modelspec
    INTEGER :: tablesVersion = 19         ! Main switch for table version, default init with 19
    INTEGER :: &                          ! Table 1.2
      & significanceOfReferenceTime       ! 0: Analysis
                                          ! 1: Start of forecast
                                          ! 2: Verifying time of forecast
                                          ! 4: ...

    INTEGER :: &
      & productDefinitionTemplateNumber   ! 1: individual ensemble forecast
                                          ! 0: standard product

    INTEGER :: &                          ! Table 1.3
      & productionStatusOfProcessedData   ! 0: Oper. products
                                          ! 1: Oper. test products
                                          ! 2: Research products
                                          ! 3: ...

    INTEGER :: &                          ! Table 1.4
      & typeOfProcessedData               ! 0: Analysis products
                                          ! 1: Forecast products
                                          ! 2: Analysis and forecast products
                                          ! 3: ...

    INTEGER :: &                          ! Table 4.3
      & typeOfGeneratingProcess           ! 0: Analysis
                                          ! 1: Initialization
                                          ! 2: Forecast
                                          ! 3: ...
                                          ! 196: invariant data

    INTEGER :: &                          ! Table: backgroundProcess
      & backgroundProcess                 ! 0: main run
                                          ! 1: pre-assimilation
                                          ! 2: assimilation
                                          ! 3: ...

    INTEGER :: &                          ! Table: generatingProcessIdentifier
      & generatingProcessIdentifier       ! 1: icogl
                                          ! 2: icrgl
                                          ! 3: icoeu
                                          ! 4: ...

    INTEGER :: &                          ! Table: local.78.254.def
      & localDefinitionNumber             ! 252: Ensemble system incl. postprocessing
                                          ! 253: Ensemble system
                                          ! 254: Deterministic system

    INTEGER :: &                          ! Table: local.78.254.def
      & localNumberOfExperiment           !


    INTEGER :: &                          ! Output generating center
      & generatingCenter                  ! DWD  : 78
                                          ! ECMWF: 98


    INTEGER :: &                          ! Output generating subcenter
      & generatingSubcenter               ! DWD  : 255
                                          ! ECMWF: 0

    INTEGER :: typeOfEnsembleForecast,        &
      &        localTypeOfEnsembleForecast,   &
      &        numberOfForecastsInEnsemble,   &
      &        perturbationNumber

  END TYPE t_grib2_modelspec

  ! Enumerators for data types in t_attr:
  INTEGER, PARAMETER, PUBLIC :: INT_ATT=1, REAL_ATT=2, CHAR_ATT=3

  TYPE t_attr
    INTEGER                          :: dtype    = -1
    CHARACTER(len=MAX_LEN_ATTR_NAME) :: name     = REPEAT(' ', MAX_LEN_ATTR_NAME)
    INTEGER                          :: int_val  = 0
    REAL(dp)                         :: real_val = 0.0_dp
    CHARACTER(len=MAX_LEN_ATTR)      :: char_val = REPEAT(' ', MAX_LEN_ATTR)
  END TYPE t_attr
  
  TYPE t_cdfin_globalatt
    ! Some fixed global attributes:
    INTEGER                     :: product_type
    CHARACTER(len=LEN_DATETIME) :: Creation_date    = REPEAT(' ', LEN_DATETIME)
    CHARACTER(len=MAX_LEN_ATTR) :: Creator          = REPEAT(' ', MAX_LEN_ATTR)
    CHARACTER(len=MAX_LEN_ATTR) :: Data_source      = REPEAT(' ', MAX_LEN_ATTR)
    CHARACTER(len=MAX_LEN_DESC) :: Data_description = REPEAT(' ', MAX_LEN_DESC)
    
    ! Some more flexible global attributes:
    INTEGER                     :: n_attr    ! total number of flexible global attributes att(:)
!!$ att is declared as allocatable, because a direct att(NATTR_MAX) fails with gfortran 7.4.1:
    TYPE(t_attr), ALLOCATABLE   :: att(:)
  END TYPE t_cdfin_globalatt

  PUBLIC ::  t_scaledMetadataList, t_grib, t_grib2_modelspec, t_grib_loc_sec_tab, t_cdfin_globalatt

END MODULE radar_data_io
