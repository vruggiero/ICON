! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Module containing utility routines for setting GRIB2 keys.

MODULE mo_grib2_util

  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_exception,          ONLY: finish
  USE mo_cdi,                ONLY: streamInqVlist, cdiInqKeyInt, cdiDefKeyInt,          &
    &                              vlistInqVarTsteptype, vlistInqTaxis, taxisInqTunit,  &
    &                              TSTEP_CONSTANT, TSTEP_AVG, TSTEP_ACCUM, TSTEP_MAX,   &
    &                              TSTEP_MIN, TUNIT_SECOND, TUNIT_MINUTE, TUNIT_HOUR,   &
    &                              vlistDefVarIntKey, CDI_KEY_TYPEOFGENERATINGPROCESS,  &
    &                              CDI_KEY_PRODUCTDEFINITIONTEMPLATE, CDI_NOERR,        &
    &                              zaxisInqLbounds, zaxisInqUbounds
  USE mo_gribout_config,     ONLY: t_gribout_config, GRIB_UNDEFVAL, GRIB_NOTINUSEVAL, &
    &                              GRIB_LIB_COMPAT_ECC_2_31_0, GRIB_MAX_NUM_MOD_COMP
  USE mo_var_metadata_types, ONLY: t_var_metadata, CLASS_TILE, CLASS_SYNSAT, &
    &                              CLASS_CHEM, CLASS_TILE_LAND,              &
    &                              CLASS_CHEM_STAT, CLASS_CHEM_OPTP,         &
    &                              CLASS_DISTR, CLASS_DISTR_STAT
  USE mo_action_types,       ONLY: ACTION_RESET
#ifndef __NO_ICON_ATMO__
  USE mo_lnd_nwp_config,     ONLY: tile_list
  USE mo_nwp_sfc_tiles,      ONLY: t_tileinfo_icon, t_tileinfo_grb2
#endif
  USE mo_grib2_tile,         ONLY: t_grib2_template_tile
  ! calendar operations
  USE mtime,                 ONLY: timedelta, datetime,                &
    &                              OPERATOR(-), getTotalSecondsTimeDelta
  USE mo_grib2,              ONLY: t_grib2_key_list

  IMPLICIT NONE

  PRIVATE

  ! Subroutines/Functions
  PUBLIC :: set_GRIB2_additional_keys
  PUBLIC :: set_GRIB2_ensemble_keys
  PUBLIC :: set_GRIB2_synsat_keys
  PUBLIC :: set_GRIB2_local_keys
  PUBLIC :: set_GRIB2_tile_keys
  PUBLIC :: set_GRIB2_chem_keys
  PUBLIC :: set_GRIB2_timedep_keys
  PUBLIC :: set_GRIB2_timedep_local_keys
  PUBLIC :: grib_lib_compatibility

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Set additional GRIB2 keys
  !
  !
  ! Set additional GRIB2 keys. These are added to each single variable, even though
  ! adding it to the vertical or horizontal grid description may be more elegant.
  !
  SUBROUTINE set_GRIB2_additional_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf
    !
    ! Local
    INTEGER :: steptype
    INTEGER :: res
    CHARACTER(len=*), PARAMETER :: routine = 'set_GRIB2_additional_keys'

    ! inquire steptype (needed below)
    steptype = vlistInqVarTsteptype(vlistID, varID)


    ! SECTION 1: Identification Section
    !
    ! Load correct tables
    !
    CALL vlistDefVarIntKey(vlistID, varID, "tablesVersion", grib_conf%tablesVersion)
    CALL vlistDefVarIntKey(vlistID, varID, "localTablesVersion", grib_conf%localTablesVersion)    
    !
    CALL vlistDefVarIntKey(vlistID, varID, "significanceOfReferenceTime",     &
      &                    grib_conf%significanceOfReferenceTime)
    CALL vlistDefVarIntKey(vlistID, varID, "productionStatusOfProcessedData", &
      &                    grib_conf%productionStatusOfProcessedData)
    CALL vlistDefVarIntKey(vlistID, varID, "typeOfProcessedData",             &
      &                    grib_conf%typeOfProcessedData)


    ! SECTION 3: Grid definition Section
    CALL vlistDefVarIntKey(vlistID, varID, "shapeOfTheEarth", 6)



    ! SECTION 4: Product definition
    CALL vlistDefVarIntKey(vlistID, varID, "backgroundProcess",               &
      &                    grib_conf%backgroundProcess)
    ! in case of lon-lat output, "1" has to be added to generatingProcessIdentifier
    CALL vlistDefVarIntKey(vlistID, varID, "generatingProcessIdentifier",     &
      &                    grib_conf%generatingProcessIdentifier)


    IF (ANY((/TSTEP_AVG,TSTEP_ACCUM,TSTEP_MAX,TSTEP_MIN/) == steptype)) THEN
      ! Always set
      !   typeOfTimeIncrement = 2 
      !   "Successive times processed have same start time of forecast, 
      !    forecast time is incremented"
      ! since this is the only type of time processing available in ICON
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfTimeIncrement", 2)
    ENDIF

    IF (grib_conf%lspecialdate_invar) THEN
      ! Use special date for invariant and climatological fields
      !
      IF ( steptype == TSTEP_CONSTANT ) THEN
        ! invariant data
        res = cdiDefKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, 196)
      ELSE
        res = cdiDefKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, grib_conf%typeOfGeneratingProcess)
      ENDIF
    ELSE
      ! no special treatment of invariant and climatological fields
      !
      res = cdiDefKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, grib_conf%typeOfGeneratingProcess)
    ENDIF
    !
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when defining typeOfGeneratingProcess")
    ENDIF

  END SUBROUTINE set_GRIB2_additional_keys



  !------------------------------------------------------------------------------------------------

  !! Set local GRIB2 keys (SECTION 2)
  !!
  !! Set DWD specific local keys (SECTION 2)
  !!
  SUBROUTINE set_GRIB2_local_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

  !----------------------------------------------------------------

    ! SECTION 2: Initialize local use section
    CALL vlistDefVarIntKey(vlistID, varID, "grib2LocalSectionPresent", 1)


    ! SECTION 2: specific settings (local use) for DWD and COSMO Partners
    ! similar modifications are done in subroutine set_GRIB2_timedep_local_keys below
    !
    ! list of centers:
    !   Moscow:     4   (rums Moscow)
    !               5   Moscow (WMC)
    !              76   Moscow (RSMC/RAFC)
    !   DWD:       78   Offenbach (RSMC)
    !   Italy:     80   Rome (RSMC)
    !   Athens:    96
    !   MCH:      215   Zuerich
    !   Poland:   220   Poland  
    !   Romania:  242   Romania
    !             250   cosmo COnsortium for Small scale MOdelling (COSMO)

    !        subcenter = 255
    !US we do not need to check the subcenter here
    !US IF ((grib_conf%generatingCenter == 78) .AND. (grib_conf%generatingSubcenter == 255)) THEN

    ! Add more centers, if required
    SELECT CASE (grib_conf%generatingCenter)
    CASE (78, 80, 215)

      CALL vlistDefVarIntKey(vlistID, varID, "localDefinitionNumber"  ,         &
        &                    grib_conf%localDefinitionNumber)

      CALL vlistDefVarIntKey(vlistID, varID, "localNumberOfExperiment",         &
        &                    grib_conf%localNumberOfExperiment)



      IF (grib_conf%localDefinitionNumber == 254) THEN
        !
        ! -------------------------------------------
        ! Local definition for deterministic forecast
        ! -------------------------------------------

      ELSE IF (grib_conf%localDefinitionNumber == 253) THEN
        !
        ! --------------------------------------
        ! Local definition for ensemble products
        ! --------------------------------------

        IF (grib_conf%localTypeOfEnsembleForecast /= GRIB_UNDEFVAL)                       &
          &   CALL vlistDefVarIntKey(vlistID, varID, "localTypeOfEnsembleForecast" ,      &
          &                          grib_conf%localTypeOfEnsembleForecast)

      ELSE IF (grib_conf%localDefinitionNumber == 230) THEN
        !
        ! --------------------------------------
        ! Local definition for model composition
        ! --------------------------------------

        CALL set_GRIB2_local_model_composition_keys(vlistID, varID, grib_conf)

      END IF ! localDefinitionNumber
    END SELECT

  END SUBROUTINE set_GRIB2_local_keys


  !------------------------------------------------------------------------------------------------

  !! Set GRIB2 ensemble keys (SECTION 4)
  !!
  !! ATTENTION: To be called AFTER set_GRIB2_additional_keys
  !!            due to its dependency on typeOfGeneratingProcess that may 
  !!            be changed for invariant fields in set_GRIB2_additional_keys
  !!
  SUBROUTINE set_GRIB2_ensemble_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    INTEGER :: typeOfGeneratingProcess
    INTEGER :: res
    CHARACTER(len=*), PARAMETER :: routine = 'set_GRIB2_ensemble_keys'
  !----------------------------------------------------------------

    ! get typeOfGeneratingProcess
    ! We do not make use of grib_conf%typeOfGeneratingProcess, since 
    ! typeOfGeneratingProcess is modified for invariant fields.
    res = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, typeOfGeneratingProcess)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when inquiring typeOfGeneratingProcess")
    ENDIF

    IF (typeOfGeneratingProcess == 4) THEN  ! Ensemble forecast

      ! SECTION 4: Product definition Section

      IF (grib_conf%typeOfEnsembleForecast /= GRIB_UNDEFVAL)                            &
        &   CALL vlistDefVarIntKey(vlistID, varID, "typeOfEnsembleForecast" ,           &
        &                          grib_conf%typeOfEnsembleForecast)
      IF (grib_conf%numberOfForecastsInEnsemble /= GRIB_UNDEFVAL)                       &
        &   CALL vlistDefVarIntKey(vlistID, varID, "numberOfForecastsInEnsemble" ,      &
        &                          grib_conf%numberOfForecastsInEnsemble)
      IF (grib_conf%perturbationNumber /= GRIB_UNDEFVAL)                                &
        &   CALL vlistDefVarIntKey(vlistID, varID, "perturbationNumber" ,               &
        &                          grib_conf%perturbationNumber)
    END IF ! typeOfGeneratingProcess

  END SUBROUTINE set_GRIB2_ensemble_keys


  !! Set GRIB2 keys which are specific to synthetic satellite products.
  !!
  SUBROUTINE set_GRIB2_synsat_keys (vlistID, varID, info)
    
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),  INTENT(IN) :: info

    ! Local
    INTEGER :: typeOfGeneratingProcess
    INTEGER :: res
    CHARACTER(len=*), PARAMETER :: routine = 'set_GRIB2_synsat_keys'
    ! ----------------------------------------------------------------

    ! Skip inapplicable fields
    IF ( info%var_class /= CLASS_SYNSAT ) RETURN

    ! get typeOfGeneratingProcess
    res = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, typeOfGeneratingProcess)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when inquiring typeOfGeneratingProcess")
    ENDIF

    ! change product definition template
    IF (typeOfGeneratingProcess == 4) THEN
      ! Ensemble forecast
      res = cdiDefKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, 33)
    ELSE
      ! Deterministic forecast
      res = cdiDefKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, 32)
    END IF
    !
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when defining productDefinitionTemplate")
    ENDIF

  END SUBROUTINE set_GRIB2_synsat_keys


  !! Set keys specific to atmospheric chemical species
  !!
  !! Set GRIB2 keys which are specific to atmospheric chemical species.
  !! Here, only the PDT will be changed. Additional Template-specific 
  !! keys will be set at the end of 
  !! mo_name_list_output_init:add_variables_to_vlist
  !!
  SUBROUTINE set_GRIB2_chem_keys (vlistID, varID, info)
    
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),  INTENT(IN) :: info

    ! Local
    INTEGER :: typeOfGeneratingProcess
    INTEGER :: productDefinitionTemplate     ! template number
    INTEGER :: res
    CHARACTER(len=*), PARAMETER :: routine = 'set_GRIB2_chem_keys'
    ! ----------------------------------------------------------------
    
    SELECT CASE (info%var_class)
    CASE (CLASS_CHEM)
      productDefinitionTemplate = 40
    CASE (CLASS_CHEM_STAT)
      productDefinitionTemplate = 42
    CASE (CLASS_CHEM_OPTP)
      productDefinitionTemplate = 48
    CASE (CLASS_DISTR)
      productDefinitionTemplate = 57
    CASE (CLASS_DISTR_STAT)
      productDefinitionTemplate = 67
    CASE DEFAULT
      ! skip inapplicable fields
      RETURN
    END SELECT

    ! get typeOfGeneratingProcess
    res = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, typeOfGeneratingProcess)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when inquiring typeOfGeneratingProcess")
    ENDIF

    ! change product definition template in case of ensemble run
    IF (typeOfGeneratingProcess == 4)  &
      &  productDefinitionTemplate = productDefinitionTemplate + 1

    ! set product definition template
    res = cdiDefKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, productDefinitionTemplate)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when defining productDefinitionTemplate")
    ENDIF

  END SUBROUTINE set_GRIB2_chem_keys


  !! Set GRIB2 keys which are specific to tile-variables.
  !!
  SUBROUTINE set_GRIB2_tile_keys (vlistID, varID, info, i_lctype, grib2_template_tile)

    INTEGER,                     INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),       INTENT(IN) :: info
    INTEGER,                     INTENT(IN) :: i_lctype  !< Tile classification
    TYPE(t_grib2_template_tile), INTENT(IN) :: grib2_template_tile ! set of allowed tile templates

#ifndef __NO_ICON_ATMO__
    ! local
    INTEGER                   :: res
    INTEGER                   :: typeOfGeneratingProcess
    INTEGER                   :: productDefinitionTemplate        ! Tile template number 
    INTEGER                   :: natt
    TYPE(t_tileinfo_grb2)     :: tileinfo_grb2
    CHARACTER(len=*), PARAMETER :: routine = 'set_GRIB2_tile_keys'
  !----------------------------------------------------------------

    ! Skip inapplicable fields
    IF ( ALL((/CLASS_TILE,CLASS_TILE_LAND/) /= info%var_class) ) RETURN

    res = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, typeOfGeneratingProcess)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when inquiring typeOfGeneratingProcess")
    ENDIF

    ! change product definition template
    !
    IF (typeOfGeneratingProcess == 4) THEN
      ! ensemble
      productDefinitionTemplate = grib2_template_tile%tpl_inst_ens
    ELSE
      ! deterministic
      productDefinitionTemplate = grib2_template_tile%tpl_inst
    ENDIF

! use the following IF condition, once the WMO tile templates for statistical processing 
! become available (validation by WMO pending).
!
!!$    IF (typeOfGeneratingProcess == 4) THEN
!!$      ! ensemble
!!$      IF (ANY((/TSTEP_MAX, TSTEP_MIN, TSTEP_AVG, TSTEP_ACCUM/) == info%isteptype)) THEN
!!$        productDefinitionTemplate = grib2_template_tile%tpl_acc_ens
!!$      ELSE
!!$        productDefinitionTemplate = grib2_template_tile%tpl_inst_ens
!!$      ENDIF
!!$    ELSE
!!$      ! deterministic
!!$      IF (ANY((/TSTEP_MAX, TSTEP_MIN, TSTEP_AVG, TSTEP_ACCUM/) == info%isteptype)) THEN
!!$        productDefinitionTemplate = grib2_template_tile%tpl_acc
!!$      ELSE
!!$        productDefinitionTemplate = grib2_template_tile%tpl_inst
!!$      ENDIF
!!$    ENDIF

    res = cdiDefKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, productDefinitionTemplate)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when defining productDefinitionTemplate")
    ENDIF

    ! Set tile classification
    CALL vlistDefVarIntKey(vlistID, varID, TRIM(grib2_template_tile%keys%tileClassification), i_lctype)

    ! Set total number of tile/attribute pairs
    CALL vlistDefVarIntKey(vlistID, varID, TRIM(grib2_template_tile%keys%totalNumberOfTileAttributePairs), &
      &                    info%maxcontained)

    ! Set number of used tiles
    CALL vlistDefVarIntKey(vlistID, varID, TRIM(grib2_template_tile%keys%numberOfUsedSpatialTiles), &
      &                    tile_list%getNumberOfTiles(varClass=info%var_class))

    ! get the following attributes:
    ! - tileIndex
    ! - numberOfTileAttributes
    ! - tileAttribute
    !
    ! Select GRIB2 tileinfo object corresponding to internal tile index info%ncontained
    tileinfo_grb2 = tile_list%getTileinfo_grb2( t_tileinfo_icon(info%ncontained) )
    ! get number of GRIB2 tile attributes corresponding to internal tile index info%ncontained
    natt          = tile_list%getNumberOfTileAttributes( t_tileinfo_icon(info%ncontained) )

    ! Set tile index
    CALL vlistDefVarIntKey(vlistID, varID, TRIM(grib2_template_tile%keys%tileIndex), tileinfo_grb2%idx)

    ! Set total number of tile attributes for given tile
    CALL vlistDefVarIntKey(vlistID, varID, TRIM(grib2_template_tile%keys%numberOfUsedTileAttributes), natt)

    ! Set tile attribute
    CALL vlistDefVarIntKey(vlistID, varID, TRIM(grib2_template_tile%keys%attributeOfTile), &
      &                    tileinfo_grb2%att)
#endif

  END SUBROUTINE set_GRIB2_tile_keys




  !------------------------------------------------------------------------------------------------
  !> Set additional, time-dependent GRIB2 keys
  !!  
  !!  This subroutine sets GRIB2 keys that may change during the
  !!  simulation. Currently this is the case for the start and end time 
  !!  of the statistical process interval.
  !!
  !!  Description of how the GRIB2 keys 'forecastTime' and 'lengthOfTimeRange' are computed.
  !!
  !!  |===================*======================|================================> time axis
  !!  ^start_time         |                      ^current_time (output)
  !!                      |                      |
  !!                      ^EventPrevTriggerTime  |
  !!                      |                      |
  !!                      |                      |
  !!  <-------------------><--------------------->
  !!     forecastTime         lengthOfTimeRange
  !!
  !!  CAVEATs
  !!  - we implicitly assume that actions are ordered according to increasing forecast time.
  !!
  SUBROUTINE set_GRIB2_timedep_keys(streamID, varID, info, start_date, cur_date)
    INTEGER                             ,INTENT(IN) :: streamID, varID
    TYPE (t_var_metadata)               ,INTENT(IN) :: info
    !> model start (initialization) time
    TYPE(datetime), INTENT(IN) :: start_date
    !> model current time (rounded)
    TYPE(datetime), INTENT(IN) :: cur_date
    ! local variables
    TYPE(timedelta)               :: mtime_lengthOfTimeRange       ! length of time range (mtime format)
    INTEGER                       :: ilengthOfTimeRange_secs       ! length of time range in seconds
    INTEGER                       :: ilengthOfTimeRange            ! length of time range in appropriate time units
    INTEGER                       :: taxis_tunit                   ! time axis units
    INTEGER                       :: vlistID, taxisID
    TYPE(timedelta)               :: forecast_delta                ! forecast time (mtime format)
    INTEGER                       :: forecast_secs                 ! forecast time in seconds
    INTEGER                       :: forecast_time                 ! forecast time in appropriate time units
    TYPE(datetime)                :: statProc_startDateTime        ! statistical process starting DateTime
    INTEGER                       :: var_actionId                  ! action from which metainfo is used
    CHARACTER(len=*), PARAMETER   :: routine = 'set_GRIB2_timedep_keys'

    ! special fields for which time-dependent metainfos should be set even though they are not of
    ! steptype TSTEP_MAX or TSTEP_MIN. These fields are special in the sense that averaging is not
    ! performed over the entire model run but over only some intervals.
    !CHARACTER(LEN=17) :: ana_avg_vars(16) = (/"rain_gsp         ",                                          &
    !                                        & "snow_gsp         ", "graupel_gsp      ", "ice_gsp          ",&
    !                                        & "hail_gsp         ", "prec_gsp         ", "rain_con         ",&
    !                                        & "snow_con         ", "prec_con         ", "tot_prec         ",&
    !                                        & "tot_prec_d       ", "prec_gsp_d       ", "prec_con_d       ",&
    !                                        & "prec_con_rate_avg", "prec_gsp_rate_avg", "tot_prec_rate_avg"/)


    !---------------------------------------------------------
    ! Set time-dependent metainfo
    !---------------------------------------------------------
    !
    ! Skip inapplicable fields
    ! Currently all TSTEP_AVG and TSTEP_ACC fields are skipped, except for special ones 
    ! listed in ana_avg_vars

    ! UB: this method does not work for the internally interpolated vars of name ana_avg_vars, because
    !     the internal info%name gets a prefix "latlon_internal_".

!    IF ((ALL((/TSTEP_MAX, TSTEP_MIN/) /= info%isteptype)) .AND. &
!      & (one_of(TRIM(info%name),ana_avg_vars) == -1) ) RETURN

    ! DR
    ! less cumbersome solution, which does no longer skip variables of type 
    ! TSTEP_AVG, TSTEP_ACCUM. By this, the manual list ana_avg_vars becomes superfluous.
    !
    ! This version differs from the current version (see above) by the fact that 
    ! forecastTime and lengthOfTimeRange are set for ALL statistically processed variables 
    ! and not only for those of type TSTEP_MAX, TSTEP_MIN and members of the list ana_avg_vars.
    !
    ! First tests showed that metadata for variables of type TSTEP_AVG, TSTEP_ACCUM are 
    ! still correct. Putting it the other way around: It is not clear 
    ! to me why with the currently active version the keys forecastTime and lengthOfTimeRange are 
    ! set correctly for variables of type TSTEP_AVG, TSTEP_ACCUM.
    IF ((ALL((/TSTEP_MAX, TSTEP_MIN, TSTEP_AVG, TSTEP_ACCUM/) /= info%isteptype))) RETURN

    ! get vlistID. Note that the stream-internal vlistID must be used. 
    ! It is obtained via streamInqVlist(streamID)
    vlistID     = streamInqVlist(streamID)
    taxisID     = vlistInqTaxis(vlistID)
    taxis_tunit = taxisInqTunit(taxisID)

    IF (info%action_list%n_actions > 0) THEN

      ! more than one RESET action may be defined for a single variable.
      ! get ID of currently active RESET action
      var_actionId = info%action_list%getActiveAction(ACTION_RESET, cur_date)

      IF (var_actionId == -1) THEN
        write(0,*) 'set_timedependent_GRIB2_keys: no active action of type ACTION_RESET found. '//&
          &             'lengthOfTimeRange may not be set correctly'
        CALL finish (routine, 'Illegal actionId')
      ENDIF

      ! get latest (intended) triggering time, which is equivalent to 
      ! the statistical process starting time
      statProc_startDateTime = info%action_list%action(var_actionId)%EventPrevTriggerDate

    ELSE
      ! If there is no RESET action available, it is assumed that the 
      ! statistical process start time is equal to the model start time
      statProc_startDateTime = start_date
    ENDIF


    ! get time interval, over which statistical process has been performed
    ! It is the time difference between the current time (rounded) cur_date and 
    ! the last time the nullify-action took place (rounded) statProc_startDateTime.
    ! mtime_lengthOfTimeRange = current time - statistical process start time
    mtime_lengthOfTimeRange = cur_date - statProc_startDateTime

    ! time interval over which statistical process has been performed (in secs)    
    ilengthOfTimeRange_secs = INT(getTotalSecondsTimeDelta(mtime_lengthOfTimeRange, &
      &                                                statProc_startDateTime)      &
      &                           )


    ! get forecast_time: forecast_time = statProc_startDateTime - model_startDateTime
    ! Note that for statistical quantities, the forecast time is the time elapsed between the 
    ! model start time and the start time of the statistical process
    forecast_delta = statProc_startDateTime - start_date

    ! forecast time in seconds
    forecast_secs = INT(getTotalSecondsTimeDelta(forecast_delta, start_date))


    SELECT CASE (taxis_tunit)
    CASE (TUNIT_SECOND)
       ilengthOfTimeRange          = ilengthOfTimeRange_secs
       forecast_time               = forecast_secs
    CASE (TUNIT_MINUTE)
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/60)
       forecast_time               = INT(forecast_secs/60)
    CASE (TUNIT_HOUR)
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/3600)
       forecast_time               = INT(forecast_secs/3600)
    CASE DEFAULT
    END SELECT

    ! set forecast time: statProc_startDateTime - model_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "forecastTime", forecast_time) 

    !
    ! set length of time range: current time - statProc_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "lengthOfTimeRange", ilengthOfTimeRange)
    ! Note that if one of the statistics templates 4.8 or 4.11 is selected, the time unit 
    ! (GRIB2 key "indicatorOfUnitForTimeRange") is set automatically by CDI.
    ! It is always set identical to "indicatorOfUnitOFTimeRange"



  END SUBROUTINE set_GRIB2_timedep_keys


  !------------------------------------------------------------------------------------------------
  !! Set time-dependent local GRIB2 keys (SECTION 2)
  !!
  !! Set DWD specific time-dependent local keys (SECTION 2)
  !!
  SUBROUTINE set_GRIB2_timedep_local_keys(streamID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: streamID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    CHARACTER(len=MAX_CHAR_LENGTH)     :: ydate, ytime
    INTEGER :: cent, year, month, day    ! date
    INTEGER :: hour, minute, second      ! time
    INTEGER :: vlistID
  !----------------------------------------------------------------


    ! get vlistID. Note that the stream-internal vlistID must be used. 
    ! It is obtained via streamInqVlist(streamID)
    vlistID = streamInqVlist(streamID)

    ! SECTION 2: specific settings (local use) for DWD and COSMO Partners
    !  for list of centers see subroutine set_GRIB2_local_keys above

    ! Only set local time information for localDefinitionNumber 253, 254, 230
    ! COSMO local section 250 does not include time information

    SELECT CASE (grib_conf%localDefinitionNumber)
    CASE (253, 254, 230)
      SELECT CASE (grib_conf%generatingCenter)
      CASE (78, 80, 215)

        IF (grib_conf%ldate_grib_act) THEN
          ! get date and time
          ! ydate : ccyymmdd, ytime : hhmmss.sss
          CALL date_and_time(ydate,ytime)
          READ(ydate,'(4i2)') cent, year, month, day
          READ(ytime,'(3i2)') hour, minute, second
        ELSE ! set date to "01010101" (for better comparability of GRIB files)
          cent  = 1
          year  = 1
          month = 1
          hour  = 1
          minute= 1
          second= 1
        ENDIF


        CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateYear"  , 100*cent+year)
        CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMonth" , month)
        CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateDay"   , day)
        CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateHour"  , hour)
        CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMinute", minute)
        CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateSecond", second)
        ! CALL vlistDefVarIntKey(vlistID, varID, "localValidityDateYear"  , 2013)

      END SELECT    ! generatingCenter
    END SELECT      ! localDefinitionNumber

  END SUBROUTINE set_GRIB2_timedep_local_keys

  !------------------------------------------------------------------------------------------------
  !! Set local GRIB2 keys for Local-Use-Section (Section 2) template 230 "Model composition":
  !! This template was implemented in the local DWD definitions for ecCodes in the first half-year of 2024.
  !! In production at DWD, ICON-NWP may be coupled together with other model components,
  !! such as the ART model for aerosol and tracer transport and physics, and the ocean model.
  !! Template 230 provides keys to specify such model components in a modular way.
  !! (As the model composition is a highly centre-specific piece of information,
  !! it is incorporated into Section 2 rather than a new official WMO template for Section 4, e.g.)
  !!
  SUBROUTINE set_GRIB2_local_model_composition_keys(vlistID, varID, grib_conf)

    ! Arguments
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local variables
    INTEGER :: localOriginOfProduct,              &
      &        localNumberOfModelComponentsInUse, &
      &        localModelComponent(8)

    CHARACTER(LEN=*), PARAMETER :: routine = 'set_GRIB2_local_model_composition_keys'

    !----------------------------------------------------------------

    ! The following settings are only available for:
    ! * Local-Use-Section (Section 2) template 230
    ! * generatingCenter = 78 / 80 / 215
    IF (grib_conf%localDefinitionNumber /= 230) THEN
      RETURN
    ELSEIF (.NOT. ANY([78, 80, 215] == grib_conf%generatingCenter)) THEN
      RETURN
    ENDIF
    
    ! The following keys of the model composition template
    ! will be set here (further keys are set in set_GRIB2_timedep_local_keys):
    !
    ! * localProductionSystem......................(eccodes local definition table: 2.230)
    ! * localOriginOfProduct.......................(eccodes local definition table: 2.231)
    ! * localNumberOfModelComponentsInUse
    ! * localDrivingModelComponent.................(eccodes local definition table: 2.231.<localOriginOfProduct> = 2.231.1)
    ! * local2nd/3rd/.../8thModelComponent.........(eccodes local definition table: 2.231.<localOriginOfProduct> = 2.231.1)

    ! localOriginOfProduct:
    ! 1 -> "Product is generated by a process that is based on a model or model composition described by the following 17 octets"
    localOriginOfProduct = 1
    ! localNumberOfModelComponentsInUse:
    localNumberOfModelComponentsInUse = COUNT(grib_conf%model_components(:) > GRIB_NOTINUSEVAL)
    ! local(Driving/2nd/3rd/.../8th)ModelComponent:
    ! The template allows to specify 8 model components at most.
    ! Default-value initialization:
    ! 0 -> "Unused = Not in use/Does not apply"
    localModelComponent(:) = GRIB_NOTINUSEVAL
    ! Currently GRIB_MAX_NUM_MOD_COMP < 8
    localModelComponent(1:GRIB_MAX_NUM_MOD_COMP) = grib_conf%model_components(1:GRIB_MAX_NUM_MOD_COMP)

    ! Hand over metadata to CDI
    CALL vlistDefVarIntKey(vlistID, varID, "localProductionSystem", grib_conf%localProductionSystem)
    CALL vlistDefVarIntKey(vlistID, varID, "localOriginOfProduct", localOriginOfProduct)
    CALL vlistDefVarIntKey(vlistID, varID, "localNumberOfModelComponentsInUse", localNumberOfModelComponentsInUse)
    CALL vlistDefVarIntKey(vlistID, varID, "localDrivingModelComponent", localModelComponent(1))
    CALL vlistDefVarIntKey(vlistID, varID, "local2ndModelComponent", localModelComponent(2))
    CALL vlistDefVarIntKey(vlistID, varID, "local3rdModelComponent", localModelComponent(3))
    CALL vlistDefVarIntKey(vlistID, varID, "local4thModelComponent", localModelComponent(4))
    CALL vlistDefVarIntKey(vlistID, varID, "local5thModelComponent", localModelComponent(5))
    CALL vlistDefVarIntKey(vlistID, varID, "local6thModelComponent", localModelComponent(6))
    CALL vlistDefVarIntKey(vlistID, varID, "local7thModelComponent", localModelComponent(7))
    CALL vlistDefVarIntKey(vlistID, varID, "local8thModelComponent", localModelComponent(8))

    IF (grib_conf%localTypeOfEnsembleForecast /= GRIB_UNDEFVAL) &
      & CALL vlistDefVarIntKey(vlistID, varID, "localTypeOfEnsembleForecast", grib_conf%localTypeOfEnsembleForecast)

  END SUBROUTINE set_GRIB2_local_model_composition_keys

  !-------------------------------------------------------------------------
  !> Adjustments for compensating non-backward-compatible
  !! metadata changes after version updates of the GRIB library
  !!
  SUBROUTINE grib_lib_compatibility(grib_lib_compat, vlistID, varID, zaxisID, additional_grib_keys)

    ! Arguments
    INTEGER,                INTENT(IN) :: grib_lib_compat
    INTEGER,                INTENT(IN) :: vlistID
    INTEGER,                INTENT(IN) :: varID
    INTEGER,                INTENT(IN) :: zaxisID
    TYPE(t_grib2_key_list), INTENT(IN) :: additional_grib_keys

    ! Local variables
    INTEGER :: typeOfSecondFixedSurface
    INTEGER :: scaleFactorOfSecondFixedSurface
    INTEGER :: scaledValueOfSecondFixedSurface
    INTEGER :: jkey
    INTEGER :: key_len
    INTEGER :: zaxisLboundsSize, zaxisUboundsSize
    INTEGER :: counter
    LOGICAL :: is_layer

    !---------------------------------

    ! Initialization of local auxiliary variables
    ! (we use GRIB_UNDEFVAL = -1 for initialization with negative integer)
    zaxisLboundsSize                = GRIB_UNDEFVAL
    zaxisUboundsSize                = GRIB_UNDEFVAL
    is_layer                        = .FALSE.
    typeOfSecondFixedSurface        = GRIB_UNDEFVAL
    scaleFactorOfSecondFixedSurface = GRIB_UNDEFVAL
    scaledValueOfSecondFixedSurface = GRIB_UNDEFVAL

    IF (grib_lib_compat == GRIB_LIB_COMPAT_ECC_2_31_0) THEN

      ! CDI uses the ecCodes sample file "GRIB2.tmpl" as a starting file for ecCodes.
      ! The SecondFixedSurface keys are assigned the following values in GRIB2.tmpl:
      !
      ! - typeOfSecondFixedSurface        = 255 (Missing)
      ! - scaleFactorOfSecondFixedSurface = MISSING (= 255)
      ! - scaledValueOfSecondFixedSurface = MISSING (= 2,147,483,647)
      !
      ! Now, if typeOfSecondFixedSurface is set to a value >= 10, 102 say,
      ! but scaleFactorOfSecondFixedSurface and scaledValueOfSecondFixedSurface
      ! are not explicitly set, the result is as follows:
      !
      ! For ecCodes version < 2.32.0:
      ! -----------------------------
      ! - typeOfSecondFixedSurface        = 102 (Specific altitude above mean sea level)
      ! - scaleFactorOfSecondFixedSurface = 0
      ! - scaledValueOfSecondFixedSurface = 0
      !
      ! For ecCodes version >= 2.32.0:
      ! ------------------------------
      ! - typeOfSecondFixedSurface        = 102 (Specific altitude above mean sea level)
      ! - scaleFactorOfSecondFixedSurface = MISSING
      ! - scaledValueOfSecondFixedSurface = MISSING
      !
      ! With grib_lib_compat = 1, we try to overwrite the behavior of
      ! ecCodes versions >= 2.32.0 with the behavior of versions < 2.32.0, in this respect.
      
      ! Check if the SecondFixedSurface keys are among the additional integer GRIB keys
      counter = 0
      KEY_LOOP: DO jkey = 1, additional_grib_keys%nint_keys
        key_len = LEN_TRIM(additional_grib_keys%int_key(jkey)%key)
        ! Avoid unnecessary string comparisons:
        ! - LEN(typeOfSecondFixedSurface) = 24
        ! - LEN(scaleFactorOfSecondFixedSurface) = LEN(scaledValueOfSecondFixedSurface) = 31
        IF (key_len == 24) THEN
          IF (additional_grib_keys%int_key(jkey)%key(1:key_len) == "typeOfSecondFixedSurface") THEN
            typeOfSecondFixedSurface = additional_grib_keys%int_key(jkey)%val
            counter = counter + 1
          END IF
        ELSEIF (key_len == 31) THEN
          IF (additional_grib_keys%int_key(jkey)%key(1:key_len) == "scaleFactorOfSecondFixedSurface") THEN
            scaleFactorOfSecondFixedSurface = additional_grib_keys%int_key(jkey)%val
            counter = counter + 1
          ELSEIF (additional_grib_keys%int_key(jkey)%key(1:key_len) == "scaledValueOfSecondFixedSurface") THEN
            scaledValueOfSecondFixedSurface = additional_grib_keys%int_key(jkey)%val
            counter = counter + 1
          ENDIF
        ENDIF
        ! If there are 3 hits, we can already leave the search loop
        IF (counter == 3) EXIT KEY_LOOP
      END DO KEY_LOOP
      
      IF (typeOfSecondFixedSurface > 9) THEN
        ! Figure out if zaxis is level- or layer-based:
        ! In case a zaxis is layer-based, the lower and upper bounds of a layer - 
        ! zaxisLbounds and zaxisUbounds - should have been specified during the creation of the zaxis.
        ! In this case, the CDI functions zaxisInqLbounds and zaxisInqUbounds should return a size > 0.
        zaxisLboundsSize = zaxisInqLbounds(zaxisID)
        zaxisUboundsSize = zaxisInqUbounds(zaxisID)
        is_layer         = ((zaxisLboundsSize > 0) .AND. (zaxisUboundsSize > 0))
        ! If a zaxis is layer-based, CDI should set the values of
        ! the GRIB keys 'scaleFactorOfSecondFixedSurface' and 'scaledValueOfSecondFixedSurface'.
        ! Therefore, we should not overwrite these values here,
        ! even though the key 'typeOfSecondFixedSurface' was explicitly set in the 'add_var'
        ! of the respective variable (for whatever reason).
        IF (.NOT. is_layer) THEN
          ! We will set the value of the following two GRIB keys to zero,
          ! only if they are not explicitly specified in 'add_var'
          IF (scaleFactorOfSecondFixedSurface < 0) CALL vlistDefVarIntKey(vlistID, varID, "scaleFactorOfSecondFixedSurface", 0)
          IF (scaledValueOfSecondFixedSurface < 0) CALL vlistDefVarIntKey(vlistID, varID, "scaledValueOfSecondFixedSurface", 0)
        ENDIF
      ENDIF ! IF (typeOfSecondFixedSurface > 9)

    ENDIF ! IF (grib_lib_compat == GRIB_LIB_COMPAT_ECC_2_31_0)
      
  END SUBROUTINE grib_lib_compatibility

END MODULE mo_grib2_util
