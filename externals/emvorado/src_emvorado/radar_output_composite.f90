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


MODULE radar_output_composite

  !------------------------------------------------------------------------------
  !
  ! Description: Methods of the radar forward operator EMVORADO for writing
  !              reflectivity composites on a rotated lat/lon grid to grib2 files.
  !
  ! Method:
  !   See subroutines below
  !
  !------------------------------------------------------------------------------
  !
  ! Declarations:
  !
  ! Modules used:
  !

  USE radar_kind, ONLY :  dp
  
  USE radar_data, ONLY :  miss_value

  USE radar_interface, ONLY : &
       abort_run,          &
       get_model_time_ddhhmmss, &
       get_datetime_act, &
       get_datetime_ini, &
       grib2_add_modelspec_info

  USE radar_utilities, ONLY : new_datetime, diff_seconds

  USE radar_data, ONLY : &
       cmaxlen, &
       ydate_ini_mod, &
       idom, &
       my_radar_id, & ! rank of this PE in the radar communicator (cart+radario)
       composite_meta_type

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, &
       comp_grib2_packingtype

  USE radar_output_utils,  ONLY : get_next_key_from_pos, replace_substr_with_value, &
                                  check_codes_err
  
  USE radar_data_io, ONLY : GRIB_TRUE, GRIB_FALSE, t_grib2_modelspec, &
                            TGP_ENSEMBLE

#ifdef GRIBAPI
  USE eccodes, ONLY : codes_clone, codes_set, codes_write, codes_release, codes_close_file, &
                      codes_grib_new_from_samples, codes_open_file, GRIB_SUCCESS, &
                      codes_get_error_string, codes_set_missing
#endif

  !================================================================================
  !================================================================================

  IMPLICIT NONE

  !================================================================================
  !================================================================================

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !==============================================================================
  ! Public and Private Subroutines

  PRIVATE

#ifdef GRIBAPI
  PUBLIC ::  write_composite_grib
#endif

  !==============================================================================
  ! Module variables

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS


#ifdef GRIBAPI
  SUBROUTINE write_composite_grib(outdir, file_pattern, cmp_meta, comptyp, ilevel, comp2d_tot, error, errmsg)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: writes the radar composite to a grib file. The composite is
    !              taken "as is" for the timestep at which this procedure is called.
    !              This means that it contains all radar stations which have
    !              an obs_time within +/- 0.5*dt_model of the actual model time.
    !
    ! Note: the actual time in the filename and in the attributes is the actual
    !       model time. It is NOT necessarily the same as the obs_time of a radar
    !       station, just closer than +/- 0.5*dt_model to such a time!
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! .. Input/Output variables:
    TYPE(composite_meta_type), INTENT(in) :: cmp_meta
    CHARACTER(*)             , INTENT(in) :: file_pattern, comptyp, outdir
    CHARACTER(*)             , INTENT(out):: errmsg
    REAL(kind=dp)            , INTENT(in) :: comp2d_tot(:,:)
    INTEGER                  , INTENT(in) :: ilevel  ! integer for level coding in the grib:
                                                     ! 0     = composite for bubble generator
                                                     ! 1...9 = 9 different elevations depending on namelist settings
                                                     ! 98    = composite of the DWD precipitation scans
                                                     ! 99    = maximum composite using all elevations
    INTEGER                  , INTENT(out) :: error

    ! .. Local variables:
    TYPE(t_grib2_modelspec) :: grib_locinfo
    INTEGER                 :: igribsample, igribclone, igribout, ierr, ilocerr
    CHARACTER(len=cmaxlen)  :: samplefile, outfile

    INTEGER            :: ie_loc, je_loc, isecdiff, ierr_diff(2,6)
    INTEGER            :: forecastminutes, idate, itime, isecond
    INTEGER            :: i, j, k, ij
    INTEGER            :: localCreationDateYear, localCreationDateMonth, localCreationDateDay, &
                          localCreationDateHour, localCreationDateMinute, localCreationDateSecond, &
                          datevalues(8)
    CHARACTER(len=32)  :: yzroutine
    CHARACTER(len=20)  :: cvar
    CHARACTER(len=LEN(comptyp))  :: zcomptyp
    CHARACTER(len=14)  :: model_starttime, model_validtime, utcdate
    CHARACTER(len=8 )  :: currdate
    CHARACTER(len=10)  :: currtime
    CHARACTER(len=1)   :: cmode
    REAL(kind=KIND(1.0)), ALLOCATABLE :: comp_tot_grb(:)
    LOGICAL :: fileexist

    yzroutine(:) = ' '
    yzroutine    = 'write_composite_grib'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine)//' on proc ', my_radar_id, TRIM(comptyp), ilevel

    error = 0

    ie_loc = SIZE(comp2d_tot, 1)
    je_loc = SIZE(comp2d_tot, 2)

    IF (ie_loc /= cmp_meta%ni .OR. je_loc /= cmp_meta%nj) THEN
      ! This is a programming error, so abort:
      errmsg(:) = ' '
      WRITE (errmsg, '(4(a,i4))') 'ERROR '//TRIM(yzroutine)//'(): '// &
           'dimension of input field different from composite definition: '// &
           'ni ', ie_loc, ' /= ', cmp_meta%ni, ';  nj ', je_loc, ' /= ', cmp_meta%nj
      error = 1
      RETURN
    END IF

    samplefile(:) = ' '
    samplefile = 'DWD_rotated_ll_7km_G_grib2'

    SELECT CASE (TRIM(comptyp))
    CASE ('obs','sim')
      cvar(:) = ' '
      cvar    = 'dbzcmp'
      zcomptyp(:) = ' '
      zcomptyp    = comptyp
    CASE ('obs_bub','sim_bub')
      cvar(:) = ' '
      cvar    = 'dbzcmpbub'
      ! Remove the suffix '_obs' from the comptyp:
      zcomptyp(:) = ' '
      zcomptyp = comptyp(1:3)
    CASE default
      errmsg(:) = ' '
      WRITE (errmsg,*) 'ERROR '//TRIM(yzroutine)//': Unknown composite type '//TRIM(comptyp)
      error = 2
      RETURN
    END SELECT

    model_starttime = ydate_ini_mod
    model_validtime = get_datetime_act()
    CALL diff_seconds ( model_starttime, model_validtime, isecdiff, ierr_diff )
    IF (isecdiff < 0) THEN
      ! HANDWAVING FIX for grib2 not beeing able to have negative forecast times:
      !  Set the model start time to the actual forecast time:
      model_starttime = model_validtime
      ! Set the forecast time to 0 minutes:
      forecastminutes = 0
    ELSE
      forecastminutes = ABS( NINT( REAL(isecdiff,dp) / 60.0_dp ) )
    END IF

    CALL DATE_AND_TIME(date=currdate, time=currtime, values=datevalues(:))  ! this in is in local time, but we need UTC
    ! Convert to UTC:
    ! datevalues(4) contains the difference to UTC in minutes. Use get_utc_date() to compute UTC time:
    utcdate = new_datetime(currdate(1:8)//currtime(1:6), -datevalues(4)*60.0_dp)
    READ(utcdate,'(i4,5i2)') localCreationDateYear,   localCreationDateMonth, &
                             localCreationDateDay,    localCreationDateHour, &
                             localCreationDateMinute, localCreationDateSecond

    outfile(:) = ' '
    CALL composite_create_filename (file_pattern, REAL(isecdiff,dp), TRIM(cvar), TRIM(zcomptyp), outfile)
    outfile = TRIM(outdir)//TRIM(outfile)

    CALL codes_grib_new_from_samples(igribsample, TRIM(samplefile), ilocerr)
    IF (ilocerr /= GRIB_SUCCESS) THEN
      errmsg(:) = ' '
      CALL codes_get_error_string(ilocerr, errmsg)
      errmsg = 'ERROR '//TRIM(yzroutine)//' reading grib sample file '//TRIM(samplefile)//': '//TRIM(errmsg)
      error = 3
      RETURN
    END IF

    INQUIRE (file=TRIM(outfile), exist=fileexist)
    IF (fileexist) THEN
      cmode = 'a'   ! append
      IF (ldebug_radsim) WRITE (*,'(A)') 'Appending to '//TRIM(outfile)
    ELSE
      cmode = 'w'   ! write
      IF (ldebug_radsim) WRITE (*,'(A)') 'Creating '//TRIM(outfile)
    END IF
    CALL codes_open_file (igribout,    TRIM(outfile),  cmode, ilocerr)
    IF (ilocerr /= GRIB_SUCCESS) THEN
      errmsg(:) = ' '
      CALL codes_get_error_string(ilocerr, errmsg)
      errmsg = 'ERROR '//TRIM(yzroutine)//' opening '//TRIM(outfile)//': '//TRIM(errmsg)
      error = 4
      RETURN
    END IF

    ALLOCATE(comp_tot_grb(cmp_meta%ni*cmp_meta%nj))
    comp_tot_grb = miss_value

    ij = 0
    DO j=1, cmp_meta%nj
      DO i=1, cmp_meta%ni
        ij = i + (j-1)*cmp_meta%ni
        comp_tot_grb(ij) = REAL(comp2d_tot(i,j), kind=KIND(1.0))
      END DO
    END DO
    CALL codes_clone(igribsample, igribclone) ! clone sample before modifying it

    ierr = 0

    ! .. This adds model specific infos like run-type, ens-information etc., with own error checking:
    CALL grib2_add_modelspec_info(idom, grib_locinfo, ierr)
    ierr = ierr + ilocerr


    ! local section has to be deleted first before changing the centre:

    CALL codes_set  (igribclone, 'tablesVersion',    grib_locinfo%tablesVersion, ilocerr)
      CALL check_codes_err(ilocerr, 'tablesVersion', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'typeOfProcessedData',    grib_locinfo%typeOfProcessedData, ilocerr)
      CALL check_codes_err(ilocerr, 'typeOfProcessedData', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'centre',    grib_locinfo%generatingCenter, ilocerr)
      CALL check_codes_err(ilocerr, 'centre', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'subCentre',    grib_locinfo%generatingSubcenter, ilocerr)
      CALL check_codes_err(ilocerr, 'subCentre', TRIM(yzroutine), increrr=ierr)

    ! local section has to be deleted, to avoid problems coming from the sample file. Then it can be re-created:
    CALL codes_set (igribclone, 'setLocalDefinition', GRIB_FALSE, ierr)   ! delete local section
      CALL check_codes_err(ilocerr, 'setLocalDefinition', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'setLocalDefinition',  GRIB_TRUE, ierr)   ! activate local section again
      CALL check_codes_err(ilocerr, 'setLocalDefinition', TRIM(yzroutine), increrr=ierr)

    CALL codes_set  (igribclone, 'localDefinitionNumber',    grib_locinfo%localDefinitionNumber, ilocerr)
      CALL check_codes_err(ilocerr, 'localDefinitionNumber', TRIM(yzroutine), increrr=ierr)

    ! for certain cases (e.g. cdo-remapped model data) we found issues with the local section (missing keys,
    !   also in the template clone). these issues likely depend on the setting of localDefinitionNumber.
    ! to circumvent them, we check for localDefinitionNumber value and in case it is set to the
    !   undetermined-flag value (-1), we skip the setting & error checking for the effected keys (no
    !   guarantee that no issues with the written file will occur downstream, e.g. in postprocessing).
    IF (grib_locinfo%localDefinitionNumber .GT. 0) THEN
      CALL codes_set  (igribclone, 'localNumberOfExperiment',    grib_locinfo%localNumberOfExperiment, ilocerr)
        CALL check_codes_err(ilocerr, 'localNumberOfExperiment', TRIM(yzroutine), increrr=ierr)

      CALL codes_set  (igribclone, 'localCreationDateYear', localCreationDateYear, ilocerr)
        CALL check_codes_err(ilocerr, 'localCreationDateYear', TRIM(yzroutine), increrr=ierr)

      CALL codes_set  (igribclone, 'localCreationDateMonth', localCreationDateMonth, ilocerr)
        CALL check_codes_err(ilocerr, 'localCreationDateMonth', TRIM(yzroutine), increrr=ierr)

      CALL codes_set  (igribclone, 'localCreationDateDay', localCreationDateDay, ilocerr)
        CALL check_codes_err(ilocerr, 'localCreationDateDay', TRIM(yzroutine), increrr=ierr)

      CALL codes_set  (igribclone, 'localCreationDateHour', localCreationDateHour, ilocerr)
        CALL check_codes_err(ilocerr, 'localCreationDateHour', TRIM(yzroutine), increrr=ierr)

      CALL codes_set  (igribclone, 'localCreationDateMinute', localCreationDateMinute, ilocerr)
        CALL check_codes_err(ilocerr, 'localCreationDateMinute', TRIM(yzroutine), increrr=ierr)

      CALL codes_set  (igribclone, 'localCreationDateSecond', localCreationDateSecond, ilocerr)
        CALL check_codes_err(ilocerr, 'localCreationDateSecond', TRIM(yzroutine), increrr=ierr)
    ELSE
      WRITE(*,*) 'WARNING: Grib keys in Local Section will be incomplete.'
      WRITE(*,*) '         Skipping: localNumberOfExperiment, localCreationDate*'
    END IF

    IF (grib_locinfo%typeOfGeneratingProcess == TGP_ENSEMBLE) THEN
      CALL codes_set  (igribclone, 'localTypeOfEnsembleForecast',    grib_locinfo%localTypeOfEnsembleForecast, ilocerr)
        CALL check_codes_err(ilocerr, 'localTypeOfEnsembleForecast', TRIM(yzroutine), increrr=ierr)
    END IF
    CALL codes_set  (igribclone, 'productionStatusOfProcessedData',    grib_locinfo%productionStatusOfProcessedData, ilocerr)
      CALL check_codes_err(ilocerr, 'productionStatusOfProcessedData', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'generatingProcessIdentifier',    grib_locinfo%generatingProcessIdentifier, ilocerr)
      CALL check_codes_err(ilocerr, 'generatingProcessIdentifier', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'backgroundProcess',    grib_locinfo%backgroundProcess, ilocerr)
      CALL check_codes_err(ilocerr, 'backgroundProcess', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'typeOfGeneratingProcess',    grib_locinfo%typeOfGeneratingProcess, ilocerr)
      CALL check_codes_err(ilocerr, 'typeOfGeneratingProcess', TRIM(yzroutine), increrr=ierr)

    IF (grib_locinfo%typeOfGeneratingProcess == TGP_ENSEMBLE) THEN
      CALL codes_set  (igribclone, 'productDefinitionTemplateNumber',    1, ilocerr)  ! 1 = individual ensemble forecast
        CALL check_codes_err(ilocerr, 'productDefinitionTemplateNumber', TRIM(yzroutine), increrr=ierr)
      CALL codes_set  (igribclone, 'typeOfEnsembleForecast',    grib_locinfo%typeOfEnsembleForecast, ilocerr)
        CALL check_codes_err(ilocerr, 'typeOfEnsembleForecast', TRIM(yzroutine), increrr=ierr)
      CALL codes_set  (igribclone, 'numberOfForecastsInEnsemble',    grib_locinfo%numberOfForecastsInEnsemble, ilocerr)
        CALL check_codes_err(ilocerr, 'numberOfForecastsInEnsemble', TRIM(yzroutine), increrr=ierr)
      CALL codes_set  (igribclone, 'perturbationNumber',    grib_locinfo%perturbationNumber, ilocerr)
        CALL check_codes_err(ilocerr, 'perturbationNumber', TRIM(yzroutine), increrr=ierr)
    ELSE
      CALL codes_set  (igribclone, 'productDefinitionTemplateNumber',    0, ilocerr)  ! 0 = standard products
        CALL check_codes_err(ilocerr, 'productDefinitionTemplateNumber', TRIM(yzroutine), increrr=ierr)
    END IF

    CALL codes_set  (igribclone, 'Ni',    cmp_meta%ni, ilocerr)  ! this also adjusts the size of memory in the sample
      CALL check_codes_err(ilocerr, 'Ni', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone, 'Nj',    cmp_meta%nj, ilocerr)  ! this also adjusts the size of memory in the sample
      CALL check_codes_err(ilocerr, 'Nj', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone,'ijDirectionIncrementGiven',       1)
      CALL check_codes_err(ilocerr ,'ijDirectionIncrementGiven', TRIM(yzroutine), increrr=ierr)


    CALL codes_set (igribclone, 'typeOfLevel',  'radarElevComposite', ilocerr)
      CALL check_codes_err(ilocerr, 'typeOfLevel', TRIM(yzroutine), increrr=ierr)

    ! .. Time coding:
    SELECT CASE (TRIM(zcomptyp))
    CASE ('obs')
      CALL codes_set (igribclone, 'significanceOfReferenceTime',     3, ilocerr)   ! 3= observation time
        CALL check_codes_err(ilocerr, 'significanceOfReferenceTime', TRIM(yzroutine), increrr=ierr)
      CALL codes_set (igribclone,'indicatorOfUnitOfTimeRange',       0, ilocerr)   ! 0=minutes in grib2
        CALL check_codes_err(ilocerr, 'indicatorOfUnitOfTimeRange', TRIM(yzroutine), increrr=ierr)

      READ(model_validtime( 1: 8),'(I8)') idate
      READ(model_validtime( 9:12),'(I4)') itime
      READ(model_validtime(13:14),'(I2)') isecond
      CALL codes_set (igribclone,'dataDate',                      idate, ilocerr)   ! yyyymmdd
        CALL check_codes_err(ilocerr,'dataDate', TRIM(yzroutine), increrr=ierr)
      CALL codes_set (igribclone,'dataTime',                      itime, ilocerr)   ! hhmm
        CALL check_codes_err(ilocerr,'dataTime', TRIM(yzroutine), increrr=ierr)
      CALL codes_set (igribclone,'second',                      isecond, ilocerr)   ! ss
        CALL check_codes_err(ilocerr,'second', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'forecastTime',0,ilocerr)
        CALL check_codes_err(ilocerr, 'forecastTime', TRIM(yzroutine), increrr=ierr)

    CASE ('sim')
      ! We code the time information in the way of a forecast. This is not strictly correct for analyses.
      CALL codes_set (igribclone, 'significanceOfReferenceTime',     1, ilocerr)   ! 1= start of forecast; 0= analysis time
        CALL check_codes_err(ilocerr, 'significanceOfReferenceTime', TRIM(yzroutine), increrr=ierr)
      CALL codes_set (igribclone,'indicatorOfUnitOfTimeRange',       0, ilocerr)   ! 0=minutes in grib2
        CALL check_codes_err(ilocerr, 'indicatorOfUnitOfTimeRange', TRIM(yzroutine), increrr=ierr)

      READ(model_starttime( 1: 8),'(I8)') idate
      READ(model_starttime( 9:12),'(I4)') itime
      READ(model_starttime(13:14),'(I2)') isecond
      CALL codes_set (igribclone,'dataDate',                      idate, ilocerr)   ! yyyymmdd
        CALL check_codes_err(ilocerr,'dataDate', TRIM(yzroutine), increrr=ierr)
      CALL codes_set (igribclone,'dataTime',                      itime, ilocerr)   ! hhmm
        CALL check_codes_err(ilocerr,'dataTime', TRIM(yzroutine), increrr=ierr)
      CALL codes_set (igribclone,'second',                      isecond, ilocerr)   ! ss
        CALL check_codes_err(ilocerr,'second', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'forecastTime',forecastminutes,ilocerr)
        CALL check_codes_err(ilocerr, 'forecastTime', TRIM(yzroutine), increrr=ierr)
      
    END SELECT

    CALL codes_set (igribclone, 'latitudeOfSouthernPoleInDegrees',  -cmp_meta%pollat, ilocerr)
      CALL check_codes_err(ilocerr, 'latitudeOfSouthernPoleInDegrees', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'longitudeOfSouthernPoleInDegrees',  cmp_meta%pollon + 180.0_dp, ilocerr)
      CALL check_codes_err(ilocerr, 'longitudeOfSouthernPoleInDegrees', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'angleOfRotationInDegrees',          cmp_meta%polgam, ilocerr)
      CALL check_codes_err(ilocerr, 'angleOfRotationInDegrees', TRIM(yzroutine), increrr=ierr)

    CALL codes_set (igribclone, 'iDirectionIncrementInDegrees',            cmp_meta%dlon, ilocerr)
      CALL check_codes_err(ilocerr, 'iDirectionIncrementInDegrees', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'jDirectionIncrementInDegrees',            cmp_meta%dlat, ilocerr)
      CALL check_codes_err(ilocerr, 'jDirectionIncrementInDegrees', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'longitudeOfFirstGridPointInDegrees',  cmp_meta%startlon, ilocerr)
      CALL check_codes_err(ilocerr, 'longitudeOfFirstGridPointInDegrees', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'latitudeOfFirstGridPointInDegrees',   cmp_meta%startlat, ilocerr)
      CALL check_codes_err(ilocerr, 'latitudeOfFirstGridPointInDegrees', TRIM(yzroutine), increrr=ierr)
    CALL codes_set (igribclone, 'longitudeOfLastGridPointInDegrees',   &
         cmp_meta%startlon+(cmp_meta%ni-1)*cmp_meta%dlon, ilocerr)
      CALL check_codes_err(ilocerr, 'longitudeOfLastGridPointInDegrees', TRIM(yzroutine), increrr=ierr)

    CALL codes_set (igribclone, 'latitudeOfLastGridPointInDegrees',    &
         cmp_meta%startlat+(cmp_meta%nj-1)*cmp_meta%dlat, ilocerr)
      CALL check_codes_err(ilocerr, 'latitudeOfLastGridPointInDegrees', TRIM(yzroutine), increrr=ierr)

    CALL codes_set  (igribclone, 'level', ilevel, ilocerr) ! 0     = composite for bubble generator
                                                          ! 1...9 = 9 different elevations depending on namelist settings
                                                          ! 98    = composite of the DWD precipitation scans
                                                          ! 99    = maximum composite using all elevations
      CALL check_codes_err(ilocerr, 'level', TRIM(yzroutine), increrr=ierr)

    CALL codes_set  (igribclone,'scaleFactorOfFirstFixedSurface',      0, ilocerr)
      CALL check_codes_err(ilocerr,'scaleFactorOfFirstFixedSurface', TRIM(yzroutine), increrr=ierr)
    CALL codes_set  (igribclone,'scaledValueOfFirstFixedSurface', ilevel, ilocerr)
      CALL check_codes_err(ilocerr,'scaledValueOfFirstFixedSurface', TRIM(yzroutine), increrr=ierr)

    ! Set Shortname and (re-)set 'typeOfGeneratingProcess' for observations:
    SELECT CASE (TRIM(zcomptyp))
    CASE ('obs')
      CALL codes_set  (igribclone,'shortName', 'DBZCMP_OBS', ilocerr)
        CALL check_codes_err(ilocerr,'shortName=DBZCMP_OBS', TRIM(yzroutine), increrr=ierr)
      CALL codes_set  (igribclone,'typeOfGeneratingProcess',    8, ilocerr)
        CALL check_codes_err(ilocerr,'typeOfGeneratingProcess', TRIM(yzroutine), increrr=ierr)
    CASE ('sim')
      CALL codes_set  (igribclone,'shortName', 'DBZCMP_SIM', ilocerr)
        CALL check_codes_err(ilocerr,'shortName=DBZCMP_SIM', TRIM(yzroutine), increrr=ierr)
    END SELECT

    ! Set bit resolution:
    CALL codes_set  (igribclone, 'bitsPerValue',       16, ilocerr)  ! 16=value from sample, but can be changed here (e.g., 24)!
      CALL check_codes_err(ilocerr, 'bitsPerValue', TRIM(yzroutine), increrr=ierr)

    ! Put data to "values":
    CALL codes_set  (igribclone, 'values', comp_tot_grb(:), ilocerr)
      CALL check_codes_err(ilocerr, 'values', TRIM(yzroutine), increrr=ierr)

    ! Set grib packing/compression method: "grid_simple", "grid_ccsds", "png"
    CALL codes_set(igribclone, 'packingType', TRIM(comp_grib2_packingtype), ilocerr)
      CALL check_codes_err(ilocerr, 'packingType', TRIM(yzroutine), increrr=ierr)

    IF (ierr /= GRIB_SUCCESS) THEN
      error = 4
      errmsg(:) = ' '
      errmsg = 'ERROR '//TRIM(yzroutine)//': codes_set failed for some keys or the values. Check earlier ERROR messages above!'
    ELSE
      CALL codes_write(igribclone, igribout, ilocerr)
      IF (ilocerr /= GRIB_SUCCESS) THEN
        errmsg(:) = ' '
        CALL codes_get_error_string(ilocerr, errmsg)
        errmsg = 'ERROR '//TRIM(yzroutine)//' writing grib: '//TRIM(errmsg)
        error = 5
      END IF
    END IF
    CALL codes_release(igribclone)

    DEALLOCATE (comp_tot_grb)
    CALL codes_release   (igribsample)
    CALL codes_close_file(igribout   , ilocerr)
    IF (ilocerr /= GRIB_SUCCESS) THEN
      errmsg(:) = ' '
      CALL codes_get_error_string(ilocerr, errmsg)
      errmsg = 'ERROR '//TRIM(yzroutine)//' closing igribout: '//TRIM(errmsg)
      error = 6
      RETURN
    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE write_composite_grib

#endif
  
  !=========================================================================
  !
  ! Procedure for creating output file names for radar composite files
  ! based on a pattern given in the file_pattern on input
  !
  !=========================================================================
  
  SUBROUTINE composite_create_filename (file_pattern, time_mod, cfieldname, comptype, ofilename)

    CHARACTER(len=*), INTENT(in)         :: file_pattern ! parameter name
    REAL(KIND=dp)   , INTENT(in)         :: time_mod   ! time in seconds since model start
    CHARACTER(len=*), INTENT(in)         :: cfieldname ! parameter name
    CHARACTER(len=*), INTENT(in)         :: comptype   ! 'obs' or 'sim'
    CHARACTER(len=*), INTENT(out)        :: ofilename  ! output file name

    CHARACTER(len=*), PARAMETER :: yzroutine = 'emvorado::composite_create_filename'
    CHARACTER(len=14) :: model_starttime, model_validtime
    CHARACTER(len=8)  :: ddhhmmss_validtime
    CHARACTER(len=LEN(file_pattern)) :: fpat
    CHARACTER(len=20) :: key
    INTEGER           :: ku, ko, pos
    LOGICAL           :: have_time


    ofilename(:) = ' '

    !--------------------------------------------------------------------------
    ! .. Preparation of various time stamp formats for the filename pattern:

    ddhhmmss_validtime = get_model_time_ddhhmmss( time_mod )
    model_starttime    = get_datetime_ini()
    model_validtime    = new_datetime(model_starttime, time_mod)
    
    IF (LEN_TRIM(file_pattern) <= 0) THEN

      !--------------------------------------------------------------------------
      ! Set the default file name pattern:

      ofilename = TRIM(cfieldname)//'_'//TRIM(comptype)//'_'//TRIM(model_starttime)//'_'//TRIM(model_validtime)//'.grb2'

    ELSE

      !--------------------------------------------------------------------------
      ! Set the user defined file name pattern

      ! .. List of valid <keys> for filename patterns:
      !
      !       - <varname>      : will result in cfieldname (optional)
      !       - <simobs>       : will result in 'sim' or 'obs', depending on comptype (optional)
      !       - <tmodelini>    : will result in YYYYMMDDhhmmss  (absolute datetime of model start) (optional)
      !       - <tact>         : will result in YYYYMMDDhhmmss  (absolute datetime of actual model time) \
      !       - <tvvzact>      : will result in DDhhmmss  (forecast lead time of actual model time)      / (one of both is mandatory)
      
      fpat = file_pattern ! copy original pattern

      have_time     = .FALSE.
      
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. and, if it is a time key, replace by its actual value:
        SELECT CASE (TRIM(key))
        CASE ('<tmodelini>')
          ! optional key, therefore no have_time = .TRUE.:
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(model_starttime))
        CASE ('<tact>')
          have_time = .TRUE.
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(model_validtime))
        CASE ('<tvvzact>')
          have_time = .TRUE.
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(ddhhmmss_validtime))
        CASE ('<varname>')
          ! optional key:
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(cfieldname))
        CASE ('<simobs>')
          ! optional key:
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(comptype))
        END SELECT
      END DO

      ! At last, check if any other keys or "<" or ">" are still remaining in the fpat. If yes, it is a wrong key:
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. if there is still a key left in fpat, it is a wrong key:
        CALL abort_run (my_radar_id, 29074, &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" contains unknown key '//TRIM(key), &
             TRIM(yzroutine)//', input parameter ''file_pattern''')
      END DO

      ! .. if still some remaining single '<' or '>' characters are in fpat, there is an error:
      IF (INDEX(fpat,'<') > 0 .OR. INDEX(fpat,'>') > 0) THEN
        CALL abort_run (my_radar_id, 29075, &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" contains garbled "<" or ">" characters', &
             TRIM(yzroutine)//', input parameter ''file_pattern''')
      END IF

      ! .. check if all mandatory keys have been specified in file_pattern:
      IF (.NOT.have_time) THEN
        CALL abort_run (my_radar_id, 29076, &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" does not contain a <key> for actual time (<tact>, <tvvzact>)', &
             TRIM(yzroutine)//', input parameter ''file_pattern''')
      END IF

      ! .. at this stage, fpat contains a valid filename:
      ofilename = TRIM(fpat)
      
    END IF

  END SUBROUTINE composite_create_filename

END MODULE radar_output_composite
