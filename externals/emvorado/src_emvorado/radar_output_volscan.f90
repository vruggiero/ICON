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


MODULE radar_output_volscan

!------------------------------------------------------------------------------
!
! Description: Core procedures of the radar forward operator EMVORADO for writing
!              volume data into netcdf (cdfin) or grib2 files.
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
  
  USE radar_data, ONLY :  miss_value, zero_value

  USE radar_interface, ONLY : &
       abort_run,          &
       grib2_add_modelspec_info

  USE radar_utilities, ONLY :  new_datetime, diff_seconds

  USE radar_output_utils, ONLY : eps_ascii, check_nc, check_codes_err, to_lower_or_upper
   
  USE radar_data, ONLY : &
       cmaxlen, &
       radar_meta_type, &
       ydate_ini_mod, &
       idom, &
       my_radar_id

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, &
       cdfin_creator_model, labort_if_problems_gribout

  USE radar_obs_meta_list, ONLY : get_elarr_precipscan

  USE radar_data_io, ONLY : PRODUCT_TYPE_OBS, PRODUCT_TYPE_SIM, &
       &                    t_scaledMetadataList, t_grib, t_grib2_modelspec, t_grib_loc_sec_tab, &
       &                    MAX_LEN_TAB_ENTRY, MAX_LEN_DESC, LEN_DATETIME, &
       &                    TGP_ENSEMBLE, TGP_ANALYSIS, TGP_FORECAST, TGP_INITRUN, &
       &                    BITMAP_PRI_NOBITMAP, BITMAP_PRI_MISSVAL, &
       &                    BITMAP_PRI_UNDEVAL, BITMAP_PRI_MISSVAL_UNDEVAL, &
       &                    GRIB_MISSING, GRIB_TRUE, GRIB_FALSE, GRIB_EDITION_2, &
       &                    INT_ATT, REAL_ATT, CHAR_ATT, t_cdfin_globalatt, NATTR_MAX

  !------------------------------------------------------------------------------

#ifdef NETCDF
  USE netcdf, ONLY :  &
       nf90_open, &
       NF90_NOERR, &
       NF90_CLOBBER, &
       NF90_WRITE, &
       NF90_NOWRITE, &
       NF90_UNLIMITED, &
       nf90_create, &
       nf90_inquire, &
       nf90_inquire_dimension, &
       nf90_def_dim, &
       nf90_inq_dimid , &
       nf90_inq_varid , &
       nf90_get_var , &
       nf90_close, &
       nf90_put_att, &
       nf90_get_att, &
       nf90_redef, &
       nf90_strerror, &
       NF90_INT, &
       nf90_global, &
       nf90_enddef, &
       nf90_put_var, &
       NF90_CHAR, &
       nf90_def_var, &
       nf90_def_var_deflate, &
       NF90_NETCDF4, &
       NF90_CLASSIC_MODEL, &
       NF90_REAL
#endif

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

#ifdef GRIBAPI
  INTERFACE get_scaling
    MODULE PROCEDURE       &
         get_scaling_real, &
         get_scaling_int
  END INTERFACE get_scaling
#endif

  !==============================================================================
  ! Public and Private Subroutines

  PRIVATE

  PUBLIC :: prepare_write_volscan

#ifdef NETCDF
  PUBLIC :: write_cdfin_from_volscan
#endif
    
#ifdef GRIBAPI
  PUBLIC :: write_grib_from_volscan, fill_grib_input_volscan, write_grib_to_file, get_loc_sec_tab
#endif

  !==============================================================================
  ! Module variables

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  ! --------------------------------------------------------------------------------
  ! ----------------------- Writer for volscans in NETCDF --------------------------
  ! --------------------------------------------------------------------------------

#ifdef NETCDF

  SUBROUTINE write_cdfin_from_volscan(rsm, itime, grib_input, grib_locinfo, &
       &                              ncdf_filename, force_overwrite, ncdf_shortname, &
       &                              config_string, unit_string, field_desc, dat, &
       &                              ind_ele_outlist, nele_out)

    !------------------------------------------------------------------------------
    !
    ! Description: creates netcdf file for modeled output data. Similar format
    !  as DWD cdfin, except for the names of dimensions/variables.
    ! No quality flags yet.
    !
    ! Call prepare_write_volscan() before!
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! INCOMING:
    TYPE(radar_meta_type)  , INTENT(IN)   :: rsm
    INTEGER                , INTENT(in)   :: itime
    TYPE(t_grib)           , INTENT(IN)   :: grib_input             ! Take only date/time info + product_type from that
    TYPE(t_grib2_modelspec), INTENT(IN)   :: grib_locinfo
    REAL(KIND=dp)          , INTENT(IN)   :: dat(:,:,:)             ! Radar data 3D (naz,nra,nel)
    CHARACTER (len=*)      , INTENT(in)   :: ncdf_filename          ! NetCDF file name
    LOGICAL                , INTENT(in)   :: force_overwrite        ! if .true., forcefully overwrite any pre-existing output file of name "grib2_filename"
    CHARACTER (len=*)      , INTENT(in)   :: ncdf_shortname         ! NetCDF short name of radar data field
    CHARACTER (len=*)      , INTENT(in)   :: config_string          ! Description of emvorado config
    CHARACTER (len=*)      , INTENT(in)   :: unit_string            ! Unit of data field (vari attribute)
    CHARACTER (len=*)      , INTENT(in)   :: field_desc             ! Description of data field (vari attribute)
    INTEGER                , INTENT(in)   :: ind_ele_outlist(:)
    INTEGER                , INTENT(in)   :: nele_out               ! length of ind_ele_outlist

    ! Local scalars:
    CHARACTER(len=*), PARAMETER :: yzroutine = 'write_cdfin_from_volscan'
    LOGICAL                     :: ncdf_file_exists             ! Check if the file exists
    INTEGER                     :: fid                          ! file IDs

    TYPE(t_cdfin_globalatt)     :: glob_attrib                  ! Container for global attributes

    INTEGER                     :: status, nobs
    INTEGER                     :: cmode


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE netcdf_create_dwd_from_mod
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    fid = -1

    IF (nele_out <= 0) RETURN

    ! Fill type glob_attrib with informations. This will allocate glob_attrib%att(NATTR_MAX):
    CALL fill_cdfin_global_attributes (rsm, grib_input, grib_locinfo, TRIM(config_string), glob_attrib)

    !------------------------------------------------------------------------------
    ! Section 2: Opening/Creating the file
    !------------------------------------------------------------------------------

    INQUIRE (file=TRIM(ncdf_filename), exist=ncdf_file_exists)

    IF ( .NOT. ncdf_file_exists .OR. force_overwrite ) THEN ! First call for this file or the name changes every call; file is created or old file overwritten

      ! Open and create File
      ! ----------------------

      cmode = NF90_CLOBBER
      cmode = IOR(cmode, NF90_NETCDF4)
      cmode = IOR(cmode, NF90_CLASSIC_MODEL)  ! no fancy netcdf4 features (new data types, compounds, multiple unlimited dims, etc.), therefore smaller files

!!$ apparently not possible together with NF90_NETCDF4:     cmode = IOR(cmode, NF90_64BIT_OFFSET)
      status = check_nc( nf90_create (TRIM(ncdf_filename), cmode, fid), 'nf90_create '//TRIM(ncdf_filename) )

      IF (status /= nf90_noerr) THEN
        WRITE(*,'(i5,1x,3a)') my_radar_id, 'ERROR creating NetCDF file ', TRIM(ncdf_filename)//' ', &
             TRIM(nf90_strerror(status))
        fid = -1
        RETURN
      END IF

      CALL cdfin_define (fid, rsm, glob_attrib, ncdf_shortname, ind_ele_outlist, nele_out)

    ELSE

      ! Successive calls just open the file
      ! -----------------------------------

      IF (ldebug_radsim) WRITE (*,'(A)') TRIM(yzroutine)//': opening '//TRIM(ncdf_filename)//' ...'

      status =  check_nc( nf90_open (TRIM(ncdf_filename), NF90_WRITE, fid), 'nf90_open '//TRIM(ncdf_filename) )

      IF (status /= nf90_noerr) THEN
        WRITE(*,'(i5,1x,3a)') my_radar_id, 'ERROR opening NetCDF file ', TRIM(ncdf_filename)//' ', &
             TRIM(nf90_strerror(status))
        fid = -1
        RETURN
      END IF


    END IF

    !------------------------------------------------------------------------------
    ! Section 3: Writing/appending data to file
    !------------------------------------------------------------------------------

    CALL cdfin_write_mod_data (fid, rsm, itime, ncdf_shortname, unit_string, field_desc, glob_attrib%product_type, &
         &                     TRIM(grib_input%dateTime_act), grib_input%forecastseconds, grib_input%forecastminutes, &
         &                     dat, ind_ele_outlist, nele_out)

    !------------------------------------------------------------------------------
    ! Section 4: Close file
    !------------------------------------------------------------------------------

    status = check_nc( nf90_close(fid), 'nf90_close '//TRIM(ncdf_filename) )
    IF (status /= nf90_noerr) WRITE (*,'(A)') "Error in closing file "//TRIM(ncdf_filename)// &
         " . Status error:  "//TRIM(nf90_strerror(status))

    !------------------------------------------------------------------------------
    ! Section 5: Clean up
    !------------------------------------------------------------------------------

    IF (ALLOCATED(glob_attrib%att)) THEN
      DEALLOCATE(glob_attrib%att)
      glob_attrib%n_attr = 0
    END IF
    
    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE write_cdfin_from_volscan

! ------------------------------------------------------------------------------

  SUBROUTINE cdfin_define (fid, rsm, glob_att, ncdf_shortname, ind_ele_outlist, nele_out)

    !------------------------------------------------------------------------------
    !
    ! INCOMING:
    INTEGER                , INTENT(in) :: fid                 ! file ID
    TYPE(radar_meta_type)  , INTENT(IN) :: rsm
    TYPE(t_cdfin_globalatt), INTENT(IN) :: glob_att
    CHARACTER (len=*)      , INTENT(in) :: ncdf_shortname      ! NetCDF short name of radar data field
    INTEGER                , INTENT(in) :: ind_ele_outlist(:)  ! List of elevation indices for output
    INTEGER                , INTENT(in) :: nele_out            ! length of ind_ele_outlist(:)

    ! Local variables
    INTEGER                  :: i, status, dummy
    ! Idenfidiers of dimensions id
    INTEGER                  :: az_dimid, ra_dimid, el_dimid, rec_dimid, strlen_dimid
    ! Metadata identifiers
    INTEGER                  :: varid_yssosn  ! Short station name
    INTEGER                  :: varid_mii     ! Station ID country
    INTEGER                  :: varid_niii    ! Local station ID
    INTEGER                  :: varid_mlah    ! Latitude of station
    INTEGER                  :: varid_mloh    ! Longitude of station
    INTEGER                  :: varid_mhp     ! height of station
    INTEGER                  :: varid_rngbs   ! range bin size (021201)
    INTEGER                  :: varid_azres   ! azimutal resolution (021202)
    INTEGER                  :: varid_rgoff   ! offset of radial (021203)
    INTEGER                  :: varid_azoff   ! azimut start (center of sampling interval)
    INTEGER                  :: varid_enyq    ! extended NYQUIST Velocity (021236)
    INTEGER                  :: varid_hnyq    ! high NYQUIST Velocity (021237)
    INTEGER                  :: varid_dprf    ! DualPRF Ratio (002194)
    INTEGER                  :: varid_mrglen  ! RANGE-GATE LENGTH
    INTEGER                  :: varid_mngave  ! NUMBER OF GATES AVERAGED
    INTEGER                  :: varid_mnipul  ! NUMBER OF INTEGRATED PULSES
    INTEGER                  :: varid_mjjj    ! Year
    INTEGER                  :: varid_mmm     ! Month
    INTEGER                  :: varid_myy     ! Day
    INTEGER                  :: varid_mgg     ! Hour
    INTEGER                  :: varid_ngg     ! Minute
    INTEGER                  :: varid_nsecm   ! Second within a minute
    INTEGER                  :: varid_vvzs    ! Forecast time in seconds
    INTEGER                  :: varid_vvzm    ! Forecast time in minutes
    INTEGER                  :: varid_time    ! observation time
    INTEGER                  :: varid_date    ! reference date
    INTEGER                  :: varid_mabaz   ! azimuth mabaz (end of sampling interval)
    INTEGER                  :: varid_mabaz0  ! azimuth mabaz0 (end of sampling interval)
    INTEGER                  :: varid_manel   ! antenna elevation
    INTEGER                  :: varid_manel0  ! elevation manel
    INTEGER                  :: varid_rbins   ! number of range bins of individual report

    ! Fill values for some variables
    REAL,    PARAMETER    :: missing_val_real =  HUGE(1.0)
    INTEGER, PARAMETER    :: missing_val_int  = -HUGE(1)
    INTEGER               :: chunksizes1D(1), chunksizes2D(2), chunksizes3D(3)

    CHARACTER(LEN=*), PARAMETER :: yzroutine = 'cdfin_define'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//': creating file ...'

    ! Define dimensions:
    !-------------------
    status = check_nc( nf90_def_dim(fid, "n_azimuth",   rsm%naz, az_dimid), 'nf90_def_dim n_azimuth')  ! Azimuhtal dimension
    status = check_nc( nf90_def_dim(fid, "n_range",     rsm%nra, ra_dimid), 'nf90_def_dim n_range')  ! Range dimension
    status = check_nc( nf90_def_dim(fid, "n_elevation", nele_out, el_dimid), 'nf90_def_dim n_elevation')  ! Elevations dimension
    status = check_nc( nf90_def_dim(fid, "strlen",      3,   strlen_dimid), 'nf90_def_dim strlen')  ! String length dimension
    status = check_nc( nf90_def_dim(fid, "records", NF90_UNLIMITED, rec_dimid), 'nf90_def_dim records') ! Records dimension


    ! Define variables:
    !------------------

    chunksizes1D(1) = nele_out    ! For variables having only the unlimited dimension,
                                ! the default chunksize in netcdf4 is 1048576 and creates horribly large files!
                                ! For our case, a value of rsm%nel resp. nele_out or even 1 is much more appropriate!
    chunksizes2D(:) = (/rsm%naz, nele_out/)

    ! 0-dim variables (these are 1d but constant)
    status = check_nc( nf90_def_var(fid, 'station_name',          NF90_CHAR, (/strlen_dimid,rec_dimid/), &
         varid_yssosn), 'nf90_def_var station_name')
    status = check_nc( nf90_def_var(fid, 'country_ID',            NF90_INT,  (/rec_dimid/), &
         varid_mii,    chunksizes=chunksizes1D), 'nf90_def_var country_ID')
    status = check_nc( nf90_def_var(fid, 'station_ID_national',   NF90_INT,  (/rec_dimid/), &
         varid_niii,   chunksizes=chunksizes1D), 'nf90_def_var station_ID_national')
    status = check_nc( nf90_def_var(fid, 'station_latitude',      NF90_REAL, (/rec_dimid/), &
         varid_mlah,   chunksizes=chunksizes1D), 'nf90_def_var station_latitude')
    status = check_nc( nf90_def_var(fid, 'station_longitude',     NF90_REAL, (/rec_dimid/), &
         varid_mloh,   chunksizes=chunksizes1D), 'nf90_def_var station_longitude')
    status = check_nc( nf90_def_var(fid, 'station_height',        NF90_REAL, (/rec_dimid/), &
         varid_mhp,    chunksizes=chunksizes1D), 'nf90_def_var station_height')
    status = check_nc( nf90_def_var(fid, 'range_resolution',      NF90_REAL, (/rec_dimid/), &
         varid_rngbs,  chunksizes=chunksizes1D), 'nf90_def_var range_resolution')
    status = check_nc( nf90_def_var(fid, 'azimuthal_resolution',  NF90_REAL, (/rec_dimid/), &
         varid_azres,  chunksizes=chunksizes1D), 'nf90_def_var azimuthal_resolution')
    status = check_nc( nf90_def_var(fid, 'range_start',           NF90_REAL, (/rec_dimid/), &
         varid_rgoff,  chunksizes=chunksizes1D), 'nf90_def_var range_start')
    status = check_nc( nf90_def_var(fid, 'azimuth_start',         NF90_REAL, (/rec_dimid/), &
         varid_azoff,  chunksizes=chunksizes1D), 'nf90_def_var azimuth_start')
    status = check_nc( nf90_def_var(fid, 'extended_nyquist',      NF90_REAL, (/rec_dimid/), &
         varid_enyq,   chunksizes=chunksizes1D), 'nf90_def_var extended_nyquist')
    status = check_nc( nf90_def_var(fid, 'high_nyquist',          NF90_REAL, (/rec_dimid/), &
         varid_hnyq,   chunksizes=chunksizes1D), 'nf90_def_var high_nyquist')
    status = check_nc( nf90_def_var(fid, 'dualPRF_ratio',         NF90_REAL, (/rec_dimid/), &
         varid_dprf,   chunksizes=chunksizes1D), 'nf90_def_var dualPRF_ratio')
    status = check_nc( nf90_def_var(fid, 'range_gate_length',     NF90_REAL, (/rec_dimid/), &
         varid_mrglen, chunksizes=chunksizes1D), 'nf90_def_var range_gate_length')
    status = check_nc( nf90_def_var(fid, 'n_ranges_averaged',     NF90_INT,  (/rec_dimid/), &
         varid_mngave, chunksizes=chunksizes1D), 'nf90_def_var n_ranges_averaged')
    status = check_nc( nf90_def_var(fid, 'n_pulses_averaged',     NF90_INT,  (/rec_dimid/), &
         varid_mnipul, chunksizes=chunksizes1D), 'nf90_def_var n_pulses_averaged')

    ! 1-dim variables
    status = check_nc( nf90_def_var (fid, 'year',  NF90_INT,  (/rec_dimid/), &    ! Year
         varid_mjjj,  chunksizes=chunksizes1D), 'nf90_def_var year')
    status = check_nc( nf90_def_var (fid, 'month', NF90_INT,  (/rec_dimid/), &    ! Month
         varid_mmm,   chunksizes=chunksizes1D), 'nf90_def_var month')
    status = check_nc( nf90_def_var (fid, 'day',   NF90_INT,  (/rec_dimid/), &    ! Day
         varid_myy,   chunksizes=chunksizes1D), 'nf90_def_var day')
    status = check_nc( nf90_def_var (fid, 'hour',  NF90_INT,  (/rec_dimid/), &    ! Hour
         varid_mgg,   chunksizes=chunksizes1D), 'nf90_def_var hour')
    status = check_nc( nf90_def_var (fid, 'minute',NF90_INT,  (/rec_dimid/), &    ! Minute
         varid_ngg,   chunksizes=chunksizes1D), 'nf90_def_var minute')
    status = check_nc( nf90_def_var (fid, 'second',NF90_INT,  (/rec_dimid/), &    ! Second within a minute
         varid_nsecm, chunksizes=chunksizes1D), 'nf90_def_var second')
    IF (glob_att%product_type == PRODUCT_TYPE_SIM) THEN
      status = check_nc( nf90_def_var (fid, 'forecast_seconds',NF90_INT,  (/rec_dimid/), &    ! Forecast time in seconds
           varid_vvzs, chunksizes=chunksizes1D), 'nf90_def_var forecast_seconds')
      status = check_nc( nf90_def_var (fid, 'forecast_minutes',NF90_INT,  (/rec_dimid/), &    ! Forecast time in minutes
           varid_vvzm, chunksizes=chunksizes1D), 'nf90_def_var forecast_minutes')
    END IF

    ! Variables from read_field_obs_dwd
    status = check_nc( nf90_def_var (fid, 'DATE',          NF90_INT,  (/rec_dimid/), &  ! Date
         varid_date,   chunksizes=chunksizes1D), 'nf90_def_var DATE')
    status = check_nc( nf90_def_var (fid, 'TIME',          NF90_INT,  (/rec_dimid/), &  ! Observation Time
         varid_time,   chunksizes=chunksizes1D), 'nf90_def_var TIME')
    status = check_nc( nf90_def_var (fid, 'ppi_azimuth',   NF90_REAL, (/rec_dimid/), &  ! Azimuth mabaz
         varid_mabaz,  chunksizes=chunksizes1D), 'nf90_def_var ppi_azimuth')
    status = check_nc( nf90_def_var (fid, 'ppi_elevation', NF90_REAL, (/rec_dimid/), &  ! Elevation manel
         varid_manel,  chunksizes=chunksizes1D), 'nf90_def_var ppi_elevation')
    status = check_nc( nf90_def_var (fid, 'n_range_bins',  NF90_INT,  (/rec_dimid/), &  ! Number of range bins of individual report
         varid_rbins,  chunksizes=chunksizes1D), 'nf90_def_var n_range_bins')

    status = check_nc( nf90_def_var (fid, 'ray_azimuth',     NF90_REAL, (/az_dimid,rec_dimid/), &    ! Azimuth mabaz0
         varid_mabaz0, chunksizes=chunksizes2D), 'nf90_def_var ray_azimuth')
    status = check_nc( nf90_def_var (fid, 'ray_elevation',   NF90_REAL, (/az_dimid,rec_dimid/), &    ! Elevation manel0
         varid_manel0, chunksizes=chunksizes2D), 'nf90_def_var ray_elevation')


    ! Enable data compression (netcdf4-feature): byte-shuffling and compression with a level between 1 and 9
    !  Note that compression is done chunk by chunk, so that the chunksize ideally has to match the data size for
    !  one "slice" of the unlimited dimension! This is done above in exactly this way.
    status = check_nc( nf90_def_var_deflate(fid, varid_mabaz0, shuffle=1, deflate=1, deflate_level=5), &
         'nf90_def_var_deflate varid_mabaz0')
    status = check_nc( nf90_def_var_deflate(fid, varid_manel0, shuffle=1, deflate=1, deflate_level=5), &
         'nf90_def_var_deflate varid_manel0')

    ! Put variable attributes (like _FVillValue)
    status = check_nc( nf90_put_att(fid,varid_mabaz0 , "_FillValue", missing_val_real), &
         'nf90_put_att _FillValue varid_mabaz0')
    status = check_nc( nf90_put_att(fid,varid_mnipul , "_FillValue", missing_val_int) , &
         'nf90_put_att _FillValue varid_mnipul')

    status = check_nc( nf90_put_att(fid,varid_mabaz0 , "_MissingValue", missing_val_real), &
         'nf90_put_att _MissingValue varid_mabaz0')
    status = check_nc( nf90_put_att(fid,varid_mnipul , "_MissingValue", missing_val_int) , &
         'nf90_put_att _MissingValue varid_mnipul')

    ! Put global attributes (like a data description etc.)
    status = check_nc( nf90_put_att(fid, NF90_GLOBAL , 'Creator', TRIM(glob_att%Creator)), &
         'nf90_put_att Creator')
    status = check_nc( nf90_put_att(fid, NF90_GLOBAL , 'Creation_date', glob_att%Creation_date), &
         'nf90_put_att Creation_date')
    ! Data source (simulation or observation)
    status = check_nc( nf90_put_att(fid, NF90_GLOBAL , 'Data_source', TRIM(glob_att%Data_source)), &
         'nf90_put_att Data_source')
    ! Description of the data (header line from ASCII files):
    status = check_nc( nf90_put_att(fid, NF90_GLOBAL , 'Data_description', TRIM(glob_att%Data_description)), &
         'nf90_put_att Data_description')
    ! Put the list of other flexible global attributes, dependig on their data type:
    DO i=1, glob_att%n_attr
      SELECT CASE (glob_att%att(i)%dtype)
      CASE (CHAR_ATT)
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL , TRIM(glob_att%att(i)%name), TRIM(glob_att%att(i)%char_val)), &
             'nf90_put_att '//TRIM(glob_att%att(i)%name))
      CASE (INT_ATT)
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL , TRIM(glob_att%att(i)%name), glob_att%att(i)%int_val), &
             'nf90_put_att '//TRIM(glob_att%att(i)%name))
      CASE (REAL_ATT)
        status = check_nc( nf90_put_att(fid, NF90_GLOBAL , TRIM(glob_att%att(i)%name), glob_att%att(i)%real_val), &
             'nf90_put_att '//TRIM(glob_att%att(i)%name))
      END SELECT
    END DO
    
    ! End definition mode of netcdf
    status = check_nc( nf90_enddef(fid) , 'nf90_enddef')

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE cdfin_define

  SUBROUTINE cdfin_write_mod_data (fid, rsm, itime, ncdf_shortname, unit_string, field_desc, product_type, &
       &                           actdate, vvz_s, vvz_m, datap, ind_ele_outlist, nele_out)

    ! INCOMING:
    TYPE(radar_meta_type), INTENT(IN) :: rsm                 ! Radar data meta structure
    INTEGER,               INTENT(in) :: itime               ! time index of current time step
    INTEGER,               INTENT(IN) :: fid                 ! file ID to write data
    CHARACTER(len=*),      INTENT(IN) :: ncdf_shortname      ! NetCDF short name of radar data field vr or z usually
    CHARACTER (len=*),     INTENT(in) :: unit_string         ! Unit of data field (vari attribute)
    CHARACTER (len=*),     INTENT(in) :: field_desc          ! Description of data field (vari attribute)
    INTEGER,               INTENT(in) :: product_type        ! either PRODUCT_TYPE_SIM or PRODUCT_TYPE_OBS
    CHARACTER(len=*),      INTENT(IN) :: actdate             ! YYYYMMDDhhmmss for date/time information
    INTEGER,               INTENT(IN) :: vvz_s               ! Forecast time in seconds
    INTEGER,               INTENT(IN) :: vvz_m               ! Forecast time in minutes
    REAL     (KIND=dp),    INTENT(IN) :: datap(:,:,:)        ! Data to write (naz,nra,nel)
    INTEGER,               INTENT(in) :: ind_ele_outlist(:)  ! List of elevation indices for output
    INTEGER,               INTENT(in) :: nele_out            ! length of ind_ele_outlist(:)



    INTEGER  :: status, dummy, ii

    ! Indices
    INTEGER             :: start_report
    ! Dimension  ID
    INTEGER             :: rec_dimid, az_dimid, ra_dimid
    INTEGER             :: len_report    ! Number of reports in formal file
    ! Variables ID
    INTEGER             :: varid_yssosn, varid_mii, varid_niii, varid_mlah, &
                           varid_mloh, varid_mhp, varid_rngbs, &
                           varid_azres, varid_rgoff, varid_azoff, &
                           varid_enyq, varid_hnyq, varid_dprf, varid_mrglen, &
                           varid_mngave, varid_mnipul, varid_mjjj, &
                           varid_mmm, varid_myy, varid_mgg, &
                           varid_ngg, varid_nsecm, varid_vvzs, varid_vvzm, varid_time, varid_date, &
                           varid_mabaz, varid_rbins, varid_manel, &
                           varid_data, varid_manel0, varid_mabaz0

    CHARACTER(len=300) :: varnames_attrib
    INTEGER            :: chunksizes3D(3)
    REAL,    PARAMETER :: missing_val_real =  HUGE(1.0)
    INTEGER, PARAMETER :: missing_val_int  = -HUGE(1)

    INTEGER            :: lastdate_int, lasttime_int, actdate_int, acttime_int
    LOGICAL            :: append_new_datetime

    ! Dummy arrays
    INTEGER             :: ele_array_int(nele_out)
    REAL     (KIND=dp)  :: ele_array_real(nele_out), ele_2Darray_real(rsm%naz,nele_out), &
                           ele_3Darray_real(rsm%nra,rsm%naz,nele_out)
    REAL (KIND=dp), ALLOCATABLE       :: el_precip(:)

    ! Check shape and size of datap:
    IF ( ANY( SHAPE(datap) /= (/ rsm%naz, rsm%nra, rsm%nel /) ) ) THEN
      CALL abort_run (my_radar_id, 17075, &
           'ERROR: problem in cdfin_write_mod_data(): wrong shape of input data field ''datap''! Fix the source code!', &
           'radar_output_methods.f90, cdfin_write_mod_data()')
    END IF

    ! Inquire variables
    status = check_nc( nf90_inq_varid (fid, 'station_name', varid_yssosn)      , 'nf90_inq_varid station_name')
    status = check_nc( nf90_inq_varid (fid, 'country_ID',   varid_mii)         , 'nf90_inq_varid country_ID')
    status = check_nc( nf90_inq_varid (fid, 'station_ID_national',  varid_niii), 'nf90_inq_varid station_ID_national')
    status = check_nc( nf90_inq_varid (fid, 'station_latitude',  varid_mlah)   , 'nf90_inq_varid station_latitude')
    status = check_nc( nf90_inq_varid (fid, 'station_longitude',  varid_mloh)  , 'nf90_inq_varid station_longitude')
    status = check_nc( nf90_inq_varid (fid, 'station_height' ,  varid_mhp)     , 'nf90_inq_varid station_height')
    status = check_nc( nf90_inq_varid (fid, 'range_resolution',varid_rngbs)    , 'nf90_inq_varid range_resolution')
    status = check_nc( nf90_inq_varid (fid, 'azimuthal_resolution',varid_azres), 'nf90_inq_varid azimuthal_resolution')
    status = check_nc( nf90_inq_varid (fid, 'range_start',varid_rgoff)         , 'nf90_inq_varid range_start')
    status = check_nc( nf90_inq_varid (fid, 'azimuth_start',varid_azoff)       , 'nf90_inq_varid azimuth_start')
    status = check_nc( nf90_inq_varid (fid, 'extended_nyquist',varid_enyq)     , 'nf90_inq_varid extended_nyquist')   ! extended Nyquist Velocity
    status = check_nc( nf90_inq_varid (fid, 'high_nyquist',varid_hnyq)         , 'nf90_inq_varid high_nyquist')       ! high Nyquist Velocity
    status = check_nc( nf90_inq_varid (fid, 'dualPRF_ratio',varid_dprf)        , 'nf90_inq_varid dualPRF_ratio')      ! DualPRF ratio
    status = check_nc( nf90_inq_varid (fid, 'range_gate_length',varid_mrglen)  , 'nf90_inq_varid range_gate_length')  ! RANGE-GATE LENGTH
    status = check_nc( nf90_inq_varid (fid, 'n_ranges_averaged',varid_mngave)  , 'nf90_inq_varid n_ranges_averaged')  ! NUMBER OF GATES AVERAGED
    status = check_nc( nf90_inq_varid (fid, 'n_pulses_averaged',varid_mnipul)  , 'nf90_inq_varid n_pulses_averaged')  ! NUMBER OF INTEGRATED PULSES
    status = check_nc( nf90_inq_varid (fid, 'year',  varid_mjjj)               , 'nf90_inq_varid year')    ! Year
    status = check_nc( nf90_inq_varid (fid, 'month',   varid_mmm)              , 'nf90_inq_varid month')   ! Month
    status = check_nc( nf90_inq_varid (fid, 'day',   varid_myy)                , 'nf90_inq_varid day')     ! Day
    status = check_nc( nf90_inq_varid (fid, 'hour',   varid_mgg)               , 'nf90_inq_varid hour')    ! Hour
    status = check_nc( nf90_inq_varid (fid, 'minute',   varid_ngg)             , 'nf90_inq_varid minute')  ! Minute
    status = check_nc( nf90_inq_varid (fid, 'second', varid_nsecm)             , 'nf90_inq_varid second')  ! Second within a minute
    IF (product_type == PRODUCT_TYPE_SIM) THEN
      status = check_nc( nf90_inq_varid (fid, 'forecast_seconds', varid_vvzs)    , 'nf90_inq_varid forecast_seconds')  ! Forecast time in seconds
      status = check_nc( nf90_inq_varid (fid, 'forecast_minutes', varid_vvzm)    , 'nf90_inq_varid forecast_minutes')  ! Forecast time in minutes
    END IF
    status = check_nc( nf90_inq_varid (fid, 'DATE', varid_date)                , 'nf90_inq_varid DATE')
    status = check_nc( nf90_inq_varid (fid, 'TIME'  , varid_time)              , 'nf90_inq_varid TIME')          ! Reference time
    status = check_nc( nf90_inq_varid (fid, 'ppi_elevation' , varid_manel)     , 'nf90_inq_varid ppi_elevation') ! Elevation nominal
    status = check_nc( nf90_inq_varid (fid, 'ppi_azimuth' , varid_mabaz)       , 'nf90_inq_varid ppi_azimuth')   ! Azimuth start (centers of bin)
    status = check_nc( nf90_inq_varid (fid, 'n_range_bins', varid_rbins)       , 'nf90_inq_varid n_range_bins')  ! Bin range

    status = check_nc( nf90_inq_varid (fid, 'ray_azimuth', varid_mabaz0)       , 'nf90_inq_varid ray_azimuth')   ! Azimuth (ends of bins)
    status = check_nc( nf90_inq_varid (fid, 'ray_elevation', varid_manel0)     , 'nf90_inq_varid ray_elevation') ! Elevation

    ! Get the number of existing reports in file:
    status = check_nc( nf90_inquire(fid, unlimitedDimId = rec_dimid)         , 'nf90_inquire rec_dimid')
    status = check_nc( nf90_inquire_dimension(fid, rec_dimid, len=len_report) , 'nf90_inquire_dimension len_report')

    ! Get the other relevant dimension IDs:
    status = check_nc( nf90_inq_dimid(fid, "n_azimuth"  , az_dimid), 'nf90_inq_dimid n_azimuth'  )  ! Azimuhtal dimension
    status = check_nc( nf90_inq_dimid(fid, "n_range"    , ra_dimid), 'nf90_inq_dimid n_range'    )  ! Range dimension

    ! Find out if the new datap is a new datetime or a new moment:
    status = nf90_inq_varid (fid, TRIM(ncdf_shortname), varid_data)  ! Data variable
    
    IF (status /= NF90_NOERR) THEN

      ! For this, first retrieve the global attribute of the list of varnames which are already in the file:
      varnames_attrib(:) = ' '
      status = nf90_get_att(fid, NF90_GLOBAL , 'Varnames', varnames_attrib)
      varnames_attrib = ADJUSTL(TRIM(varnames_attrib)//' '//TRIM(ncdf_shortname))

     ! Now, re-define cdfin file to add a new variable with name ncdf_shortname:
      status = check_nc ( nf90_redef( fid ), 'nf90_redef fid')
          
      ! Add 3D data variable:
      chunksizes3D(:) = (/rsm%nra, rsm%naz, nele_out/)
      status = check_nc( nf90_def_var (fid, TRIM(ncdf_shortname) ,NF90_REAL, (/ra_dimid,az_dimid,rec_dimid/), &
           varid_data, chunksizes=chunksizes3D), 'nf90_def_var ')

      ! Enable compression:
      status = check_nc( nf90_def_var_deflate(fid, varid_data,   shuffle=1, deflate=1, deflate_level=5), &
           'nf90_def_var_deflate varid_data')

      ! Put variable attributes (like _FVillValue)
      status = check_nc( nf90_put_att(fid,varid_data   , 'Unit', TRIM(unit_string)), &
           'nf90_put_att _FillValue varid_data unit_string')
      status = check_nc( nf90_put_att(fid,varid_data   , 'Description', TRIM(field_desc)), &
           'nf90_put_att _FillValue varid_data field_desc')
      status = check_nc( nf90_put_att(fid,varid_data   , "_FillValue", missing_val_real), &
           'nf90_put_att _FillValue varid_data _FillValue')
      status = check_nc( nf90_put_att(fid,varid_data   , "_MissingValue", missing_val_real), &
           'nf90_put_att _MissingValue varid_data _MissingValue')

      status = check_nc( nf90_put_att(fid, NF90_GLOBAL , 'Varnames', TRIM(varnames_attrib)), &
           'nf90_put_att Varnames append')

      ! End definition mode of netcdf
      status = check_nc( nf90_enddef(fid) , 'nf90_enddef fid')

      ! The actual data will be written to file below!
      
    END IF

    ! Check if we have a new datetime and if we have to write all the additional attributes and metadata for this datetime:
    IF ( len_report < 1 ) THEN

      ! This is an empty file with no data yet written to. A new datetime with all necessary infos has to be added:
      append_new_datetime = .TRUE.
      ! Where to start writing:
      start_report = len_report  + 1

    ELSE
      
      ! This is a file where data of any moment(s) already have been stored.
      ! If the actual datetime is the datetime of the last record in the file,
      !  we don't have to add all the additional attributes and metadata a second time:
      status = check_nc( nf90_get_var (fid, varid_date, lastdate_int, (/len_report/)), 'nf90_get_var lastdate' )
      status = check_nc( nf90_get_var (fid, varid_time, lasttime_int, (/len_report/)), 'nf90_get_var lasttime' )

      READ (actdate(1:14),'(I8,I6)') actdate_int, acttime_int  !YearMonthDay, HourMinuteSecond

      IF (actdate_int /= lastdate_int .OR. acttime_int /= lasttime_int) THEN

        append_new_datetime = .TRUE.
        start_report = len_report  + 1
      
      ELSE
      
        append_new_datetime = .FALSE.
        start_report = len_report - nele_out + 1
      
      END IF

    END IF

    IF ( append_new_datetime ) THEN

      ! Put station name into the file:
      DO ii = 1, nele_out
        status =  check_nc( nf90_put_var(fid,varid_yssosn, rsm%station_name, &
             start=(/1,start_report+ii-1/), count=(/LEN(rsm%station_name),1/)), 'nf90_put_var varid_yssosn')
      END DO

      ! Set constant variables
      ele_array_int(:) = FLOOR(rsm%station_id/1000.0)
      status =  check_nc( nf90_put_var(fid,varid_mii,ele_array_int, start=(/start_report/) ), 'nf90_put_var varid_mii')

      ele_array_int(:) = rsm%station_id - ele_array_int(:)*1000
      status =  check_nc( nf90_put_var(fid,varid_niii,ele_array_int,start=(/start_report/) ), 'nf90_put_var varid_niii')

      ele_array_real(:) = rsm%lat
      status =  check_nc( nf90_put_var(fid,varid_mlah,ele_array_real,start=(/start_report/) ), 'nf90_put_var varid_mlah')

      ele_array_real(:) = rsm%lon
      status =  check_nc( nf90_put_var(fid,varid_mloh,ele_array_real,start=(/start_report/) ), 'nf90_put_var varid_mloh')

      ele_array_real(:) =  rsm%alt_msl
      status =  check_nc( nf90_put_var(fid,varid_mhp,ele_array_real,start=(/start_report/) ), 'nf90_put_var varid_mhp')

      ele_array_real(:) = rsm%ra_inc
      status =  check_nc( nf90_put_var(fid,varid_rngbs,ele_array_real,start=(/start_report/)) , 'nf90_put_var varid_rngbs')

      ele_array_real(:) = rsm%az_inc
      status =  check_nc( nf90_put_var(fid,varid_azres,ele_array_real,start=(/start_report/)), 'nf90_put_var varid_azres')

      ele_array_real(:) = rsm%ra_inc  ! ra_start = ra_inc!
      status =  check_nc( nf90_put_var(fid,varid_rgoff,ele_array_real,start=(/start_report/)) , 'nf90_put_var varid_rgoff')

      ele_array_real(:) = rsm%az_start
      status =  check_nc( nf90_put_var(fid,varid_azoff,ele_array_real,start=(/start_report/)), 'nf90_put_var varid_azoff')

      ele_array_real(:) = rsm%ext_nyq(ind_ele_outlist(1:nele_out),itime)
      status =  check_nc( nf90_put_var(fid,varid_enyq,ele_array_real,start=(/start_report/)) , 'nf90_put_var varid_enyq')

!!$      ele_array_real(:) = rsm%high_nyq(ind_ele_outlist(1:nele_out))
      ele_array_real(:) = rsm%ext_nyq(ind_ele_outlist(1:nele_out),itime)
      status =  check_nc( nf90_put_var(fid,varid_hnyq,ele_array_real,start=(/start_report/)), 'nf90_put_var varid_hnyq')

      ele_array_real(:) = rsm%dualprf_ratio(ind_ele_outlist(1:nele_out))
      status =  check_nc( nf90_put_var(fid,varid_dprf,ele_array_real,start=(/start_report/)), 'nf90_put_var varid_dprf')

      ele_array_real(:) = rsm%rngate_len
      status =  check_nc( nf90_put_var(fid,varid_mrglen,ele_array_real,start=(/start_report/)), 'nf90_put_var varid_mrglen')

      ele_array_int(:) = rsm%num_gates
      status =  check_nc( nf90_put_var(fid,varid_mngave,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_mngave')

      ele_array_int(:) = rsm%num_pulses
      status =  check_nc( nf90_put_var(fid,varid_mnipul,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_mnipul')

      READ (actdate(1:4),'(I4)') dummy   !Year
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_mjjj,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_mjjj')

      READ (actdate(5:6),'(I2)') dummy   !Month
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_mmm,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_mmm')

      READ (actdate(7:8),'(I2)') dummy    !Day
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_myy,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_myy')


      READ (actdate(9:10),'(I2)') dummy    !Hour
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_mgg,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_mgg')

      READ (actdate(11:12),'(I2)') dummy    !Minute
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_ngg,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_ngg')

      READ (actdate(13:14),'(I2)') dummy     !Second within a minute
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_nsecm,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_nsecm')

      IF (product_type == PRODUCT_TYPE_SIM) THEN
        ele_array_int(:) = vvz_s
        status = check_nc( nf90_put_var(fid,varid_vvzs,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_vvzs')

        ele_array_int(:) = vvz_m
        status = check_nc( nf90_put_var(fid,varid_vvzm,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_vvzm')
      END IF
    
      ! Variables from read_field_obs_dwd

      ! Get time and date from actdate string
      READ (actdate(1:8),'(I8)') dummy
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_date,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_date')
      READ (actdate(9:14),'(I6)') dummy
      ele_array_int(:) = dummy
      status = check_nc( nf90_put_var(fid,varid_time,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_time')

      ele_array_real(:) = rsm%el_arr(ind_ele_outlist(1:nele_out))
      status = check_nc( nf90_put_var(fid,varid_manel,ele_array_real,start=(/start_report/)), &
           'nf90_put_var varid_manel')

      ! the azimut start values represent the center of the first azimuth bin of the data:
      ele_array_real(:) = rsm%az_start
      status = check_nc( nf90_put_var(fid,varid_mabaz,ele_array_real,start=(/start_report/)), 'nf90_put_var varid_mabaz')

      ele_array_int(:) = rsm%nra
      status = check_nc( nf90_put_var(fid,varid_rbins,ele_array_int,start=(/start_report/)), 'nf90_put_var varid_rbins')

      ! Elevations for each ray, taking into account DWD precipitation scans:
      IF (TRIM(rsm%scanname) == 'PRECIP') THEN
        CALL get_elarr_precipscan (rsm, el_precip)
        ele_2Darray_real(:,1) = el_precip
      ELSE
        DO ii = 1, rsm%naz
          ele_2Darray_real(ii,:) = rsm%el_arr(ind_ele_outlist(1:nele_out))
        END DO
      END IF
      status =check_nc(  nf90_put_var(fid,varid_manel0,ele_2Darray_real,start=(/1,start_report/),count=(/rsm%naz,nele_out/)), &
           'nf90_put_var varid_manel0')

      ! Azimuths for each ray: values should represent the end of the averaging interval, not the center, therefore the +0.5*az_inc:
      DO ii = 1,rsm%naz
        ele_2Darray_real(ii,:) = rsm%az_start + (ii-1+0.5)* rsm%az_inc   ! end of azimuth interval
!        ele_2Darray_real(ii,:) = rsm%az_start + (ii-1)* rsm%az_inc       ! center of azimuth interval
      END DO
      status = check_nc( nf90_put_var(fid,varid_mabaz0,ele_2Darray_real,start=(/1,start_report/),count=(/rsm%naz,nele_out/)), &
           'nf90_put_var varid_mabaz0')

    END IF

    !=====================================================================================
    ! Write data of this moment to the variable named ncdf_shortname:
    
    ! Reshape data to (nra, naz, nel):
    DO ii = 1, nele_out
      ele_3Darray_real(:,:,ii) = TRANSPOSE(datap(:,:,ind_ele_outlist(ii)))
    END DO
    ! Put data into the file
    status = check_nc( nf90_put_var(fid,varid_data,ele_3Darray_real,start=(/1,1,start_report/), &
         count=(/rsm%nra,rsm%naz,nele_out/)), 'nf90_put_var varid_data')
    
    !=====================================================================================
    ! Clean up:
    IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)
    
    
  END SUBROUTINE cdfin_write_mod_data

#endif

  !==================================================================================================
  !==================================================================================================

  SUBROUTINE fill_cdfin_global_attributes (rsm, grib_input, grib_locinfo, config_string, glob_att)

    !------------------------------------------------------------------------------
    !
    ! Description: prepares the writing of a radar volume scan, either as cdfin or grib2.
    !              Some derived types are filled with informations on date/time
    !              and attributes of the model run (model type, run IDs, ensemble codings, etc.)
    !
    ! Call this routine before you call
    !              write_cdfin_from_volscan()
    !              write_grib_from_volscan()
    !
    !------------------------------------------------------------------------------

    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! INCOMING:
    TYPE(radar_meta_type)  , INTENT(IN)    :: rsm
    TYPE(t_grib)           , INTENT(in)    :: grib_input
    TYPE(t_grib2_modelspec), INTENT(in)    :: grib_locinfo
    CHARACTER(len=*)       , INTENT(in)    :: config_string
    TYPE(t_cdfin_globalatt), INTENT(inout) :: glob_att
    
    ! LOCAL:
    ! -none-
    
    
    ! First define the fixed global attributes:
    glob_att%product_type = grib_input%product_type
    
    IF (glob_att%product_type == PRODUCT_TYPE_OBS) THEN
      glob_att%Data_source = 'Observation'
    ELSE
      glob_att%Data_source = 'Model simulation'
    END IF

    glob_att%Creator           = TRIM(cdfin_creator_model)
    glob_att%Data_description  = TRIM(config_string)
    glob_att%Creation_Date     = grib_input%creationDate

    ! Now fill the variable global attributes:
    IF (.NOT.ALLOCATED(glob_att%att)) ALLOCATE(glob_att%att(NATTR_MAX))
   
    glob_att%n_attr = 0

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'alt_msl_true'
    glob_att%att(glob_att%n_attr)%dtype   = REAL_ATT
    glob_att%att(glob_att%n_attr)%real_val = rsm%alt_msl_true

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'datetime_ini_model'
    glob_att%att(glob_att%n_attr)%dtype   = CHAR_ATT
    glob_att%att(glob_att%n_attr)%char_val = grib_input%dateTime_ini
    
    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'tablesVersion'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%tablesVersion

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'generatingCenter'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%generatingCenter

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'generatingSubcenter'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%generatingSubcenter

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'significanceOfReferenceTime'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%significanceOfReferenceTime

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'Timestamp_meaning'
    glob_att%att(glob_att%n_attr)%dtype   = CHAR_ATT
    IF (glob_att%product_type == PRODUCT_TYPE_OBS) THEN
      glob_att%att(glob_att%n_attr)%char_val = 'Observation time'
    ELSE
      glob_att%att(glob_att%n_attr)%char_val = 'Forecast valid time'
    END IF

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'productDefinitionTemplateNumber'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%productDefinitionTemplateNumber

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'productionStatusOfProcessedData'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%productionStatusOfProcessedData

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'typeOfProcessedData'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%typeOfProcessedData

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'typeOfGeneratingProcess'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%typeOfGeneratingProcess

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'backgroundProcess'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%backgroundProcess

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'generatingProcessIdentifier'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%generatingProcessIdentifier

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'localDefinitionNumber'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%localDefinitionNumber

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'localNumberOfExperiment'
    glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
    glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%localNumberOfExperiment

    CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
    glob_att%att(glob_att%n_attr)%name    = 'model_run_type'
    glob_att%att(glob_att%n_attr)%dtype   = CHAR_ATT
    IF (grib_locinfo%numberOfForecastsInEnsemble > 0) THEN
      glob_att%att(glob_att%n_attr)%char_val = 'Ensemble'
    ELSE
      glob_att%att(glob_att%n_attr)%char_val = 'Deterministic'
    END IF

    IF (grib_locinfo%numberOfForecastsInEnsemble > 0) THEN

      ! This is an ensemble run:

      CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
      glob_att%att(glob_att%n_attr)%name    = 'typeOfEnsembleForecast'
      glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
      glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%typeOfEnsembleForecast

      CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
      glob_att%att(glob_att%n_attr)%name    = 'localTypeOfEnsembleForecast'
      glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
      glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%localTypeOfEnsembleForecast

      CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
      glob_att%att(glob_att%n_attr)%name    = 'numberOfForecastsInEnsemble'
      glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
      glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%numberOfForecastsInEnsemble

      CALL incr_and_check_counter (counter=glob_att%n_attr, incr=1)
      glob_att%att(glob_att%n_attr)%name    = 'perturbationNumber'
      glob_att%att(glob_att%n_attr)%dtype   = INT_ATT
      glob_att%att(glob_att%n_attr)%int_val = grib_locinfo%perturbationNumber
      
    END IF

  CONTAINS

    SUBROUTINE incr_and_check_counter (counter, incr)
      INTEGER, INTENT(inout) :: counter
      INTEGER, INTENT(in)    :: incr
      CHARACTER(len=cmaxlen) :: errormsg
      
      IF (counter+incr > NATTR_MAX) THEN
        errormsg(:) = ' '
        errormsg    = 'ERROR fill_cdfin_global_attributes(): too much attributes, increase NATTR_MAX in radar_data_io.f90!'
        CALL abort_run(my_radar_id, 321, TRIM(errormsg), 'fill_cdfin_global_attributes()')
      ELSE
        counter = counter + incr
      END IF
    END SUBROUTINE incr_and_check_counter
    
  END SUBROUTINE fill_cdfin_global_attributes

  !==================================================================================================
  !==================================================================================================

  SUBROUTINE prepare_write_volscan (rsm, itime, varname, grib_input, grib_locinfo, error, errmsg)

    !------------------------------------------------------------------------------
    !
    ! Description: prepares the writing of a radar volume scan, either as cdfin or grib2.
    !              Some derived types are filled with informations on date/time
    !              and attributes of the model run (model type, run IDs, ensemble codings, etc.)
    !
    ! Call this routine before you call
    !              write_cdfin_from_volscan()
    !              write_grib_from_volscan()
    !
    ! Note: the actual time returned in the attributes is the nominal obs_time of the
    !       scan as given in the list of obs_times in the input rsm type.
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! .. Input/Output variables:
    TYPE(radar_meta_type)  , INTENT(IN)    :: rsm
    INTEGER                , INTENT(in)    :: itime
    CHARACTER(len=*)       , INTENT(in)    :: varname
    TYPE(t_grib)           , INTENT(out)   :: grib_input
    TYPE(t_grib2_modelspec), INTENT(out)   :: grib_locinfo
    INTEGER                , INTENT(out)   :: error
    CHARACTER(len=*)       , intent(out)   :: errmsg

   ! Local:
    CHARACTER(len=*), PARAMETER :: yzroutine = 'prepare_write_volscan'
    INTEGER                     :: ierr, ierr_diff(2,6)

    ! Date/time array
    CHARACTER(len=8 )           :: currdate
    CHARACTER(len=10)           :: currtime
    CHARACTER(len=LEN_DATETIME) :: utcdate
    INTEGER                     :: year, month, day, hour, minute, second, isecdiff, datevalues(8)

    error = 0
    errmsg(:) = ' '

    IF (INDEX(TRIM(varname),'obs') > 0) THEN
      grib_input%product_type           = PRODUCT_TYPE_OBS
    ELSE
      grib_input%product_type           = PRODUCT_TYPE_SIM
    END IF

    !--------------------------------------------------------
    ! Sort out date and time informations:
    
    grib_input%dateTime_ini = ydate_ini_mod
    ! Time of the nearest obs time to the actual model time:
    grib_input%dateTime_act = new_datetime( grib_input%dateTime_ini, rsm%obs_times(itime) )
    IF (grib_input%product_type == PRODUCT_TYPE_OBS) THEN
      ! The time information is the actual time:
      READ (grib_input%dateTime_act, '(i4,5i2,i2,i2)') year, month, day, hour, minute, second
      grib_input%forecastseconds = 0
      grib_input%forecastminutes = 0
    ELSE
      ! The forecast lead time is given in minutes and seconds since forecast init time:
      CALL diff_seconds ( grib_input%dateTime_ini, grib_input%dateTime_act, isecdiff, ierr_diff )
      IF (isecdiff < 0) THEN
        ! HANDWAVING FIX for grib2 not beeing able to have negative forecast times:
        !  Set the model start time to the actual forecast time:
        grib_input%dateTime_ini = grib_input%dateTime_act
        ! Set the forecast time to 0 minutes:
        grib_input%forecastseconds = 0
        grib_input%forecastminutes = 0
      ELSE
        grib_input%forecastseconds =  ABS( isecdiff )
        grib_input%forecastminutes = ABS( NINT ( REAL(isecdiff,dp) / 60.0_dp ) )
      END IF
      ! The time is the forecast init time:
      READ (grib_input%dateTime_ini, '(i4,5i2,i2,i2)') year, month, day, hour, minute, second
    END IF

    grib_input%year   = year
    grib_input%month  = month
    grib_input%day    = day
    grib_input%hour   = hour
    grib_input%minute = minute
    grib_input%second = second

    CALL DATE_AND_TIME(date=currdate, time=currtime, values=datevalues(:))  ! this in is in local time, but we need UTC
    ! Convert to UTC:
    ! datevalues(4) contains the difference to UTC in minutes. Use get_utc_date() to compute UTC time:
    utcdate = new_datetime(currdate(1:8)//currtime(1:6), -datevalues(4)*60.0_dp)
    grib_input%creationDate = TRIM(utcdate)

    ! Get model specific informations (center, ensemble, process IDs, ...)
    CALL grib2_add_modelspec_info(idom, grib_locinfo, ierr)

    error = ierr
    
  END SUBROUTINE prepare_write_volscan
  
  !==================================================================================================
  !==================================================================================================
  
  SUBROUTINE write_grib_from_volscan(rsm, itime, grib_input_in, grib_locinfo_in, &
                                     grib2_filename, force_overwrite, grib2_shortname, grib2_packingtype, &
                                     config_string, dat, &
                                     ind_ele_outlist, nele_out, error, errmsg)

    !------------------------------------------------------------------------------
    !
    ! Description: creates a grib file of a radar volume scan.
    !
    ! Call prepare_write_volscan() before!
    !
    !------------------------------------------------------------------------------

    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! INCOMING:
    TYPE(radar_meta_type)  , INTENT(IN)    :: rsm
    INTEGER                , INTENT(in)    :: itime
    TYPE(t_grib)           , INTENT(IN)    :: grib_input_in
    TYPE(t_grib2_modelspec), INTENT(IN)    :: grib_locinfo_in
    REAL(KIND=dp)          , INTENT(IN)    :: dat(:,:,:)             ! Radar data 3D (naz,nra,nel)
    CHARACTER (len=*)      , INTENT(in)    :: grib2_filename         ! grib2 filename
    LOGICAL                , INTENT(in)    :: force_overwrite        ! if .true., forcefully overwrite any pre-existing output file of name "grib2_filename"
    CHARACTER (len=*)      , INTENT(in)    :: grib2_shortname        ! short name of radar data field
    CHARACTER (len=*)      , INTENT(in)    :: grib2_packingtype      ! 'grid_simple', 'grid_ccsds', 'png'
    CHARACTER (len=*)      , INTENT(in)    :: config_string          ! Description of data set
    INTEGER                , INTENT(in)    :: ind_ele_outlist(:)
    INTEGER                , INTENT(in)    :: nele_out               ! length of ind_ele_outlist
    INTEGER                , INTENT(out)   :: error
    CHARACTER(len=*)       , intent(out)   :: errmsg

    ! Local:
    CHARACTER(len=*), PARAMETER :: yzroutine = 'write_grib_from_volscan'
    INTEGER                     :: istat, jrecord, ind_ele, ierr
    CHARACTER(len=cmaxlen)      :: errormsg, stationstr
    LOGICAL                     :: lclobber

    TYPE(t_grib)                :: grib_input
    TYPE(t_grib2_modelspec)     :: grib_locinfo

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    error = 0
    errmsg(:) = ' '

    !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    ! Copy input types to local variables for safety:
    grib_input   = grib_input_in
    grib_locinfo = grib_locinfo_in
    
    !--------------------------------------------------------------------------------
    ! Put some few elevation-invariant meta data to the grib_input-type.
    ! The rest will be filled below by fill_grib_input_volscan() for each elevation.
    !--------------------------------------------------------------------------------

    grib_input%varname(:) = ' '
    grib_input%varname    = TRIM(grib2_shortname)

    SELECT CASE (TRIM(grib2_shortname))
    CASE ('zrobs')
      
      grib_input%product_type           = PRODUCT_TYPE_OBS

      grib_input%sec0%discipline        = 0   ! 0.0.table:      0  Meteorological products
      grib_input%sec4%parameterCategory = 15
      grib_input%sec4%parameterNumber   = 200 ! DBZSCAN_OBS

      ! Bitmap: if bitmap, and if yes, which value to use for the bitmap?
      !   BITMAP_PRI_NOBITMAP           = no bitmap at all
      !   BITMAP_PRI_MISSVAL            = miss_value  (-999.99)
      !   BITMAP_PRI_UNDEVAL            = zero_value  ( -99.99)
      !   BITMAP_PRI_MISSVAL_UNDEVAL    = miss_value if miss_values present, otherwise zero_value
      grib_input%bitmap_priority = BITMAP_PRI_UNDEVAL   ! good for reflectivity
      
    CASE ('zrsim')
      
      grib_input%product_type           = PRODUCT_TYPE_SIM
      
      grib_input%sec0%discipline        = 0   ! 0.0.table:      0  Meteorological products
      grib_input%sec4%parameterCategory = 16
      grib_input%sec4%parameterNumber   = 4   ! DBZSCAN_SIM

      ! Bitmap: if bitmap, and if yes, which value to use for the bitmap?
      !   BITMAP_PRI_NOBITMAP           = no bitmap at all
      !   BITMAP_PRI_MISSVAL            = miss_value  (-999.99)
      !   BITMAP_PRI_UNDEVAL            = zero_value  ( -99.99)
      !   BITMAP_PRI_MISSVAL_UNDEVAL    = miss_value if miss_values present, otherwise zero_value
      grib_input%bitmap_priority = BITMAP_PRI_UNDEVAL   ! good for reflectivity
      
    CASE default
      
      errormsg(:) = ' '
      errormsg    = 'ERROR: grib2 output of '//TRIM(grib2_shortname)//' not implemented!'
      CALL abort_run(my_radar_id, 123, TRIM(errormsg), TRIM(yzroutine)//', fill_grib_input_volscan()')

    END SELECT

    IF (TRIM(rsm%scanname) == 'PRECIP') THEN
      grib_input%localOperatingMode_in = 2
    ELSE
      grib_input%localOperatingMode_in = 1
    END IF

    SELECT CASE (TRIM(grib2_packingtype))
    CASE ('grid_simple', 'grid_ccsds', 'png')
      grib_input%packingtype(:) = ' '
      grib_input%packingtype    = TRIM(grib2_packingtype)
    CASE default
      errormsg(:) = ' '
      errormsg    = 'ERROR: grib2_packingtype '//TRIM(grib2_packingtype)//' not implemented!'
      CALL abort_run(my_radar_id, 124, TRIM(errormsg), TRIM(yzroutine)//', fill_grib_input_volscan()')      
    END SELECT
    
    !--------------------------------------------------------------------------------
    ! Write GRIB-File
    !--------------------------------------------------------------------------------

    lclobber = force_overwrite
    
    DO jrecord = 1, nele_out

      ind_ele = ind_ele_outlist(jrecord)

      stationstr(:) = ' '
      WRITE (stationstr, '("Ele ",i0," for ",a," of station ",i0)') ind_ele, TRIM(grib2_shortname), rsm%station_id

      istat = 0
      errormsg(:) = ' '
      CALL fill_grib_input_volscan ( rsm     = rsm,         &
           &                         itime   = itime,       &
           &                         emvo_config_str = TRIM(config_string), &
           &                         dat     = dat,         &
           &                         ind_ele = ind_ele,     &
           &                         grib_locinfo = grib_locinfo, &
           &                         grib    = grib_input,  &
           &                         lverbose= ldebug_radsim, &
           &                         error   = istat,       &
           &                         errmsg  = errormsg     )

      IF (istat /= 0) THEN
        CALL abort_run(my_radar_id, istat, TRIM(errormsg)//' '//TRIM(stationstr), TRIM(yzroutine)//', fill_grib_input_volscan()')
      END IF

      ! NOTE 1: write_grib_to_file(grib_outfilename, ...) will create a new file if file "grib_outfilename"
      !         does not exist yet, or it will append to it if it exists.
      !
      ! NOTE 2: fill_grib_input_volscan() has allocated the following arrays in type(grib):
      ! ---------------------------------------------------------------------------------
      ! - grib_input%sec2%localIndexOfScaledMetadata
      ! - grib_input%sec2%localScaledMetadataTableNumber
      ! - grib_input%sec2%localScaledMetadataId
      ! - grib_input%sec2%localScaleFactorOfMetadata
      ! - grib_input%sec2%localScaledValueOfMetadata
      ! - grib_input%sec3%startingAzimuth
      ! - grib_input%sec3%azimuthalWidth
      ! - grib_input%sec7%values
      !
      ! write_grib_to_file() will deallocate these after having transferred the data to grib:

      istat = 0
      errormsg(:) = ' '
      CALL write_grib_to_file ( grib_outfilename = TRIM(grib2_filename), &
           &                    force_overwrite  = lclobber, &
           &                    grib             = grib_input,  &
           &                    lverbose         = ldebug_radsim, &
           &                    error            = istat,       &
           &                    errmsg           = errormsg     )

      IF (istat /= 0) THEN
        IF (labort_if_problems_gribout) THEN
          CALL abort_run(my_radar_id, istat, TRIM(errormsg)//' '//TRIM(stationstr), &
               TRIM(yzroutine)//', write_grib_to_file()')
        ELSE
          WRITE (*,'(a)') 'WARNING emvorado::'//TRIM(yzroutine)// &
               '::write_grib_to_file: writing volscan grib: '//TRIM(errormsg)//' '//TRIM(stationstr)
        END IF
      END IF

      ! After first record switch lclobber to false so that subsequent elevations go in the same file:
      lclobber = .FALSE.
      
    END DO

    !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE write_grib_from_volscan

  !===============================================================================

  SUBROUTINE fill_grib_input_volscan ( rsm,   &
       &                               itime, &
       &                               emvo_config_str, &
       &                               dat, &
       &                               ind_ele,   &
       &                               grib_locinfo, &
       &                               grib,  &
       &                               lverbose, &
       &                               error, &
       &                               errmsg )

    ! in/out
    TYPE(radar_meta_type), INTENT(IN)  :: rsm
    INTEGER,               INTENT(in)  :: itime
    CHARACTER (len=*),     INTENT(in)  :: emvo_config_str  ! Info-string on global config of EMVORADO, so far not used for grib output
    REAL(KIND=dp),         INTENT(IN)  :: dat(:,:,:)       ! Radar data 3D (naz,nra,nel)
    INTEGER,               INTENT(in)  :: ind_ele
    TYPE(t_grib2_modelspec)            :: grib_locinfo     ! INPUT: general grib metadata for the model run
    TYPE(t_grib),          INTENT(inout) :: grib             ! (mainly) OUTPUT: all resulting grib metadata and values
    LOGICAL,               INTENT(in)  :: lverbose
    INTEGER,               INTENT(OUT) :: error
    CHARACTER(len=*),      INTENT(OUT) :: errmsg

    ! local
    CHARACTER(len=*), PARAMETER :: yzroutine = 'fill_grib_input_volscan'

!!$    CHARACTER(LEN=MAX_LEN_TAB_ENTRY), ALLOCATABLE :: scaled_metadata_list(:)
!!$    REAL(dp), ALLOCATABLE                         :: value_list(:)

    INTEGER                     :: n_azimuth
    INTEGER                     :: n_range
    INTEGER                     :: n_elevation

    REAL(kind=dp), ALLOCATABLE  :: values(:)

    INTEGER                     :: jvalue, nvalue
    INTEGER                     :: jazimuth, jrange

    LOGICAL                     :: lopened, lexist

    REAL(dp)                    :: lower_bound, upper_bound, bitmap_slack
    INTEGER                     :: nbitmap, jhrchy, nhrchy

    !------------------------------------------------------

    error     = 0
    errmsg(:) = ' '

    IF (lverbose) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !----------------------------------------------------------------
    !                         Preparations
    !----------------------------------------------------------------

    ! Read contents of local grib table 2.120.table:
    CALL get_loc_sec_tab(tab2=grib%tab2, verbose=lverbose)

    ! Name of grib2 sample file:
    grib%sample_file(:) = ' '
    grib%sample_file    = 'DWD_GRIB2'

    ! Dimensions of data set:
    n_azimuth   = SIZE(dat,1)
    n_range     = SIZE(dat,2)
    n_elevation = SIZE(dat,3)

    ! Local creation date:
    READ(grib%creationDate,'(i4,5i2)') grib%sec2%localCreationDateYear, grib%sec2%localCreationDateMonth, &
                                       grib%sec2%localCreationDateDay, grib%sec2%localCreationDateHour, &
                                       grib%sec2%localCreationDateMinute, grib%sec2%localCreationDateSecond


    ! Informations for bitmap usage:
    !------------------------------

    nvalue = n_range * n_azimuth
    ALLOCATE( grib%sec7%values(nvalue) )

    jvalue = 0
    DO jazimuth=1, n_azimuth
      DO jrange=1, n_range
        jvalue = jvalue + 1
        grib%sec7%values(jvalue) = dat(jazimuth,jrange,ind_ele)
      ENDDO  !jrange
    ENDDO  !jazimuth

    ! Vorbelegung
    grib%sec6%bitmapPresent = GRIB_FALSE

    SELECT CASE(grib%bitmap_priority)
    CASE(BITMAP_PRI_NOBITMAP)

      grib%sec6%bitmapValue = 0.0_dp        ! dummy

      grib%sec2%localTypeOfDataMaskedByBitmap = GRIB_MISSING  ! 2.120.table: 255 Missing

      nhrchy = 0

    CASE(BITMAP_PRI_MISSVAL)

      grib%sec6%bitmapValue = miss_value   ! missing  -999.99

      grib%sec2%localTypeOfDataMaskedByBitmap = 1  ! 2.120.table: 1 Missing value

      nhrchy = 1

    CASE(BITMAP_PRI_UNDEVAL)

      grib%sec6%bitmapValue =zero_value    ! undetected refl  -99.99

      grib%sec2%localTypeOfDataMaskedByBitmap = 0  ! 2.120.table: 0 Undetected value

      nhrchy = 1

    CASE(BITMAP_PRI_MISSVAL_UNDEVAL)

      grib%sec6%bitmapValue = miss_value

      grib%sec2%localTypeOfDataMaskedByBitmap = 1  ! 2.120.table: 1 Missing value

      nhrchy = 2

    END SELECT

    bitmap_slack = 0.01_dp
    LOOP_HRCHY : DO jhrchy=1, nhrchy
      lower_bound = grib%sec6%bitmapValue - bitmap_slack
      upper_bound = grib%sec6%bitmapValue + bitmap_slack
      nbitmap     = 0
      LOOP_DATA : DO jvalue=1, nvalue
        IF ( grib%sec7%values(jvalue) >= lower_bound .AND. &
             &  grib%sec7%values(jvalue) <= upper_bound       ) THEN
          ! Der Datenpunkt muss auf bitmap_value gesetzt werden,
          ! sonst wird EcCodes den Datenpunkt nicht ausblenden,
          ! wenn er nicht genau den Wert von bitmap_value hat
          grib%sec7%values(jvalue) = grib%sec6%bitmapValue
          nbitmap = nbitmap + 1
        ENDIF
      ENDDO LOOP_DATA

      IF(grib%bitmap_priority == BITMAP_PRI_MISSVAL_UNDEVAL) THEN
        IF(nbitmap > 0) THEN
          EXIT LOOP_HRCHY
        ELSE
          grib%sec6%bitmapValue = zero_value
          grib%sec2%localTypeOfDataMaskedByBitmap = 0  ! 2.120.table: 0 Undetected value
        ENDIF
      ENDIF

    ENDDO LOOP_HRCHY

    IF(nbitmap > 0 .AND. nhrchy > 0) THEN
      ! Switch on gitmap
      grib%sec6%bitmapPresent = GRIB_TRUE
      ! Coding of the bitmap'ed value
      CALL get_scaling( val              = grib%sec6%bitmapValue,                   &
           &            scale_factor     = grib%sec2%localScaleFactorOfBitmapValue, &
           &            scaled_val       = grib%sec2%localScaledValueOfBitmapValue, &
           &            max_scale_factor = 8                                        )
    ELSE
      grib%sec2%localScaleFactorOfBitmapValue = 0
      grib%sec2%localScaledValueOfBitmapValue = 0
    ENDIF

    IF(lverbose) THEN
      WRITE(*,'(a,a,i4,a,i7,a,i7)') 'INFO '//TRIM(yzroutine)//': ', &
           'GRIB-Message for elevation: ', ind_ele, &
           ', number of data points: ', nvalue, &
           ', number of bitmap points: ', nbitmap
    ENDIF

    !----------------------------------------------------------------
    !                     INDICATOR SECTION 0
    !----------------------------------------------------------------

    grib%sec0%discipline    = 0               ! 0.0.table:      0  Meteorological products
    grib%sec0%editionNumber = GRIB_EDITION_2


    !----------------------------------------------------------------
    !                   IDENTIFICATION SECTION 1
    !----------------------------------------------------------------

    grib%sec1%tablesVersion = grib_locinfo%tablesVersion

    grib%sec1%productionStatusOfProcessedData = grib_locinfo%productionStatusOfProcessedData

    SELECT CASE(grib%product_type)
    CASE(PRODUCT_TYPE_OBS)
      grib%sec1%typeOfProcessedData = 7           ! 1.4.table: 7  Processed radar observations
    CASE(PRODUCT_TYPE_SIM)
      grib%sec1%typeOfProcessedData = grib_locinfo%typeOfProcessedData
    END SELECT

    grib%sec1%centre = grib_locinfo%generatingCenter

    grib%sec1%subCentre = grib_locinfo%generatingsubCenter

    SELECT CASE(grib%product_type)
    CASE(PRODUCT_TYPE_OBS)
      grib%sec1%significanceOfReferenceTime = 3   ! 1.1.table: 3  Observation time
    CASE(PRODUCT_TYPE_SIM)
      grib%sec1%significanceOfReferenceTime = 1   ! 1.1.table: 1  Start of forecast
    END SELECT
    grib%sec1%year   = grib%year
    grib%sec1%month  = grib%month
    grib%sec1%day    = grib%day
    grib%sec1%hour   = grib%hour
    grib%sec1%minute = grib%minute
    grib%sec1%second = grib%second

    !----------------------------------------------------------------
    !                     LOCAL USE SECTION 2
    !----------------------------------------------------------------

    SELECT CASE(grib%product_type)
    CASE(PRODUCT_TYPE_OBS)
      ! Observation
      grib%sec2%setLocalDefinition    = GRIB_TRUE    ! Lokale Sektion anschalten
      grib%sec2%localDefinitionNumber = 125  ! grib2LocalSectionNumber.78.table: 125 Observational radar volume scan
      ! Kind of radar device:
      grib%sec2%localTypeOfRadar = 226  ! Local 4.3.table: 226 Dual Pol Doppler Radar
    CASE(PRODUCT_TYPE_SIM)

      ! Synthetisches Produkt
      grib%sec2%setLocalDefinition = GRIB_TRUE         ! Lokale Sektion anschalten

      SELECT CASE (grib_locinfo%typeOfGeneratingProcess)
      CASE (TGP_ENSEMBLE)
        ! Synthetisches Produkt im Ensemble-System
        grib%sec2%localDefinitionNumber       = 123  ! grib2LocalSectionNumber.78.table: 123 Synthetic radar volume scan in ensemble system
        grib%sec2%localTypeOfEnsembleForecast = grib_locinfo%localTypeOfEnsembleForecast
      CASE (TGP_ANALYSIS, TGP_FORECAST, TGP_INITRUN)
        ! Synthetisches Produkt im deterministischen System
        grib%sec2%localDefinitionNumber = 124  ! grib2LocalSectionNumber.78.table: 124 Synthetic radar volume scan in deterministic system
        grib%sec2%localVersionNumber    = 255
      CASE default
        error = 5
        WRITE (errmsg, '(a,i5,a)')  'ERROR '//TRIM(yzroutine)//': unknown grib%sec2%typeOfGeneratingProcess = ', &
             grib_locinfo%typeOfGeneratingProcess, '! Cannot set up grib2 records!'
        RETURN
      END SELECT

      grib%sec2%localHostIdentifier = 0
      grib%sec2%localValidityDateYear   = 0
      grib%sec2%localValidityDateMonth  = 0
      grib%sec2%localValidityDateDay    = 0
      grib%sec2%localValidityDateHour   = 0
      grib%sec2%localValidityDateMinute = 0
      grib%sec2%localValidityDateSecond = 0
      grib%sec2%localInformationNumber  = 0
      grib%sec2%localNumberOfExperiment = grib_locinfo%localNumberOfExperiment
      grib%sec2%localTypeOfRadar = 225  ! Local 4.3.table: 225 EMVORADO
    END SELECT

    ! Lokale Metadaten für Radar-Volumenscans
    grib%sec2%localNumberOfRadarStations   = 1
    grib%sec2%localIndexOfRadarStation     = 1
    grib%sec2%localStationName(:)          = ' '
    grib%sec2%localStationName             = TRIM(rsm%station_name)
    grib%sec2%localCountryId               = rsm%station_id / 1000
    grib%sec2%localNationalStationId       = MODULO(rsm%station_id, 1000)
    ! Breite und Länge müssen in 10^-6 Grad umgerechnet werden
    if (rsm%lat < 0.0_dp) then
      ! set the most significant bit to 1, as required by the grib2 standard for
      ! negative latitudes
      grib%sec2%localStationLatitude = IBSET(NINT( rsm%lat * 1.0E+6_dp ), BIT_SIZE(grib%sec2%localStationLatitude)-1)
    else
      grib%sec2%localStationLatitude       = NINT( rsm%lat * 1.0E+6_dp )
    endif
    grib%sec2%localStationLongitude        = NINT( MODULO(rsm%lon, 360.0_dp) * 1.0E+6_dp )
    ! Höhe der Radar-Station über Meeresspiegel:
    CALL get_scaling( val              = rsm%alt_msl,                                       &
         &            scale_factor     = grib%sec2%localScaleFactorOfStationHeightAboveMSL, &
         &            scaled_val       = grib%sec2%localScaledValueOfStationHeightAboveMSL, &
         &            max_scale_factor = 8                                                  )
    ! Antenna elevation angle:
    CALL get_scaling( val              = rsm%el_arr(ind_ele),                               &
         &            scale_factor     = grib%sec2%localScaleFactorOfAntennaElevationAngle, &
         &            scaled_val       = grib%sec2%localScaledValueOfAntennaElevationAngle, &
         &            max_scale_factor = 8                                                  )
    ! Die Dauer des Radar-Volumenscans wird in Minuten angegeben
    grib%sec2%localIndicatorOfUnitOfTimeRange = 0 ! 4.4.table: 0 Minute
    grib%sec2%localAccumulationInterval       = 5
    grib%sec2%localOperatingMode              = grib%localOperatingMode_in
    grib%sec2%localQualityControlIndicator    = 1
    grib%sec2%localClutterFilterIndicator     = 1

    ! Liste zusätzlicher skalierter Metadaten (siehe lokale Tabelle: 2.120.table)
    grib%sec2%localNumberOfScaledMetadata = 14


    ALLOCATE( grib%sec2%localIndexOfScaledMetadata(grib%sec2%localNumberOfScaledMetadata),     &
         &    grib%sec2%localScaledMetadataTableNumber(grib%sec2%localNumberOfScaledMetadata), &
         &    grib%sec2%localScaledMetadataId(grib%sec2%localNumberOfScaledMetadata),          &
         &    grib%sec2%localScaleFactorOfMetadata(grib%sec2%localNumberOfScaledMetadata),     &
         &    grib%sec2%localScaledValueOfMetadata(grib%sec2%localNumberOfScaledMetadata),     &
         &    grib%sec2%scaledMetadataList(grib%sec2%localNumberOfScaledMetadata)              )

    ! List of descriptions of the additional local scaled metadata (DO NOT CHANGE!):
    grib%sec2%scaledMetadataList(:) %meaning(:) = ' '
    grib%sec2%scaledMetadataList(:) %codeFigure = 0
    grib%sec2%scaledMetadataList(1) %meaning    = "Undetected value                          "
    grib%sec2%scaledMetadataList(1) %value      = zero_value
    grib%sec2%scaledMetadataList(2) %meaning    = "Missing value                             "
    grib%sec2%scaledMetadataList(2) %value      = miss_value
    grib%sec2%scaledMetadataList(3) %meaning    = "Station height above mean sea level - mod "
    grib%sec2%scaledMetadataList(3) %value      = rsm%alt_msl
    grib%sec2%scaledMetadataList(4) %meaning    = "Station height above mean sea level - true"
    grib%sec2%scaledMetadataList(4) %value      = rsm%alt_msl_true
    grib%sec2%scaledMetadataList(5) %meaning    = "Reflectivity calibration constant         "
    grib%sec2%scaledMetadataList(5) %value      = 0.0_dp
    grib%sec2%scaledMetadataList(6) %meaning    = "Reference reflectivity for echo top       "
    grib%sec2%scaledMetadataList(6) %value      = 0.0_dp
    grib%sec2%scaledMetadataList(7) %meaning    = "Extended Nyquist                          "
    grib%sec2%scaledMetadataList(7) %VALUE      = rsm%ext_nyq(ind_ele,itime)
    grib%sec2%scaledMetadataList(8) %meaning    = "High Nyquist                              "
!!$    grib%sec2%scaledMetadataList(8) %value      = rsm%high_nyq(ind_ele)
    grib%sec2%scaledMetadataList(8) %value      = rsm%ext_nyq(ind_ele,itime)
    grib%sec2%scaledMetadataList(9) %meaning    = "Dual PRF ratio                            "
    grib%sec2%scaledMetadataList(9) %value      = rsm%dualprf_ratio(ind_ele)
    grib%sec2%scaledMetadataList(10)%meaning    = "Range gate length                         "
    grib%sec2%scaledMetadataList(10)%value      = rsm%ra_inc
    grib%sec2%scaledMetadataList(11)%meaning    = "Number of ranges averaged over            "
    grib%sec2%scaledMetadataList(11)%value      = REAL(rsm%num_gates, dp)
    grib%sec2%scaledMetadataList(12)%meaning    = "Number of pulses averaged over            "
    grib%sec2%scaledMetadataList(12)%value      = REAL(rsm%num_pulses, dp)
    grib%sec2%scaledMetadataList(13)%meaning    = "ppiStartAzimuth                           "
    grib%sec2%scaledMetadataList(13)%value      = rsm%az_start
    grib%sec2%scaledMetadataList(14)%meaning    = "ppiConstantElevation                      "
    grib%sec2%scaledMetadataList(14)%value      = rsm%el_arr(ind_ele)


    ! For later, set the codeFigure (ID in the 2.120.table) explicitly in grib%sec2%scaledMetadataList()%codeFigure above
    ! and eliminate the search over meanings in set_scaled_metadata_list(), which is based on the descriptions
    ! in grib%tab2, which contains the descriptions and the codeFigures.

    CALL set_scaled_metadata_list( number           = grib%sec2%localNumberOfScaledMetadata,    &
         &                         meta_list        = grib%sec2%scaledMetadataList,             &
         &                         table            = grib%tab2,                                &
         &                         md_index         = grib%sec2%localIndexOfScaledMetadata,     &
         &                         md_table_number  = grib%sec2%localScaledMetadataTableNumber, &
         &                         md_id            = grib%sec2%localScaledMetadataId,          &
         &                         md_scale_factor  = grib%sec2%localScaleFactorOfMetadata,     &
         &                         md_scaled_value  = grib%sec2%localScaledValueOfMetadata      )

    DEALLOCATE( grib%sec2%scaledMetadataList )

    !----------------------------------------------------------------
    !                  GRID DEFINITION SECTION 3
    !----------------------------------------------------------------

    ALLOCATE( grib%sec3%startingAzimuth(n_azimuth), &
         &    grib%sec3%azimuthalWidth(n_azimuth)   )

    grib%sec3%gridDefinitionTemplateNumber = 120  ! 3.1.table: Azimuth-range projection
    ! Keys für Template 120:
    ! * Nb - number of data bins along radials (A data bin is a data point representing the volume centred on it)
    grib%sec3%numberOfDataBinsAlongRadials = n_range
    ! * Nr - number of radials
    grib%sec3%numberOfRadials = n_azimuth
    ! * La1 - latitude of centre point
    !  (92.1.6 Latitude, longitude and angle values shall be in units of 10-6 degree, except for specific
    !          cases explicitly stated in some grid definitions.)
    !  (92.1.7 The latitude values shall be limited to the range 0 to 90 degrees inclusive.
    !          The orientation shall be north latitude positive, south latitude negative.
    !          Bit 1 is set to 1 to indicate south latitude.)
    if (rsm%lat < 0.0_dp) then
      grib%sec3%latitudeOfCenterPoint = IBSET(NINT( rsm%lat * 1.0E+6_dp ), BIT_SIZE(grib%sec3%latitudeOfCenterPoint)-1)
    else
      grib%sec3%latitudeOfCenterPoint = NINT( rsm%lat * 1.0E+6_dp )
    endif
    ! * Lo1 - longitude of centre point
    !  (92.1.8 The longitude values shall be limited to the range 0 to 360 degrees inclusive.
    !          The orientation shall be east longitude positive, with only positive values being used.)
    grib%sec3%longitudeOfCenterPoint = NINT( MODULO(rsm%lon, 360.0_dp) * 1.0E+6_dp )
    ! * Dx - spacing of bins along radials
    !  (Im Template wird eine Einheit nicht explizit angegeben.
    !   Wir nehmen hier die Einheit "Meter" an.)
    grib%sec3%spacingOfBinsAlongRadials = NINT( rsm%ra_inc )
    ! * Dstart - offset from origin to inner bound
    !  (Einheit: siehe Dx.)
    grib%sec3%offsetFromOriginToInnerBound = NINT( rsm%ra_inc )  ! ra_start = ra_inc!
    ! * Scanning mode flag:
    !   Wir gehen von folgender Situation aus: reflectivity(Ni=range,Nj=azimuth,jrecord).
    grib%sec3%iScansNegatively      = GRIB_FALSE  ! range läuft in positiver Richtung
    grib%sec3%jScansPositively      = GRIB_TRUE   ! azimuth läuft in positiver Richtung = Uhrzeigersinn!!!
    grib%sec3%jPointsAreConsecutive = GRIB_FALSE  ! ranges sind zusammenhängend
    ! * Octets 40-(39+4Nr) : For each of Nr radials:
    DO jazimuth=1, n_azimuth
      ! * Azi - starting azimuth, degrees x 10 (degrees as north) -- WE TAKE THE CENTER OF EACH AZIMUT INTERVAL, NOT THE START!
      grib%sec3%startingAzimuth(jazimuth) = NINT( (rsm%az_start+(jazimuth-1)*rsm%az_inc) * 10._dp)
      ! * Adelta - azimuthal width, degrees x 100 (+ clockwise, - counterclockwise)
      grib%sec3%azimuthalWidth(jazimuth) = NINT( rsm%az_inc * 100._dp )
    ENDDO  !jazimuth

    !----------------------------------------------------------------
    !                 PRODUCT DEFINITION SECTION 4
    !----------------------------------------------------------------

    ! (Spätestens) ab hier müssen wir zwischen einer Beobachtung und einer Simulation unterscheiden
    SELECT CASE(grib%product_type)
    CASE(PRODUCT_TYPE_SIM)

      ! Synthetic radar product:

      SELECT CASE (grib_locinfo%typeOfGeneratingProcess)
      CASE (TGP_ENSEMBLE)
        grib%sec4%productDefinitionTemplateNumber = 1   ! 4.1.table:............ Individual ensemble forecast, control and perturbed,
                                                          !                        at a horizontal level or in a horizontal layer at a point in time
        grib%sec4%typeOfEnsembleForecast          = grib_locinfo%typeOfEnsembleForecast
        grib%sec4%perturbationNumber              = grib_locinfo%perturbationNumber
        grib%sec4%numberOfForecastsInEnsemble     = grib_locinfo%numberOfForecastsInEnsemble
      CASE (TGP_ANALYSIS, TGP_FORECAST, TGP_INITRUN)
        grib%sec4%productDefinitionTemplateNumber = 0   ! 4.0.table:............ Analysis or forecast at a horizontal level or in a horizontal layer at a point in time
      END SELECT

!      grib%sec4%parameterCategory       = 15    ! already set in calling program
!      grib%sec4%parameterNumber         = 1     ! already set in calling program
      grib%sec4%typeOfGeneratingProcess = grib_locinfo%typeOfGeneratingProcess
      grib%sec4%backgroundProcess       = grib_locinfo%backgroundProcess
      grib%sec4%generatingProcessIdentifier = grib_locinfo%generatingProcessIdentifier
      grib%sec4%hoursAfterDataCutoff   = 0  ! Not used
      grib%sec4%minutesAfterDataCutoff = 0  ! Not used

      grib%sec4%indicatorOfUnitOfTimeRange = 0 ! 4.4.table: 0 Minute
      grib%sec4%forecastTime               = grib%forecastminutes

      grib%sec4%typeOfFirstFixedSurface         = 198  ! Lokale Tabelle 4.5.table: 198 radarElev radarElev Radar antenna elevation angle (degree)
      grib%sec4%scaleFactorOfFirstFixedSurface  = 1
      ! No negative values allowed, because these keys are unsigned ints! Negative elevations will therefore show up as ele+360.
      grib%sec4%scaledValueOfFirstFixedSurface  = NINT( MODULO(rsm%el_arr(ind_ele),360.0_dp) * &
                                                    10 ** grib%sec4%scaleFactorOfFirstFixedSurface )
      grib%sec4%typeOfSecondFixedSurface        = GRIB_MISSING  ! Keine Schicht, daher MISSING = 255


    CASE(PRODUCT_TYPE_OBS)

      ! Radar observation:

      grib%sec4%productDefinitionTemplateNumber = 20  ! 4.20.table:............. radar product
!      grib%sec4%parameterCategory               = 15    ! already set in calling program
!      grib%sec4%parameterNumber                 = 1     ! already set in calling program
      grib%sec4%typeOfGeneratingProcess         = 8   ! 4.3.table:.............. Observation
      grib%sec4%numberOfRadarSitesUsed          = grib%sec2%localNumberOfRadarStations
      grib%sec4%indicatorOfUnitOfTimeRange      = 0   ! 4.4.table:.............. Minute
      if (rsm%lat < 0.0_dp) then
        grib%sec4%siteLatitude = IBSET( NINT( rsm%lat * 1.0E+6_dp ), BIT_SIZE(grib%sec4%siteLatitude)-1 )
      else 
        grib%sec4%siteLatitude                  = NINT( rsm%lat * 1.0E+6_dp )  ! Site latitude (in 10-6 degree)
      endif
      grib%sec4%siteLongitude                   = NINT( MODULO(rsm%lon, 360.0_dp) * 1.0E+6_dp )  ! Site longitude (in 10-6 degree)
      grib%sec4%siteElevation                   = NINT( rsm%alt_msl )          ! Site elevation (meters)
      grib%sec4%siteId                          = rsm%station_id
      ! Unfortunately no negative values possible, therefore folding with 360:
      ! Unfortunately unsigned byte data type, so folding with 25.5:
      grib%sec4%constantAntennaElevationAngle   = NINT( MODULO(rsm%el_arr(ind_ele),25.5_dp) * 10._dp ) ! Constant antenna elevation angle (tenths of degree true)
      grib%sec4%rangeBinSpacing                 = NINT( rsm%ra_inc )          ! Range bin spacing (meters)
      grib%sec4%radialAngularSpacing            = NINT( rsm%az_inc * 10._dp ) ! Radial angular spacing (tenths of degree true)
      grib%sec4%operatingMode                   = grib%localOperatingMode_in
      grib%sec4%reflectivityCalibrationConstant = 0        ! tenth's of dB
      grib%sec4%qualityControlIndicator         = 1
      grib%sec4%clutterFilterIndicator          = 1
      grib%sec4%accumulationInterval            = 5
      grib%sec4%referenceReflectivityForEchoTop = 0  ! instead of GRIB_MISSING to be consistent with local section value


    END SELECT

    !----------------------------------------------------------------
    !                DATA REPRESENTATION SECTION 5
    !----------------------------------------------------------------

    grib%sec5%bitsPerValue              = 16
    grib%sec5%packingType(:)            = ' '
    grib%sec5%packingType               = TRIM(grib%packingtype)
    grib%sec5%typeOfOriginalFieldValues = 0  ! 5.1.table: Floating point

    !----------------------------------------------------------------
    !                      BIT-MAP SECTION 6
    !----------------------------------------------------------------

    ! Nothing to set here

    !----------------------------------------------------------------
    !                       DATA SECTION 7
    !----------------------------------------------------------------

    ! Nothing to prepare here, data values have been put to grib-type already above!

    IF (lverbose) WRITE (*,*) 'done with '//TRIM(yzroutine), ' on proc ', my_radar_id

  END SUBROUTINE fill_grib_input_volscan

  !===============================================================================
  !===============================================================================

  !===============================================================================
  !
  ! Put the content of the local section definition file for single-site radar scans
  !
  !    /grib2/tables/local/edzw/1/2.78.120.table
  !
  ! to the tab2 structure
  !
  !===============================================================================


  SUBROUTINE get_loc_sec_tab(tab2, verbose)

    IMPLICIT NONE

    ! in/out
    TYPE(t_grib_loc_sec_tab),   INTENT(INOUT) :: tab2
    LOGICAL,          OPTIONAL, INTENT(IN)    :: verbose

    INTEGER, PARAMETER  :: size_loc_tab = 15
    INTEGER :: jvalid

    !------------------------------------------------------

    ! Content of the file ${ECCODES_DEFINITION_PATH}/grib2/tables/local/edzw/1/2.120.table
    ! =======================================================================================
    !
    ! # This file is automatically generated, don't edit!
    ! #
    ! # Table 2.78.120 contains (additional) scaled metadata
    ! # for single radar volume scans for local use
    ! #
    ! # IMPORTANT "Meaning" must not contain numbers and parentheses, except in the "(unit)"!!!
    ! #
    ! # 2*(Code figure) Meaning
    ! #
    ! 0 0 Undetected value
    ! 1 1 Missing value
    ! # 2-5 Reserved (for further indicators of special data points)
    ! 6 6 Station height above mean sea level - mod (m)
    ! 7 7 Station height above mean sea level - true (m)
    ! 8 8 Reflectivity calibration constant (dB)
    ! 9 9 Reference reflectivity for echo top (dB)
    ! 10 10 Extended Nyquist (s-1)
    ! 11 11 High Nyquist (s-1)
    ! 12 12 Dual PRF ratio (proportion)
    ! 13 13 Range gate length (m)
    ! 14 14 Number of ranges averaged over (numeric)
    ! 15 15 Number of pulses averaged over (numeric)
    ! 16 16 ppiStartAzimuth (degree)
    ! 17 17 ppiConstantElevation (degree)
    ! 255 255 Missing

    IF (ALLOCATED(tab2%codeFigure)) DEALLOCATE(tab2%codeFigure)
    IF (ALLOCATED(tab2%meaning))    DEALLOCATE(tab2%meaning)

    ALLOCATE( tab2%codeFigure(size_loc_tab), &
      &       tab2%meaning(size_loc_tab)     )

    tab2%nentry      = size_loc_tab
    tab2%tableNumber = 120

    tab2%codeFigure(:)  = 255;    tab2%meaning(:)(:) = " "
    tab2%codeFigure(1)  = 0  ;    tab2%meaning(1)  = "Undetected value"
    tab2%codeFigure(2)  = 1  ;    tab2%meaning(2)  = "Missing value"
    tab2%codeFigure(3)  = 6  ;    tab2%meaning(3)  = "Station height above mean sea level - mod"
    tab2%codeFigure(4)  = 7  ;    tab2%meaning(4)  = "Station height above mean sea level - true"
    tab2%codeFigure(5)  = 8  ;    tab2%meaning(5)  = "Reflectivity calibration constant"
    tab2%codeFigure(6)  = 9  ;    tab2%meaning(6)  = "Reference reflectivity for echo top"
    tab2%codeFigure(7)  = 10 ;    tab2%meaning(7)  = "Extended Nyquist"
    tab2%codeFigure(8)  = 11 ;    tab2%meaning(8)  = "High Nyquist"
    tab2%codeFigure(9)  = 12 ;    tab2%meaning(9)  = "Dual PRF ratio"
    tab2%codeFigure(10) = 13 ;    tab2%meaning(10) = "Range gate length"
    tab2%codeFigure(11) = 14 ;    tab2%meaning(11) = "Number of ranges averaged over"
    tab2%codeFigure(12) = 15 ;    tab2%meaning(12) = "Number of pulses averaged over"
    tab2%codeFigure(13) = 16 ;    tab2%meaning(13) = "ppiStartAzimuth"
    tab2%codeFigure(14) = 17 ;    tab2%meaning(14) = "ppiConstantElevation"
    tab2%codeFigure(15) = 255;    tab2%meaning(15) = "Missing"

!!$    IF(PRESENT(verbose)) THEN
!!$      IF(verbose) THEN
!!$        ! nvalid sollte an dieser Stelle noch die Anzahl
!!$        ! der gültigen Einträge enthalten
!!$        WRITE(*,*) " "
!!$        WRITE(*,*) "Anzahl gueltiger Eintraege in Tabelle 120: ", tab2%nentry
!!$        WRITE(*,*) "Tabellennummer: ", tab2%tableNumber
!!$        WRITE(*,*) "----------------------------------------------------------"
!!$        DO jvalid=1, size_loc_tab
!!$          WRITE(*,*) tab2%codeFigure(jvalid), TRIM(tab2%meaning(jvalid))
!!$        ENDDO
!!$        WRITE(*,*) "----------------------------------------------------------"
!!$        WRITE(*,*) " "
!!$      ENDIF
!!$    ENDIF

  END SUBROUTINE get_loc_sec_tab

  !===============================================================================

  SUBROUTINE get_scaling_real(val, scale_factor, scaled_val, max_scale_factor)

    IMPLICIT NONE

    ! in/out
    REAL(dp),          INTENT(IN)    :: val
    INTEGER,           INTENT(INOUT) :: scale_factor
    INTEGER,           INTENT(INOUT) :: scaled_val
    INTEGER, OPTIONAL, INTENT(IN)    :: max_scale_factor

    ! local
    REAL(dp) :: scaled_val_real
    INTEGER  :: scale_factor_max

    REAL(dp), PARAMETER :: eps = 0.0001_dp       ! Z.B. 0.800012232
                                                 !         |---------> wegschneiden
    REAL(dp), PARAMETER :: threshold = 1.0E-8_dp
    INTEGER,  PARAMETER :: def_max_scale_factor = 8

    !------------------------------------------------------

    !-------------------------------------------
    ! Vorbereitung
    !-------------------------------------------

    IF(PRESENT(max_scale_factor)) THEN
      scale_factor_max = max_scale_factor
    ELSE
      scale_factor_max = def_max_scale_factor
    ENDIF

    scale_factor    = 0
    scaled_val_real = val

    IF(ABS(val) < threshold) THEN

      scaled_val = INT( val )

    ELSE

      !-------------------------------------------
      ! Skalierung
      !-------------------------------------------

      SCALE_LOOP : DO

        IF (scaled_val_real > REAL(HUGE(1),dp)) THEN
          IF (scale_factor == 0) THEN
            scale_factor = -(INT(LOG10(scaled_val_real))-def_max_scale_factor)
            scaled_val_real = scaled_val_real * 10.0_dp**scale_factor
          ELSE
            scale_factor    = scale_factor - 1
            scaled_val_real = scaled_val_real * 0.1_dp
          END IF
          EXIT SCALE_LOOP
        END IF

        IF( ( ABS(scaled_val_real) > 0._dp .AND. ABS( scaled_val_real - INT(scaled_val_real) ) < eps ) .OR. &
             & scale_factor == scale_factor_max ) THEN
          EXIT SCALE_LOOP
        END IF

        scale_factor = scale_factor + 1

        scaled_val_real = scaled_val_real * 10._dp


      END DO SCALE_LOOP

      IF (scaled_val_real > REAL(HUGE(1),dp)) THEN
        WRITE (*,'(a,es14.5,a)') 'Problem in get_scaling_real(): scaled value ', &
             scaled_val_real, ' too large for standard 4-byte integer!'
      END IF

      scaled_val = INT( MIN(scaled_val_real, REAL(HUGE(1),dp)) )

    END IF

!    WRITE (*,*) 'get_scaling_real(): unscaled value: ', val, &
!         '   scaled value and scale factor: ', scaled_val, scale_factor

  END SUBROUTINE get_scaling_real

  !===============================================================================

  SUBROUTINE get_scaling_int(val, scale_factor, scaled_val, max_scale_factor)

    IMPLICIT NONE

    ! in/out
    INTEGER,           INTENT(IN)    :: val
    INTEGER,           INTENT(INOUT) :: scale_factor
    INTEGER,           INTENT(INOUT) :: scaled_val
    INTEGER, OPTIONAL, INTENT(IN)   :: max_scale_factor

    !------------------------------------------------------

    scale_factor = 0
    scaled_val   = val

  END SUBROUTINE get_scaling_int

  !===============================================================================

  SUBROUTINE set_scaled_metadata_list(number, meta_list, table, md_index, md_table_number, &
    &                                 md_id, md_scale_factor, md_scaled_value)

    IMPLICIT NONE

    ! in/out
    INTEGER,                          INTENT(IN)    :: number               ! Anzahl zusätzlicher skalierter Metadaten
    type(t_scaledMetaDataList),       INTENT(IN)    :: meta_list(:)         ! Liste mit den Metadaten
    TYPE(t_grib_loc_sec_tab),         INTENT(IN)    :: table                ! Tabelle mit den Metadaten
    INTEGER,                          INTENT(INOUT) :: md_index(:)          ! Listen-Zähler der Metadaten
    INTEGER,                          INTENT(INOUT) :: md_table_number(:)   ! Tabellen-Nummer
    INTEGER,                          INTENT(INOUT) :: md_id(:)             ! Code-Figure der Metadaten in der Metadaten-Tabelle
    INTEGER,                          INTENT(INOUT) :: md_scale_factor(:)   ! Skalen-Faktor der Metadaten
    INTEGER,                          INTENT(INOUT) :: md_scaled_value(:)   ! Skalierter Wert der Metadaten

    ! local
    INTEGER :: j, k, idx
    LOGICAL :: found

    CHARACTER(LEN=MAX_LEN_TAB_ENTRY) :: description, tab_entry

    !------------------------------------------------------

    METADATA_LOOP : DO j=1, number

      md_index(j)        = j
      md_table_number(j) = table%tableNumber

      ! Suche Eintrag in Tabelle
      description = TRIM(ADJUSTL(meta_list(j)%meaning))
      description = to_lower_or_upper(description, "lower")
      found       = .FALSE.
      TABLE_LOOP : DO k=1, table%nentry

        tab_entry = to_lower_or_upper(table%meaning(k), "lower")
        idx = INDEX(tab_entry, TRIM(description))
        IF(idx > 0) THEN
          found    = .TRUE.
          md_id(j) = table%codeFigure(k)
          EXIT TABLE_LOOP
        ENDIF

      ENDDO TABLE_LOOP

      IF (.NOT. found) THEN
        WRITE(*,*) 'ERROR set_scaled_metadata_list: difference between |'//tab_entry//'| and |'//TRIM(description)//'|'
        CALL abort_run(my_radar_id, 10377, &
             '"'//TRIM(description)//'" in Tabelle nicht gefunden!', &
             'set_scaled_metadata_list')
      END IF

      CALL get_scaling(val=meta_list(j)%value, scale_factor=md_scale_factor(j), scaled_val=md_scaled_value(j))

    ENDDO METADATA_LOOP

  END SUBROUTINE set_scaled_metadata_list

  !===============================================================================
  !===============================================================================

#ifdef GRIBAPI

  SUBROUTINE write_grib_to_file ( grib_outfilename,  &
       &                          force_overwrite,   &
       &                          grib,              &
       &                          lverbose,          &
       &                          error,             &
       &                          errmsg             )

    ! in/out
    CHARACTER(LEN=*), INTENT(in)    :: grib_outfilename
    LOGICAL,          INTENT(in)    :: force_overwrite
    TYPE(t_grib),     INTENT(inout) :: grib
    LOGICAL,          INTENT(in)    :: lverbose
    INTEGER         , INTENT(OUT)   :: error
    CHARACTER(len=*), INTENT(OUT)   :: errmsg

    ! local
    CHARACTER(len=*), PARAMETER :: yzroutine = 'write_grib_to_file'


    INTEGER                     :: igribsample
    INTEGER                     :: igribclone
    INTEGER                     :: ioutfile

    INTEGER                     :: jvalue, nvalue

    INTEGER                     :: iunit, istat

    LOGICAL                     :: lopened, lexist
    INTEGER                     :: ierr, ilocerr
    CHARACTER(len=1)            :: cmode
    INTEGER, PARAMETER          :: GRIB_MISSING   = 255
    INTEGER, PARAMETER          :: GRIB_TRUE      = 1
    INTEGER, PARAMETER          :: GRIB_FALSE     = 0
    INTEGER, PARAMETER          :: GRIB_EDITION_2 = 2

    !------------------------------------------------------

    ierr      = 0
    error     = 0
    errmsg(:) = ' '

    IF (lverbose) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    INQUIRE (file=TRIM(grib_outfilename), exist=lexist)
    IF (.NOT.lexist .OR. force_overwrite) THEN
      cmode = 'w'   ! write
      IF (lverbose) WRITE (*,'(A)') 'Creating '//TRIM(grib_outfilename)
    ELSE
      cmode = 'a'   ! append
      IF (lverbose) WRITE (*,'(A)') 'Appending to '//TRIM(grib_outfilename)
    END IF

    ! Open grib file:
    CALL codes_open_file(ioutfile, grib_outfilename, cmode)

    IF (ilocerr /= GRIB_SUCCESS) THEN
      errmsg(:) = ' '
      CALL codes_get_error_string(ilocerr, errmsg)
      errmsg = 'ERROR '//TRIM(yzroutine)//' opening '//TRIM(grib_outfilename)//': '//TRIM(errmsg)
      error = 4
      RETURN
    END IF

    ! Load GRIB-Template from eccodes tables:
    CALL codes_grib_new_from_samples(igribsample, TRIM(grib%sample_file), ilocerr)

    IF (ilocerr /= GRIB_SUCCESS) THEN
      errmsg(:) = ' '
      CALL codes_get_error_string(ilocerr, errmsg)
      errmsg = 'ERROR '//TRIM(yzroutine)//' reading grib sample file '//TRIM(grib%sample_file)//': '//TRIM(errmsg)
      error = 3
      RETURN
    END IF

    ! Clone template:
    CALL codes_clone(igribsample,igribclone)

    !----------------------------------------------------------------
    !                     INDICATOR SECTION 0
    !----------------------------------------------------------------

    CALL codes_set(igribclone,'discipline',grib%sec0%discipline,ilocerr)
      CALL check_codes_err(ilocerr, 'discipline', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'editionNumber',grib%sec0%editionNumber,ilocerr)
      CALL check_codes_err(ilocerr, 'editionNumber', TRIM(yzroutine), increrr=ierr)

    !----------------------------------------------------------------
    !                   IDENTIFICATION SECTION 1
    !----------------------------------------------------------------

    CALL codes_set(igribclone,'tablesVersion',grib%sec1%tablesVersion,ilocerr)
      CALL check_codes_err(ilocerr, 'tablesVersion', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'productionStatusOfProcessedData',grib%sec1%productionStatusOfProcessedData,ilocerr)
      CALL check_codes_err(ilocerr, 'productionStatusOfProcessedData', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'typeOfProcessedData',grib%sec1%typeOfProcessedData,ilocerr)
      CALL check_codes_err(ilocerr, 'typeOfProcessedData', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'centre',grib%sec1%centre,ilocerr)
      CALL check_codes_err(ilocerr, 'centre', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'subCentre',grib%sec1%subCentre,ilocerr)
      CALL check_codes_err(ilocerr, 'subCentre', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'significanceOfReferenceTime',grib%sec1%significanceOfReferenceTime,ilocerr)
      CALL check_codes_err(ilocerr, 'significanceOfReferenceTime', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'year',grib%sec1%year,ilocerr)
      CALL check_codes_err(ilocerr, 'year', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'month',grib%sec1%month,ilocerr)
      CALL check_codes_err(ilocerr, 'month', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'day',grib%sec1%day,ilocerr)
      CALL check_codes_err(ilocerr, 'day', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'hour',grib%sec1%hour,ilocerr)
      CALL check_codes_err(ilocerr, 'hour', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'minute',grib%sec1%minute,ilocerr)
      CALL check_codes_err(ilocerr, 'minute', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'second',grib%sec1%second,ilocerr)
      CALL check_codes_err(ilocerr, 'second', TRIM(yzroutine), increrr=ierr)

    !----------------------------------------------------------------
    !                     LOCAL USE SECTION 2
    !----------------------------------------------------------------

    ! The local section has to be deleted and re-created,
    ! to avoid difficulties with some remnants in the sample file (e.g. section2Padding)
    CALL codes_set(igribclone,'setLocalDefinition',GRIB_FALSE,ilocerr)
      CALL check_codes_err(ilocerr, 'setLocalDefinition', TRIM(yzroutine), increrr=ierr)

    CALL codes_set(igribclone,'setLocalDefinition',grib%sec2%setLocalDefinition,ilocerr)
      CALL check_codes_err(ilocerr, 'setLocalDefinition', TRIM(yzroutine), increrr=ierr)

    IF(grib%sec2%setLocalDefinition == GRIB_TRUE) THEN

      CALL codes_set(igribclone,'localDefinitionNumber',grib%sec2%localDefinitionNumber,ilocerr)
        CALL check_codes_err(ilocerr, 'localDefinitionNumber', TRIM(yzroutine), increrr=ierr)

      IF(grib%product_type == PRODUCT_TYPE_SIM) THEN

        CALL codes_set(igribclone,'localHostIdentifier',grib%sec2%localHostIdentifier,ilocerr)
          CALL check_codes_err(ilocerr, 'localHostIdentifier', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localCreationDateYear',grib%sec2%localCreationDateYear,ilocerr)
          CALL check_codes_err(ilocerr, 'localCreationDateYear', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localCreationDateMonth',grib%sec2%localCreationDateMonth,ilocerr)
          CALL check_codes_err(ilocerr, 'localCreationDateMonth', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localCreationDateDay',grib%sec2%localCreationDateDay,ilocerr)
          CALL check_codes_err(ilocerr, 'localCreationDateDay', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localCreationDateHour',grib%sec2%localCreationDateHour,ilocerr)
          CALL check_codes_err(ilocerr, 'localCreationDateHour', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localCreationDateMinute',grib%sec2%localCreationDateMinute,ilocerr)
          CALL check_codes_err(ilocerr, 'localCreationDateMinute', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localCreationDateSecond',grib%sec2%localCreationDateSecond,ilocerr)
          CALL check_codes_err(ilocerr, 'localCreationDateSecond', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localValidityDateYear',grib%sec2%localValidityDateYear,ilocerr)
         CALL check_codes_err(ilocerr, 'localValidityDateYear', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localValidityDateMonth',grib%sec2%localValidityDateMonth,ilocerr)
          CALL check_codes_err(ilocerr, 'localValidityDateMonth', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localValidityDateDay',grib%sec2%localValidityDateDay,ilocerr)
          CALL check_codes_err(ilocerr, 'localValidityDateDay', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localValidityDateHour',grib%sec2%localValidityDateHour,ilocerr)
          CALL check_codes_err(ilocerr, 'localValidityDateHour', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localValidityDateMinute',grib%sec2%localValidityDateMinute,ilocerr)
          CALL check_codes_err(ilocerr, 'localValidityDateMinute', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localValidityDateSecond',grib%sec2%localValidityDateSecond,ilocerr)
          CALL check_codes_err(ilocerr, 'localValidityDateSecond', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localInformationNumber',grib%sec2%localInformationNumber,ilocerr)
          CALL check_codes_err(ilocerr, 'localInformationNumber', TRIM(yzroutine), increrr=ierr)
        CALL codes_set(igribclone,'localNumberOfExperiment',grib%sec2%localNumberOfExperiment,ilocerr)
          CALL check_codes_err(ilocerr, 'localNumberOfExperiment', TRIM(yzroutine), increrr=ierr)

        SELECT CASE (grib%sec4%typeOfGeneratingProcess)
        CASE (TGP_ENSEMBLE)
          CALL codes_set(igribclone,'localTypeOfEnsembleForecast',grib%sec2%localTypeOfEnsembleForecast,ilocerr)
          CALL check_codes_err(ilocerr, 'localTypeOfEnsembleForecast', TRIM(yzroutine), increrr=ierr)
        CASE (TGP_ANALYSIS, TGP_FORECAST, TGP_INITRUN)
          CALL codes_set(igribclone,'localVersionNumber',grib%sec2%localVersionNumber,ilocerr)
          CALL check_codes_err(ilocerr, 'localVersionNumber', TRIM(yzroutine), increrr=ierr)
        END SELECT

      ENDIF  !IF(product_type == PRODUCT_TYPE_SIM)

      CALL codes_set(igribclone,'localNumberOfRadarStations',grib%sec2%localNumberOfRadarStations,ilocerr)
        CALL check_codes_err(ilocerr, 'localNumberOfRadarStations', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localIndexOfRadarStation',grib%sec2%localIndexOfRadarStation,ilocerr)
        CALL check_codes_err(ilocerr, 'localIndexOfRadarStation', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localCountryId',grib%sec2%localCountryId,ilocerr)
        CALL check_codes_err(ilocerr, 'localCountryId', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localNationalStationId',grib%sec2%localNationalStationId,ilocerr)
        CALL check_codes_err(ilocerr, 'localNationalStationId', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localStationName',grib%sec2%localStationName(1:4),ilocerr)
        CALL check_codes_err(ilocerr, 'localStationName', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localStationLatitude',grib%sec2%localStationLatitude,ilocerr)
        CALL check_codes_err(ilocerr, 'localStationLatitude', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localStationLongitude',grib%sec2%localStationLongitude,ilocerr)
        CALL check_codes_err(ilocerr, 'localStationLongitude', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaleFactorOfStationHeightAboveMSL', &
                       grib%sec2%localScaleFactorOfStationHeightAboveMSL,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaleFactorOfStationHeightAboveMSL', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaledValueOfStationHeightAboveMSL', &
                       grib%sec2%localScaledValueOfStationHeightAboveMSL,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaledValueOfStationHeightAboveMSL', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaleFactorOfAntennaElevationAngle', &
                       grib%sec2%localScaleFactorOfAntennaElevationAngle,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaleFactorOfAntennaElevationAngle', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaledValueOfAntennaElevationAngle', &
                       grib%sec2%localScaledValueOfAntennaElevationAngle,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaledValueOfAntennaElevationAngle', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localIndicatorOfUnitOfTimeRange',grib%sec2%localIndicatorOfUnitOfTimeRange,ilocerr)
        CALL check_codes_err(ilocerr, 'localIndicatorOfUnitOfTimeRange', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localAccumulationInterval',grib%sec2%localAccumulationInterval,ilocerr)
        CALL check_codes_err(ilocerr, 'localAccumulationInterval', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localTypeOfRadar',grib%sec2%localTypeOfRadar,ilocerr)
        CALL check_codes_err(ilocerr, 'localTypeOfRadar', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localOperatingMode',grib%sec2%localOperatingMode,ilocerr)
        CALL check_codes_err(ilocerr, 'localOperatingMode', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localQualityControlIndicator',grib%sec2%localQualityControlIndicator,ilocerr)
        CALL check_codes_err(ilocerr, 'localQualityControlIndicator', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localClutterFilterIndicator',grib%sec2%localClutterFilterIndicator,ilocerr)
        CALL check_codes_err(ilocerr, 'localClutterFilterIndicator', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localTypeOfDataMaskedByBitmap',grib%sec2%localTypeOfDataMaskedByBitmap,ilocerr)
        CALL check_codes_err(ilocerr, 'localTypeOfDataMaskedByBitmap', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaleFactorOfBitmapValue',grib%sec2%localScaleFactorOfBitmapValue,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaleFactorOfBitmapValue', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaledValueOfBitmapValue',grib%sec2%localScaledValueOfBitmapValue,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaledValueOfBitmapValue', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localNumberOfScaledMetadata',grib%sec2%localNumberOfScaledMetadata,ilocerr)
        CALL check_codes_err(ilocerr, 'localNumberOfScaledMetadata', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localIndexOfScaledMetadata',grib%sec2%localIndexOfScaledMetadata,ilocerr)
        CALL check_codes_err(ilocerr, 'localIndexOfScaledMetadata', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaledMetadataTableNumber',grib%sec2%localScaledMetadataTableNumber,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaledMetadataTableNumber', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaledMetadataId',grib%sec2%localScaledMetadataId,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaledMetadataId', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaleFactorOfMetadata',grib%sec2%localScaleFactorOfMetadata,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaleFactorOfMetadata', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'localScaledValueOfMetadata',grib%sec2%localScaledValueOfMetadata,ilocerr)
        CALL check_codes_err(ilocerr, 'localScaledValueOfMetadata', TRIM(yzroutine), increrr=ierr)

    ENDIF  !IF(grib%sec2%setLocalDefinition == GRIB_TRUE)

    DEALLOCATE( grib%sec2%localIndexOfScaledMetadata,     &
         &      grib%sec2%localScaledMetadataTableNumber, &
         &      grib%sec2%localScaledMetadataId,          &
         &      grib%sec2%localScaleFactorOfMetadata,     &
         &      grib%sec2%localScaledValueOfMetadata      )

    !----------------------------------------------------------------
    !                  GRID DEFINITION SECTION 3
    !----------------------------------------------------------------


    CALL codes_set(igribclone,'gridType','azimuth_range',ilocerr)  ! = grib%sec3%gridDefinitionTemplateNumber = 120
      CALL check_codes_err(ilocerr, 'gridType', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'numberOfDataBinsAlongRadials',grib%sec3%numberOfDataBinsAlongRadials,ilocerr)
      CALL check_codes_err(ilocerr, 'numberOfDataBinsAlongRadials', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'numberOfRadials',grib%sec3%numberOfRadials,ilocerr)
      CALL check_codes_err(ilocerr, 'numberOfRadials', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'latitudeOfCentrePoint',grib%sec3%latitudeOfCenterPoint,ilocerr)
    IF (ilocerr /= GRIB_SUCCESS) THEN
      ! eccodes is older than 2.20.0 where the key was renamed from 'center' to 'centre':
      CALL codes_set(igribclone,'latitudeOfCenterPoint',grib%sec3%latitudeOfCenterPoint,ilocerr)
        CALL check_codes_err(ilocerr, 'latitudeOfCenterPoint', TRIM(yzroutine), increrr=ierr)
    END IF
    CALL codes_set(igribclone,'longitudeOfCentrePoint',grib%sec3%longitudeOfCenterPoint,ilocerr)
    IF (ilocerr /= GRIB_SUCCESS) THEN
      ! eccodes is older than 2.20.0 where the key was renamed from 'center' to 'centre':
      CALL codes_set(igribclone,'longitudeOfCenterPoint',grib%sec3%longitudeOfCenterPoint,ilocerr)
        CALL check_codes_err(ilocerr, 'longitudeOfCenterPoint', TRIM(yzroutine), increrr=ierr)
    END IF
    CALL codes_set(igribclone,'spacingOfBinsAlongRadials',grib%sec3%spacingOfBinsAlongRadials,ilocerr)
      CALL check_codes_err(ilocerr, 'spacingOfBinsAlongRadials', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'offsetFromOriginToInnerBound',grib%sec3%offsetFromOriginToInnerBound,ilocerr)
      CALL check_codes_err(ilocerr, 'offsetFromOriginToInnerBound', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'iScansNegatively',grib%sec3%iScansNegatively,ilocerr)
      CALL check_codes_err(ilocerr, 'iScansNegatively', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'jScansPositively',grib%sec3%jScansPositively,ilocerr)
      CALL check_codes_err(ilocerr, 'jScansPositively', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'jPointsAreConsecutive',grib%sec3%jPointsAreConsecutive,ilocerr)
      CALL check_codes_err(ilocerr, 'jPointsAreConsecutive', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'startingAzimuth',grib%sec3%startingAzimuth,ilocerr)
      CALL check_codes_err(ilocerr, 'startingAzimuth', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'azimuthalWidth',grib%sec3%azimuthalWidth,ilocerr)
      CALL check_codes_err(ilocerr, 'azimuthalWidth', TRIM(yzroutine), increrr=ierr)

    ! Deallozierung von Feldern
    DEALLOCATE( grib%sec3%startingAzimuth, &
         &      grib%sec3%azimuthalWidth   )


    !----------------------------------------------------------------
    !                 PRODUCT DEFINITION SECTION 4
    !----------------------------------------------------------------

    SELECT CASE(grib%product_type)
    CASE(PRODUCT_TYPE_SIM)

      CALL codes_set(igribclone,'deletePV',GRIB_TRUE,ilocerr)
        CALL check_codes_err(ilocerr, 'deletePV', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'productDefinitionTemplateNumber',grib%sec4%productDefinitionTemplateNumber,ilocerr)
        CALL check_codes_err(ilocerr, 'productDefinitionTemplateNumber', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'parameterCategory',grib%sec4%parameterCategory,ilocerr)
        CALL check_codes_err(ilocerr, 'parameterCategory', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'parameterNumber',grib%sec4%parameterNumber,ilocerr)
        CALL check_codes_err(ilocerr, 'parameterNumber', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'typeOfGeneratingProcess',grib%sec4%typeOfGeneratingProcess,ilocerr)
        CALL check_codes_err(ilocerr, 'typeOfGeneratingProcess', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'backgroundProcess',grib%sec4%backgroundProcess,ilocerr)
        CALL check_codes_err(ilocerr, 'backgroundProcess', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'generatingProcessIdentifier',grib%sec4%generatingProcessIdentifier,ilocerr)
        CALL check_codes_err(ilocerr, 'generatingProcessIdentifier', TRIM(yzroutine), increrr=ierr)

      IF (grib%sec4%typeOfGeneratingProcess == TGP_ENSEMBLE) THEN
        CALL codes_set  (igribclone, 'typeOfEnsembleForecast', grib%sec4%typeOfEnsembleForecast, ilocerr)
          CALL check_codes_err(ilocerr, 'typeOfEnsembleForecast', TRIM(yzroutine), increrr=ierr)
        CALL codes_set  (igribclone, 'numberOfForecastsInEnsemble', grib%sec4%numberOfForecastsInEnsemble, ilocerr)
          CALL check_codes_err(ilocerr, 'numberOfForecastsInEnsemble', TRIM(yzroutine), increrr=ierr)
        CALL codes_set  (igribclone, 'perturbationNumber', grib%sec4%perturbationNumber, ilocerr)
          CALL check_codes_err(ilocerr, 'perturbationNumber', TRIM(yzroutine), increrr=ierr)
      END IF

      CALL codes_set(igribclone,'hoursAfterDataCutoff',grib%sec4%hoursAfterDataCutoff,ilocerr)
        CALL check_codes_err(ilocerr, 'hoursAfterDataCutoff', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'minutesAfterDataCutoff',grib%sec4%minutesAfterDataCutoff,ilocerr)
        CALL check_codes_err(ilocerr, 'minutesAfterDataCutoff', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'indicatorOfUnitOfTimeRange',grib%sec4%indicatorOfUnitOfTimeRange,ilocerr)
        CALL check_codes_err(ilocerr, 'indicatorOfUnitOfTimeRange', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'forecastTime',grib%sec4%forecastTime,ilocerr)
        CALL check_codes_err(ilocerr, 'forecastTime', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'typeOfFirstFixedSurface',grib%sec4%typeOfFirstFixedSurface)
        CALL check_codes_err(ilocerr, 'typeOfFirstFixedSurface', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'scaleFactorOfFirstFixedSurface',grib%sec4%scaleFactorOfFirstFixedSurface,ilocerr)
        CALL check_codes_err(ilocerr, 'scaleFactorOfFirstFixedSurface', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'scaledValueOfFirstFixedSurface',grib%sec4%scaledValueOfFirstFixedSurface,ilocerr)
        CALL check_codes_err(ilocerr, 'scaledValueOfFirstFixedSurface', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'typeOfSecondFixedSurface',grib%sec4%typeOfSecondFixedSurface,ilocerr)
        CALL check_codes_err(ilocerr, 'typeOfSecondFixedSurface', TRIM(yzroutine), increrr=ierr)
      !
      CALL codes_set_missing(igribclone,'scaleFactorOfSecondFixedSurface',ilocerr)
        CALL check_codes_err(ilocerr, 'scaleFactorOfSecondFixedSurface', TRIM(yzroutine), increrr=ierr)
      CALL codes_set_missing(igribclone,'scaledValueOfSecondFixedSurface',ilocerr)
        CALL check_codes_err(ilocerr, 'scaledValueOfSecondFixedSurface', TRIM(yzroutine), increrr=ierr)

    CASE(PRODUCT_TYPE_OBS)

      CALL codes_set(igribclone,'deletePV',GRIB_TRUE,ilocerr)
        CALL check_codes_err(ilocerr, 'deletePV', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'productDefinitionTemplateNumber',grib%sec4%productDefinitionTemplateNumber,ilocerr)
        CALL check_codes_err(ilocerr, 'productDefinitionTemplateNumber', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'parameterCategory',grib%sec4%parameterCategory,ilocerr)
        CALL check_codes_err(ilocerr, 'parameterCategory', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'parameterNumber',grib%sec4%parameterNumber,ilocerr)
        CALL check_codes_err(ilocerr, 'parameterNumber', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'typeOfGeneratingProcess',grib%sec4%typeOfGeneratingProcess,ilocerr)
        CALL check_codes_err(ilocerr, 'typeOfGeneratingProcess', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'numberOfRadarSitesUsed',grib%sec4%numberOfRadarSitesUsed,ilocerr)
        CALL check_codes_err(ilocerr, 'numberOfRadarSitesUsed', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'indicatorOfUnitOfTimeRange',grib%sec4%indicatorOfUnitOfTimeRange,ilocerr)
        CALL check_codes_err(ilocerr, 'indicatorOfUnitOfTimeRange', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'siteLatitude',grib%sec4%siteLatitude,ilocerr)
        CALL check_codes_err(ilocerr, 'siteLatitude', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'siteLongitude',grib%sec4%siteLongitude,ilocerr)
        CALL check_codes_err(ilocerr, 'siteLongitude', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'siteElevation',grib%sec4%siteElevation,ilocerr)
        CALL check_codes_err(ilocerr, 'siteElevation', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'siteId',grib%sec4%siteId,ilocerr)
        CALL check_codes_err(ilocerr, 'siteId', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'constantAntennaElevationAngle',grib%sec4%constantAntennaElevationAngle,ilocerr)
        CALL check_codes_err(ilocerr, 'constantAntennaElevationAngle', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'rangeBinSpacing',grib%sec4%rangeBinSpacing,ilocerr)
        CALL check_codes_err(ilocerr, 'rangeBinSpacing', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'radialAngularSpacing',grib%sec4%radialAngularSpacing,ilocerr)
        CALL check_codes_err(ilocerr, 'radialAngularSpacing', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'operatingMode',grib%sec4%operatingMode,ilocerr)
        CALL check_codes_err(ilocerr, 'operatingMode', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'qualityControlIndicator',grib%sec4%qualityControlIndicator,ilocerr)
        CALL check_codes_err(ilocerr, 'qualityControlIndicator', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'clutterFilterIndicator',grib%sec4%clutterFilterIndicator,ilocerr)
        CALL check_codes_err(ilocerr, 'clutterFilterIndicator', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'reflectivityCalibrationConstant',grib%sec4%reflectivityCalibrationConstant,ilocerr)
        CALL check_codes_err(ilocerr, 'reflectivityCalibrationConstant', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'accumulationInterval',grib%sec4%accumulationInterval,ilocerr)
        CALL check_codes_err(ilocerr, 'accumulationInterval', TRIM(yzroutine), increrr=ierr)
      CALL codes_set(igribclone,'referenceReflectivityForEchoTop',grib%sec4%referenceReflectivityForEchoTop,ilocerr)
        CALL check_codes_err(ilocerr, 'referenceReflectivityForEchoTop', TRIM(yzroutine), increrr=ierr)

    END SELECT

    !----------------------------------------------------------------
    !                DATA REPRESENTATION SECTION 5
    !----------------------------------------------------------------

    CALL codes_set(igribclone,'bitsPerValue',grib%sec5%bitsPerValue,ilocerr)
      CALL check_codes_err(ilocerr, 'bitsPerValue', TRIM(yzroutine), increrr=ierr)
    CALL codes_set(igribclone,'typeOfOriginalFieldValues',grib%sec5%typeOfOriginalFieldValues,ilocerr)
      CALL check_codes_err(ilocerr, 'typeOfOriginalFieldValues', TRIM(yzroutine), increrr=ierr)

    !----------------------------------------------------------------
    !                      BIT-MAP SECTION 6
    !----------------------------------------------------------------

    CALL codes_set(igribclone,'bitmapPresent',grib%sec6%bitmapPresent,ilocerr)
      CALL check_codes_err(ilocerr, 'bitmapPresent', TRIM(yzroutine), increrr=ierr)

    IF(grib%sec6%bitmapPresent == GRIB_TRUE) THEN
      CALL codes_set(igribclone, 'missingValue', grib%sec6%bitmapValue,ilocerr)
        CALL check_codes_err(ilocerr, 'missingValue', TRIM(yzroutine), increrr=ierr)
    ENDIF

    !----------------------------------------------------------------
    !                       DATA SECTION 7
    !----------------------------------------------------------------

    CALL codes_set(igribclone,'values',grib%sec7%values,ilocerr)
      CALL check_codes_err(ilocerr, 'values', TRIM(yzroutine), increrr=ierr)

    ! Grib packing, which has to be defined after the values have been written to grib record:
    CALL codes_set(igribclone,'packingType',TRIM(grib%sec5%packingType),ilocerr)  ! = grib%sec5%dataRepresentationTemplateNumber
      CALL check_codes_err(ilocerr, 'packingType', TRIM(yzroutine), increrr=ierr)

    ! Clean up:
    DEALLOCATE( grib%sec7%values )

    !----------------------------------------------------------------

    IF (ierr /= GRIB_SUCCESS) THEN
      error = 5
      errmsg(:) = ' '
      errmsg = 'ERROR '//TRIM(yzroutine)//': codes_set failed for some keys or the values. Check earlier ERROR messages above!'
    ELSE
      CALL codes_write(igribclone,ioutfile, ilocerr)
      IF (ilocerr /= GRIB_SUCCESS) THEN
        errmsg(:) = ' '
        CALL codes_get_error_string(ilocerr, errmsg)
        errmsg = 'ERROR '//TRIM(yzroutine)//' writing grib: '//TRIM(errmsg)
        error = 6
      END IF
    END IF

    CALL codes_release(igribclone)
    CALL codes_release(igribsample)

    CALL codes_close_file(ioutfile, ilocerr)
    IF (ilocerr /= GRIB_SUCCESS) THEN
      errmsg(:) = ' '
      CALL codes_get_error_string(ilocerr, errmsg)
      errmsg = 'ERROR '//TRIM(yzroutine)//' closing '//TRIM(grib_outfilename)//': '//TRIM(errmsg)
      error = 7
    END IF

    IF (lverbose) WRITE (*,*) 'done with '//TRIM(yzroutine), ' on proc ', my_radar_id

  END SUBROUTINE write_grib_to_file

  !===============================================================================

#endif

END MODULE radar_output_volscan
