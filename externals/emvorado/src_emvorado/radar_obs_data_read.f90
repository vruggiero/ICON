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


MODULE radar_obs_data_read

!------------------------------------------------------------------------------
!
! Description:
!   This module provides methods to read the radar observations
!   for the radar forward operator EMVORADO and to print some diagnostic
!   output to some files.
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp, sp
  
#ifdef NETCDF
  USE radar_data, ONLY :  &
       miss_threshold, miss_value, obsfile_missingname, zero_value, Z_crit_radar, &
       missthr_int, missval_int, &
       radar_data_type, radar_meta_type, ndatakind, cndatakind, cmaxlen, cobsflen, &
       rs_meta, rs_data, missing_obs, miss_threshold, &
       my_radar_id, &
       i_vrad, i_qualvrad, i_dbzh, i_qualdbzh, &
       i_zdr, i_kdp, i_phidp, i_rhv, i_ldr, i_cflags, &
       i_dwd, i_meteoswiss, i_arpasim, i_belgium, i_denmark, i_france, i_poland, i_czech, i_netherlands, i_slovakia, &
       c_format_cdfin, &
       c_format_h5_native_smss, c_format_h5_native_smms, c_format_h5_native_mmss, c_format_h5_native_mmms,    & 
       c_format_h5_2_opera_smss, c_format_h5_2_opera_smms, c_format_h5_2_opera_mmss, c_format_h5_2_opera_mmms

  USE radar_data_namelist, ONLY : ldebug_radsim, &
       loutdbz, loutradwind, lqc_flag, loutpolstd, loutpolall, &
       itype_obserr_vr, ydirradarin, &
       lequal_azi_alldatasets

  USE radar_interface, ONLY : get_obstime_ind_of_currtime

  USE radar_utilities, ONLY : split_string, dbz_to_linear, linear_to_dbz, ind2sub3D, sub2ind3D, bubblesort

  USE netcdf, ONLY :  &
       nf90_open, &
       nf90_noerr, &
       NF90_NOWRITE, &
       nf90_inq_varid , &
       nf90_get_var, &
       nf90_close, &
       nf90_get_att, &
       nf90_strerror

#ifdef HDF5_RADAR_INPUT
  ! requires HDF5 version 1.8.X and higher!!!
  USE hdf5
  USE, INTRINSIC :: iso_c_binding
  USE radar_obs_meta_read, ONLY : read_attr_odim, read_dataset_odim, get_names_of_datasets_odim_h5, &
       reconstruct_file_series_opera, reconstruct_file_series_mch
#endif
#endif

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !================================================================================
  !================================================================================

  IMPLICIT NONE

  PRIVATE

  REAL(kind=dp), PARAMETER           :: eps = 1e-20_dp

  !==============================================================================
  ! Public Subroutines:

#ifdef NETCDF
  PUBLIC :: read_obs_rad
#endif

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS


#ifdef NETCDF

#ifdef HDF5_RADAR_INPUT

  SUBROUTINE read_field_obs_opera_h5(ista, idata, itime, hdf5_shortname, &
       datap, zind_intp_obs, ldata_avail)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine reads observed data of OPERA radars in hdf5 format for each
    !              radar station for each time in a rather generic fashion.
    !
    !    However, it expects the input format to be
    !    HDF5-files of the type that is sent to OPERA by each country. The files
    !    contain only one timestep. Concerning moments and sweeps, there are the
    !    following possibilities:
    !
    !    - single-moment single-sweep  "/dataset1/data1/"
    !    - multi-moment single-sweep   "/dataset1/data1:N/"
    !    - single-moment multi-sweep   "/dataset1:M/data1/"
    !    - multi-moment multi-sweep    "/dataset1:M/data1:N/"
    !
    !    IMPORTANT: This reader expects that "/dataset1:M/dataYY/"" represent M different elevations
    !               and "/datasetXX/data1:N" the N different radar moments!
    !
    !    TO BE CHECKED FOR EACH COUNTRY!!!
    !
    !    The naming convention of the files is:
    !
    !    T_PA<char1-var-ident><char1-ele-ident><int2-station-discr>_<char4-center-ident>_<YYYYMMDDhhmmss>.<ext>
    !
    !       <char1-scan-ident> : one character to identify the scantype
    !       <char1-ele-ident> : one character to identify the elevation number, typically 'A', 'B', 'C', ...
    !       <int2-station-discr> : 2-digit integer to discriminate different stations
    !       <char4-center-ident> : 4 characters to identify the NMS center which has sent the data, e.g.
    !                               'EDZW' for Offenbach (Germany), 'LSSW' for Switzerland (Swiss)
    !
    !       <YYYYMMDDhhmmss> : date and time string for the end time of the volume scan to which this file belongs
    !       <ext> : file extention, typically 'h5' or 'hdf'
    !
    !    E.g. a set from MeteoSwiss (center LSSW):
    !
    !     "T_PAGA52_C_LSSW_20220104101500.hdf"
    !     "T_PAGB52_C_LSSW_20220104101500.hdf"
    !     "T_PAGC52_C_LSSW_20220104101500.hdf"
    !     "T_PAGD52_C_LSSW_20220104101500.hdf"
    !     "T_PAGE52_C_LSSW_20220104101500.hdf"
    !
    !     "T_PAGA41_C_LSSW_20220104102500.hdf"
    !     "T_PAGB41_C_LSSW_20220104102500.hdf"
    !     "T_PAGC41_C_LSSW_20220104102500.hdf"
    !     "T_PAGD41_C_LSSW_20220104102500.hdf"
    !     "T_PAGE41_C_LSSW_20220104102500.hdf"
    !
    !
    !
    ! Method: The HDF5 files for one or all variable(s) and one entire volume scan
    !         are opened within this routine.
    !         "datap" has to be a pointer to a 1D data vector in
    !         the rs_data(ista) structure, which points to the desired
    !         input field. "ncdf_shortname" denotes the NetCDF shortname
    !         of the input field as determined from ncdump.
    !         "ldata_avail" is a flag indicating failure or success of reading.
    !         "missing_data" is the file-internal missing value for the data.
    !
    ! IS CALLED FROM "READ_OBS_RAD".
    !
    ! Caveat: At the moment, this subroutine has some side-effects through
    !         some global variables / data vectors, which are only correctly
    !         integrated if this subroutine is called from "read_obs_rad"!
    !
    !------------------------------------------------------------------------------


    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(in):

    INTEGER, INTENT(IN)    :: &
         ista     , & ! Station nr.
         idata    , & ! data type index: 1=vr, 2=qv, 3=z, 4=qz
         itime        ! time index, from 1 ... rs_meta(ista)%nobs_times

    CHARACTER (len=*), INTENT(in) :: hdf5_shortname ! HDF5 short name of radar data field

    !.. INTENT(out):
    REAL (KIND=dp), POINTER  :: datap(:) ! storage vector for the data
    INTEGER,        POINTER  :: zind_intp_obs(:) ! storage vector for the data index
    LOGICAL, INTENT(out)     :: ldata_avail   ! flag for success of reading

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER    :: yzroutine = 'read_field_obs_opera_h5'
    CHARACTER (LEN=80)              :: yzerrmsg

    INTEGER                         :: i, ii, k, m, n, o, ierr, nel, nel_obs, nobs_max, nobs
    REAL(KIND=dp), ALLOCATABLE      :: manel(:), indata3d(:,:,:)

    CHARACTER(len=cobsflen)         :: infilename
    CHARACTER(len=cobsflen), ALLOCATABLE :: infilenames(:)
    INTEGER                         :: ninfiles
    CHARACTER(len=cobsflen), ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    INTEGER                         :: nwords
    CHARACTER(len=60)               :: vardataset
    CHARACTER(len=20), ALLOCATABLE  :: dsetname(:)
    CHARACTER(len=60)               :: cpath

    INTEGER(kind=hid_t)             :: file_id

    INTEGER, ALLOCATABLE            :: rbins(:),nazis(:)
    LOGICAL, ALLOCATABLE            :: ldata_avail_dataset(:)

    LOGICAL :: is_multimoment, is_multisweep

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,'(a,i4,a,i6.6)') TRIM(yzroutine)//' on proc ', my_radar_id, ' for ', rs_meta(ista)%station_id

    ! .. Check index idata:
    IF (idata < 1 .OR. idata > ndatakind) THEN
      WRITE (*,'(a,/,a,i2,a)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           hdf5_shortname//') from an ODIM-HDF5 file:', &
           '  idatakind = ', idata, ' is not between 1 and '//cndatakind//'!'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    ! .. Return silently if there is no obs file for this quantity:
    IF (TRIM(rs_meta(ista)%obsfile(itime,idata)) == obsfile_missingname) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    !=====================================================================================================================
    ! Re-construct the filename series of the PPI-files from the given input filename,
    !  which contains the basic filename, the number of PPI sweeps and the list of elevation identifiers:
    !   (This should work also for ninfiles = 1)
    !=====================================================================================================================


    ! 1) Check for multi-/single-moment or multi-/single-sweep infile:
    SELECT CASE (rs_meta(ista)%obsfile_format(itime))
    CASE (c_format_h5_2_opera_mmms)
      is_multimoment = .TRUE.
      is_multisweep  = .TRUE.
    CASE (c_format_h5_2_opera_smms)
      is_multimoment = .FALSE.
      is_multisweep  = .TRUE.
    CASE (c_format_h5_2_opera_mmss)
      is_multimoment = .TRUE.
      is_multisweep  = .FALSE.
    CASE (c_format_h5_2_opera_smss)
      is_multimoment = .FALSE.
      is_multisweep  = .FALSE.
    CASE (c_format_h5_native_mmms)
      is_multimoment = .TRUE.
      is_multisweep  = .TRUE.
    CASE (c_format_h5_native_smms)
      is_multimoment = .FALSE.
      is_multisweep  = .TRUE.
    CASE (c_format_h5_native_mmss)
      is_multimoment = .TRUE.
      is_multisweep  = .FALSE.
    CASE (c_format_h5_native_smss)
      is_multimoment = .FALSE.
      is_multisweep  = .FALSE.
    CASE default
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': File format '// TRIM(rs_meta(ista)%obsfile_format(itime)) // &
           ' not supported for reading: ' // TRIM(ydirradarin)//TRIM(rs_meta(ista)%obsfile(itime,idata))
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN      
    END SELECT

    ! 2) Determine if it is a multi- or single-moment file, and define the vardataset to read from: 
    infilename(:) = ' '
    CALL split_string (TRIM(rs_meta(ista)%obsfile(itime,idata)), '/', LEN(words), words, nwords)
    IF (nwords < 1 .OR. nwords > 2) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': Strange filename problem in '// &
           TRIM(ydirradarin)//TRIM(rs_meta(ista)%obsfile(itime,idata))
      IF (ALLOCATED(words)) DEALLOCATE(words)
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    infilename = TRIM(words(1))
    vardataset(:) = ' '
    IF(nwords == 2) THEN
      IF (is_multimoment) THEN
        vardataset = TRIM(words(2))
      ELSE
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': File is not multi-moment as expected from the name: '// &
             TRIM(ydirradarin)//TRIM(infilename)
        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF
    ELSE
      IF (.NOT. is_multimoment) THEN
        vardataset = 'data1'
      ELSE
        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': File is not single-moment as expected from the name: '// &
             TRIM(ydirradarin)//TRIM(infilename)
        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF
    END IF

    ! 3) Reconstruct the file name series, if we have single-sweep files:

    IF ( rs_meta(ista)%icountry == i_meteoswiss .AND. &
         ( TRIM(rs_meta(ista)%obsfile_format(itime)) == TRIM(c_format_h5_native_smss) .OR. &
           TRIM(rs_meta(ista)%obsfile_format(itime)) == TRIM(c_format_h5_native_mmss) ) ) THEN

      CALL reconstruct_file_series_mch (yzroutine, infilename, infilenames, ninfiles, ierr)
      IF (ierr /= 0) THEN
        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF

    ELSE

      CALL reconstruct_file_series_opera (yzroutine, infilename, infilenames, ninfiles, ierr)
      IF (ierr /= 0) THEN
        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF

    END IF
      
    ! Check consistency of the format specifier and the actual number of files in the series:
    IF (is_multisweep .AND. ninfiles > 1) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': File format spec is multi-sweep, but filename suggests single-sweeps: '// &
           TRIM(ydirradarin)//TRIM(infilename)
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    !=====================================================================================================================
    ! Initialize HDF5 FORTRAN interface:
    !=====================================================================================================================

    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(infilenames)
      RETURN
    ENDIF

    !=====================================================================================================================
    ! .. Now all filenames of the respective timestep / volume scan are known, so read their data:
    !=====================================================================================================================

    !    Define three dimensional start and count arrays
    nel      = rs_meta(ista)%nel
    !    Define max. number of possible observations:
    nobs_max = rs_meta(ista)%nra_obs * rs_meta(ista)%naz * rs_meta(ista)%nel

    !    Allocate temporary storage array for the data fields from the HDF5-files:
    ALLOCATE( indata3d(rs_meta(ista)%nra_obs, rs_meta(ista)%naz, nel) )
    indata3d = miss_value

    ! .. Set %check - fields for checking of equivalence of
    !    available azimuts and elevations for different data sets
    !    from the same station from the same time:
    ALLOCATE( manel(nel) , &
              nazis(nel) , &
              rbins(nel)   )

    manel(:) = missing_obs
    nazis(:) = 0
    rbins(:) = 0


    !=====================================================================================================================
    ! Read data in a loop over the infilenames from path "/datasetXX/dataYY":
    !  If there is only one file, it is:
    !   - either a volume scan with a single elevation in "/dataset1" and vardataset contains the correct dataYY string.
    !   - or it is a multi-sweep file where we have to read the single elevations from the different "/datasetXX".
    !  If it is a series of single-sweep files, there is one elevation per file in "/dataset1" and 
    !   vardataset contains the correct dataYY string.
    !=====================================================================================================================

    IF (is_multisweep) THEN

      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(1)), H5F_ACC_RDONLY_F, file_id, ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5 file '//TRIM(ydirradarin)//TRIM(infilenames(1))
        ldata_avail = .FALSE.
      ENDIF
        
      ! .. For file "file_id", detect the number of "dataset"s in the root group (nel_obs)
      !    and allocate/build the list of dataset names (dsetname).
      !    This will allocate dsetname. yzroutine and infilename are only for annotation of error messages:
      CALL get_names_of_datasets_odim_h5(file_id, TRIM(yzroutine), TRIM(infilename), '/', 'dataset', &
           LEN(dsetname), nel_obs, dsetname, ierr)

      ALLOCATE( ldata_avail_dataset(nel_obs) ) 
      ldata_avail_dataset(:) = .TRUE.

      ele_loop_ms: DO ii=1, nel_obs

        CALL get_dataset_odim_h5 (TRIM(yzroutine), TRIM(ydirradarin)//TRIM(infilenames(1)), file_id, hdf5_shortname, &
             &                    idata, TRIM(dsetname(ii)), TRIM(vardataset), rs_meta(ista), ldata_avail_dataset(ii), &
             &                    indata3d, rbins, nazis, manel, &
             &                    lcheck_azi_sorting=.FALSE., lwarn_on_missing_startazA=.FALSE., lele_allowed_to_be_wrong=.TRUE.)

      END DO ele_loop_ms

      DEALLOCATE(dsetname)

      CALL h5fclose_f(file_id, ierr)

    ELSE

      ALLOCATE( ldata_avail_dataset(ninfiles) ) 
      ldata_avail_dataset(:) = .TRUE.

      ele_loop_ss: DO ii=1, ninfiles

        CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(ii)), H5F_ACC_RDONLY_F, file_id, ierr)
        IF (ierr /= 0) THEN
          WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5 file '//TRIM(ydirradarin)//TRIM(infilenames(ii))
          ldata_avail = .FALSE.
        ENDIF
        
        CALL get_dataset_odim_h5 (TRIM(yzroutine), TRIM(ydirradarin)//TRIM(infilenames(ii)), file_id, hdf5_shortname, &
             &                    idata, 'dataset1', TRIM(vardataset), rs_meta(ista), ldata_avail_dataset(ii), &
             &                    indata3d, rbins, nazis, manel, &
             &                    lcheck_azi_sorting=.FALSE., lwarn_on_missing_startazA=.FALSE., lele_allowed_to_be_wrong=.TRUE.)
        
        CALL h5fclose_f(file_id, ierr)

      END DO ele_loop_ss
      
    END IF

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. If any of the elevations could be sucessfully read, we set ldata_avail=.true.:
    ldata_avail = ANY(ldata_avail_dataset)
    
    IF (.NOT. ldata_avail) THEN
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,infilenames,ldata_avail_dataset)
      RETURN
    END IF

    CALL check_and_complete_dataset ( TRIM(yzroutine), rs_data(ista), rs_meta(ista), &
         &                            itime, idata, nel, nobs_max, manel, nazis, rbins, &
         &                            zind_intp_obs, ierr)
    IF (ierr /= 0) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,infilenames,ldata_avail_dataset)
      RETURN
    END IF
    

    DEALLOCATE(manel,rbins,nazis,infilenames,ldata_avail_dataset)

    ! .. Allocate memory for the 1D output vector:
    IF (ASSOCIATED(datap)) DEALLOCATE(datap)
    NULLIFY(datap)
    ALLOCATE(datap(nobs_max))

    ! .. Put the data into the 1D output vector:
    nobs = 0
    DO o = 1, rs_meta(ista)%nel        ! Loop over Reports (one elevation per record)
      DO m = 1, rs_meta(ista)%naz      ! Loop over Azimut
        DO n = 1,rs_meta(ista)%nra_obs ! Loop over Range
          nobs = nobs + 1
          datap(nobs) = indata3d(n,m,o)
        END DO
      END DO
    END DO

    DEALLOCATE(indata3d)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_field_obs_opera_h5

  !==============================================================================

  SUBROUTINE read_field_obs_dwd_h5(ista, idata, itime, hdf5_shortname, &
       datap, zind_intp_obs, ldata_avail)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine reads observed data of DWD radars in hdf5 format for each
    !              radar station for each time in a rather generic fashion.
    !              However, it expects the input format to be ODIM-HDF5 (dialect of DWD)
    !              and that there is one file per PPI-elevation per station per
    !              time, i.e., the way that DWD or MeteoSwiss store their radar obs.
    !
    ! Method: The HDF5 files for one or all variable(s) and one volume scan
    !         are opened within this routine.
    !         "datap" has to be a pointer to a 1D data vector in
    !         the rs_data(ista) structure, which points to the desired
    !         input field. "ncdf_shortname" denotes the NetCDF shortname
    !         of the input field as determined from ncdump.
    !         "ldata_avail" is a flag indicating failure or success of reading.
    !         "missing_data" is the file-internal missing value for the data.
    !
    ! IS CALLED FROM "READ_OBS_RAD".
    !
    ! Caveat: At the moment, this subroutine has some side-effects through
    !         some global variables / data vectors, which are only correctly
    !         integrated if this subroutine is called from "read_obs_rad"!
    !
    !------------------------------------------------------------------------------


    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(in):

    INTEGER, INTENT(IN)    :: &
         ista     , & ! Station nr.
         idata    , & ! data type index: 1=vr, 2=qv, 3=z, 4=qz
         itime        ! time index, from 1 ... rs_meta(ista)%nobs_times

    CHARACTER (len=*), INTENT(in) :: hdf5_shortname ! HDF5 short name of radar data field

    !.. INTENT(out):
    REAL (KIND=dp), POINTER  :: datap(:) ! storage vector for the data
    INTEGER,        POINTER  :: zind_intp_obs(:) ! storage vector for the data index
    LOGICAL, INTENT(out)     :: ldata_avail   ! flag for success of reading

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER    :: yzroutine = 'read_field_obs_dwd_h5'
    CHARACTER (LEN=80)              :: yzerrmsg

    INTEGER                         :: ii, k, kk, m, n, o, ierr, nel, nobs_max, nobs
    REAL(KIND=dp), ALLOCATABLE      :: manel(:), indata3d(:,:,:)

    CHARACTER(len=cobsflen)         :: infilename
    CHARACTER(len=cobsflen), ALLOCATABLE :: infilenames(:)
    INTEGER                         :: ninfiles
    CHARACTER(len=cobsflen), ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    INTEGER                         :: nwords
    CHARACTER(len=60)               :: vardataset

    INTEGER(kind=hid_t)             :: file_id

    INTEGER, ALLOCATABLE            :: rbins(:),nazis(:)
    LOGICAL, ALLOCATABLE            :: ldata_avail_dataset(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! .. Check index idata:
    IF (idata < 1 .OR. idata > ndatakind) THEN
      WRITE (*,'(a,/,a,i2,a)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           hdf5_shortname//') from an ODIM-HDF5 file:', &
           '  idatakind = ', idata, ' is not between 1 and '//cndatakind//'!'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    ! .. Return silently if there is no obs file for this quantity:
    IF (TRIM(rs_meta(ista)%obsfile(itime,idata)) == obsfile_missingname) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    ! .. Re-construct the filename series of the PPI-files from the given input filename,
    !     which contains only the basic filename and the number of PPI sweeps:

    ! 1) Check for multi-moment or single-moment infile:
    infilename(:) = ' '
    CALL split_string (TRIM(rs_meta(ista)%obsfile(itime,idata)), '/', LEN(words), words, nwords)
    IF (nwords < 1 .OR. nwords > 2) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Strange filename problem in '// &
           TRIM(ydirradarin)//TRIM(rs_meta(ista)%obsfile(itime,idata))
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    infilename = TRIM(words(1))
    vardataset(:) = ' '
    IF(nwords == 2) THEN
      IF (INDEX(TRIM(infilename),'sweeph5allm') > 0) THEN
        vardataset = TRIM(words(2))
      ELSE
        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': File is not sweeph5allm: '//TRIM(ydirradarin)//TRIM(infilename)
        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF
    ELSE
      IF (INDEX(TRIM(infilename),'sweeph5onem') > 0) THEN
        vardataset = 'data1'
      ELSE
        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': File is not sweeph5onem: '//TRIM(ydirradarin)//TRIM(infilename)
        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF
    END IF

    ! 1) Check for precipitation- or volume scan infile(s) and separate the elevation-date strings between "|":
    k = INDEX(infilename,'_',back=.TRUE.)   ! position of the last '_' in the string
    IF (infilename(7:13) == 'vol5min') THEN
      ! .. This is a volume scan distributed over several files:
      READ (infilename(k+1:k+2), '(i2.2)') ninfiles
      kk = INDEX(infilename,'|') - 1   ! length of the actual filename (anything before the first '|'
      CALL split_string (TRIM(infilename), '|', LEN(words), words, nwords)
      IF (nwords /= ninfiles+1) THEN
        WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Inconsistent infilename '//TRIM(infilename)
        ldata_avail = .FALSE.
        IF (ALLOCATED(words))  DEALLOCATE(words)
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        RETURN
      END IF
    ELSE IF (infilename(7:9) == 'pcp') THEN
      ! .. This is a precipitation scan in one file:
      ninfiles = 1
      kk       = LEN_TRIM(infilename)
    END IF
    ALLOCATE(infilenames(ninfiles))
    DO ii=1, ninfiles
      infilenames(ii) = infilename(1:kk)
      IF (infilename(7:13) == 'vol5min') THEN
        infilenames(ii)(k+1:k+2)   = words(ii+1)(1:2)
        infilenames(ii)(k+14:k+17) = words(ii+1)(4:7)
      END IF
    END DO
    IF (ALLOCATED(words)) DEALLOCATE(words)


    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(infilenames)
      RETURN
    ENDIF

    !=====================================================================================================================
    ! .. Now all filenames of the respective timestep / volume scan are known, so read their data:
    !=====================================================================================================================

    !    Define three dimensional start and count arrays
    nel      = rs_meta(ista)%nel
    !    Define max. number of possible observations:
    nobs_max = rs_meta(ista)%nra_obs * rs_meta(ista)%naz * rs_meta(ista)%nel

    !    Allocate temporary storage array for the data fields from the HDF5-files:
    ALLOCATE( indata3d(rs_meta(ista)%nra_obs, rs_meta(ista)%naz, nel) )
    indata3d = miss_value

    ! .. Set %check - fields for checking of equivalence of
    !    available azimuts and elevations for different data sets
    !    from the same station from the same time:
    ALLOCATE( manel(nel) , &
              nazis(nel) , &
              rbins(nel) , &
              ldata_avail_dataset(ninfiles) )

    manel(:) = missing_obs
    nazis(:) = 0
    rbins(:) = 0
    ldata_avail_dataset(:) = .TRUE.

    read_loop: DO ii=1, ninfiles

      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(ii)), H5F_ACC_RDONLY_F, file_id, ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5 file '//TRIM(ydirradarin)//TRIM(infilenames(ii))
        ldata_avail = .FALSE.
      ENDIF

      CALL get_dataset_odim_h5 (TRIM(yzroutine), TRIM(ydirradarin)//TRIM(infilenames(ii)), file_id, hdf5_shortname, &
           &                    idata, 'dataset1', TRIM(vardataset), rs_meta(ista), ldata_avail_dataset(ii), &
           &                    indata3d, rbins, nazis, manel, &
           &                    lcheck_azi_sorting=.TRUE., lwarn_on_missing_startazA=.TRUE., lele_allowed_to_be_wrong=.FALSE.)


      CALL h5fclose_f(file_id, ierr)

    END DO read_loop

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. If any of the elevations could be sucessfully read, we set ldata_avail=.true.:
    ldata_avail = ANY(ldata_avail_dataset)

    IF (.NOT. ldata_avail) THEN
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,infilenames,ldata_avail_dataset)
      RETURN
    END IF

    CALL check_and_complete_dataset ( TRIM(yzroutine), rs_data(ista), rs_meta(ista), &
         &                            itime, idata, nel, nobs_max, manel, nazis, rbins, &
         &                            zind_intp_obs, ierr)
    IF (ierr /= 0) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,infilenames,ldata_avail_dataset)
      RETURN
    END IF
    

    DEALLOCATE(manel,rbins,nazis,infilenames,ldata_avail_dataset)

    ! .. Allocate memory for the 1D output vector:
    IF (ASSOCIATED(datap)) DEALLOCATE(datap)
    NULLIFY(datap)
    ALLOCATE(datap(nobs_max))

    ! .. Put the data into the 1D output vector:
    nobs = 0
    DO o = 1, rs_meta(ista)%nel        ! Loop over Reports (one elevation per record)
      DO m = 1, rs_meta(ista)%naz      ! Loop over Azimut
        DO n = 1,rs_meta(ista)%nra_obs ! Loop over Range
          nobs = nobs + 1
          datap(nobs) = indata3d(n,m,o)
        END DO
      END DO
    END DO

    DEALLOCATE(indata3d)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_field_obs_dwd_h5

#endif

  !==============================================================================

  !==============================================================================
  !+ Module procedure in radar_src for the reading of radar data from DWD CDFIN
  !------------------------------------------------------------------------------

  SUBROUTINE read_field_obs_dwd(ista, idata, itime, ncdf_shortname, &
       datap, zind_intp_obs, ldata_avail, missing_data)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine reads observed data for each
    !              radar station for each time in a generic fashion.
    !              However, it expects the data to be in DWD's NetCDF
    !              format that is produced by "burfx2netcdf" software.
    !
    ! Method: The NetCDF file id "fid" has to point to an open file.
    !         "datap" has to be a pointer to a 1D data vector in
    !         the rs_data(ista) structure, which points to the desired
    !         input field. "ncdf_shortname" denotes the NetCDF shortname
    !         of the input field as determined from ncdump.
    !         "ldata_avail" is a flag indicating failure or success of reading.
    !         "missing_data" is the correct missing value for the data.
    !
    ! IS CALLED FROM "READ_OBS_RAD".
    !
    ! Caveat: At the moment, this subroutine has some side-effects through
    !         some global variables / data vectors, which are only correctly
    !         integrated if this subroutine is called from "read_obs_rad"!
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(in):

    INTEGER, INTENT(IN)    :: &
         ista     , & ! Station nr.
         idata    , & ! data type index: 1=vr, 2=qv, 3=z, 4=qz
         itime        ! time index, from 1 ... rs_meta(ista)%nobs_times

    CHARACTER (len=*), INTENT(in) :: ncdf_shortname ! NetCDF short name of radar data field

    !.. INTENT(out):
    REAL (KIND=dp), POINTER      :: datap(:) ! storage vector for the data
    INTEGER       , POINTER      :: zind_intp_obs(:) ! storage vector for the data index
    LOGICAL, INTENT(out)         :: ldata_avail   ! flag for success of reading
    REAL (KIND=dp), INTENT(out)  :: missing_data  ! correct missing value for the data

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'read_field_obs_dwd'
    CHARACTER (LEN=80)              :: yzerrmsg
    CHARACTER (len=500)             :: ncdffile

    INTEGER                         :: fid          ! file ID

    ! NetCDF variable id for
    INTEGER                         :: varid_mabaz  ! azimuth mabaz
    INTEGER                         :: varid_mabaz0 ! azimuth mabaz0
    INTEGER                         :: varid_manel  ! elevation manel
    INTEGER                         :: varid_manel0 ! elevation manel
    INTEGER                         :: varid_data   ! radar field shortname
    INTEGER                         :: varid_time   ! observation time
    INTEGER                         :: varid_rbins  ! number of range bins of individual report

    INTEGER                         :: m,n,o,status,status1,status2,AziDim
    REAL    (KIND=dp)               :: mdmvr, missing_az
    REAL    (KIND=dp), ALLOCATABLE  :: manel0(:,:),mabaz0(:,:),manel(:),mabaz(:)

    INTEGER          , ALLOCATABLE  :: time(:),rbins(:),azind(:,:),eleind(:)
    INTEGER                         :: nobs
    INTEGER                         :: ierr,nsta(3),ncount(3),nrep, nobs_max, &
                                       start_rec, end_rec, tol, zcheck(4)
    CHARACTER (len=cmaxlen)         :: info_in_case_of_error

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !    Check index idata:
    IF (idata < 1 .OR. idata > ndatakind) THEN
      WRITE (*,'(a,/,a,i2,a)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           ncdf_shortname//') from a NetCDF file:', &
           '  idatakind = ', idata, ' is not between 1 and '//cndatakind//'!'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    ! .. Return silently if there is no obs file for this quantity:
    IF (TRIM(rs_meta(ista)%obsfile(itime,idata)) == obsfile_missingname) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    !    Get the NetCDF start- and end record numbers from the meta-data-structure for
    !    the data type with index idata (1,2,3, ...):
    start_rec = rs_meta(ista)%obs_startrec(itime,idata)
    end_rec   = rs_meta(ista)%obs_endrec(itime,idata)

    !    Define three dimensional start and count arrays
    nrep      = end_rec - start_rec + 1
    nsta(1:3) = (/1,1,start_rec/)
    ncount(1:3) = (/rs_meta(ista)%nra_obs,rs_meta(ista)%naz_ncdf(itime,idata),nrep/)


    !    Define max. number of possible observations:
    nobs_max = rs_meta(ista)%naz_ncdf(itime,idata) * MAX(rs_meta(ista)%nel,nrep) * rs_meta(ista)%nra_obs

    ! .. NOTE: if one or several elevations are doubled, we read the doubled data, compute their
    !          correct address (ind_intp) in the (ra,az,el)-space and continue as normal. In the course
    !          of the emvorado computations, data will be simulated twice for doubled elevations, but in the step
    !          where we sort the data from the 1-D vector into the 3-D (ra,az,el) arrays for output,
    !          the last of the data chunks for a doubled elevation will overwrite previous data chunks,
    !          so that the last record of observations will make it into the output, while earlier
    !          records will be forgotten.

    ! .. Set %check - fields for checking of equivalence of
    !    available azimuts and elevations for different data sets
    !    from the same station from the same time:

    IF (ALLOCATED(time))    DEALLOCATE(time,manel,rbins,mabaz)
    ALLOCATE( time (start_rec:end_rec) , &
              manel(start_rec:end_rec) , &
              mabaz(start_rec:end_rec) , &
              rbins(start_rec:end_rec) )
    time(:) = -99
    manel(:) = missing_obs
    mabaz(:) = missing_obs
    rbins(:) = 0

    ncdffile(:) = ' '
    ncdffile = TRIM(ydirradarin)//TRIM(rs_meta(ista)%obsfile(itime,idata))

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//': opening '//TRIM(ncdffile)//' ...'
    status = nf90_open (TRIM(ncdffile), NF90_NOWRITE, fid)
    IF (status /= NF90_NOERR) THEN
      WRITE(*,'(i5,1x,3a)') my_radar_id, 'WARNING '//TRIM(yzroutine)//': problem opening '//' NetCDF file ', TRIM(ncdffile)//' ', &
           TRIM(nf90_strerror(status))
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    status = nf90_inq_varid (fid, 'TIME'  , varid_time)          ! Reference time
    status = nf90_get_var (fid, varid_time, time, (/start_rec/), (/nrep/))
    status = nf90_inq_varid (fid, 'MANEL' , varid_manel)         ! Elevation nominal
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (fid, 'ppi_elevation',varid_manel)
    status = nf90_get_var (fid, varid_manel,manel,(/start_rec/), (/nrep/))
    status = nf90_inq_varid (fid, 'MABAZ' , varid_mabaz)         ! Azimuth start of each record (end of sampling interval)
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (fid, 'ppi_azimuth',varid_mabaz)
    status = nf90_get_var (fid, varid_mabaz,mabaz,(/start_rec/), (/nrep/))
    status = nf90_inq_varid (fid, '030194', varid_rbins)
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (fid, 'n_range_bins',varid_rbins)
    status = nf90_get_var (fid, varid_rbins, rbins, (/start_rec/), (/nrep/))

    IF (.NOT.lequal_azi_alldatasets .OR. MAXVAL(rs_data(ista)%check) == -99) THEN
      ! .. Either the %check - fields have not been initialized so far, i.e., first call
      !    for the actual time, or lequal_azi_alldatasets = .false., i.e., which means
      !    that the following is called for every dataset and we do not assume
      !    that azimuts are exactly equal for all datasets.

      IF (ALLOCATED(mabaz0))           DEALLOCATE(mabaz0,manel0,azind,eleind)
      ALLOCATE( mabaz0(rs_meta(ista)%naz_ncdf(itime,idata),start_rec:end_rec) , &
           manel0(rs_meta(ista)%naz_ncdf(itime,idata),start_rec:end_rec) , &
           azind (rs_meta(ista)%naz_ncdf(itime,idata),start_rec:end_rec) , &
           eleind(start_rec:end_rec)                           )

      status = nf90_inq_varid (fid, 'MABAZ0', varid_mabaz0)        ! Azimuth (end of sampling interval)
        IF (status /= NF90_NOERR ) status = nf90_inq_varid (fid, 'ray_azimuth',varid_mabaz0)
      status = nf90_get_att (fid, varid_mabaz0, "_FillValue", missing_az)
      mabaz0 = missing_az
      status = nf90_get_var (fid, varid_mabaz0, mabaz0, (/nsta(2),start_rec/), (/ncount(2),nrep/))
      IF ( TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
        ! .. Exception for DWD precipitation scan: to bypass elevation checking, preset
        !     manel0 with the elevation-constant nominal dummy values of the precip scan:
        DO o = start_rec,end_rec
          manel0(:,o) = manel(o)
        END DO
      ELSE
        ! .. If not precip scan, read the true elevations for each ray:
        status = nf90_inq_varid (fid, 'MANEL0', varid_manel0)       ! Elevation
          IF (status /= NF90_NOERR ) status = nf90_inq_varid (fid, 'ray_elevation',varid_manel0)
        status = nf90_get_var (fid, varid_manel0, manel0, (/nsta(2),start_rec/), (/ncount(2),nrep/))
      END IF
      rs_data(ista)%check(1) = SUM(time(start_rec:end_rec))
      rs_data(ista)%check(2) = SUM(NINT(manel(start_rec:end_rec)*100.))
      rs_data(ista)%check(3) = SUM(NINT(mabaz(start_rec:end_rec)*100.))
      rs_data(ista)%check(4) = SUM(rbins(start_rec:end_rec))
      ! .. Get azimut- and elevation- indices for sorting of the radar data into a regular
      !    3D data set with nominal range, azimut and elevation:

      ! .. Azimut: assign azimut index by rounding to the next regular azimut value
      ! .. Elevation: assign index from next regular elevation in rs_meta(ista)%el_arr(:)
      ! .. mabaz0 is the end of the azimut interval in sweep direction, so, in the following
      !    subroutine, add/subtract half of az_inc (dep. on antenna rotation direction) to
      !    center it on the center of the azimut interval

      info_in_case_of_error (:) = ' '
      WRITE (info_in_case_of_error, '(a,1x,i6.6,1x,a,1x,a,1x,a)') &
           TRIM(rs_meta(ista)%station_name), &
           rs_meta(ista)%station_id, &
           TRIM(ydirradarin)//TRIM(rs_meta(ista)%obsfile(1,3)), &
           TRIM(rs_meta(ista)%obs_cdate(itime)), &
           TRIM(ncdf_shortname)

      CALL get_posind_rounded_fullrev (TRIM(info_in_case_of_error), &
           mabaz0(:,:), manel0(:,:), &
           rs_meta(ista)%naz_ncdf(itime,idata), nrep, missing_az, &
           rs_meta(ista)%az_start, rs_meta(ista)%az_inc, &
           rs_meta(ista)%el_arr  , rs_meta(ista)%nel,    &
           azind(:,:), eleind(:) , rs_meta(ista)%naz)

      ! .. Allocate and calculate the vector for storing the running index from
      !    which the indices for range, azimut and elevation can be restored:
      IF (ASSOCIATED(zind_intp_obs)) DEALLOCATE(zind_intp_obs)
      NULLIFY(zind_intp_obs)

      IF (lequal_azi_alldatasets) THEN
        IF (ASSOCIATED(rs_data(ista)%ind_intp_obs)) DEALLOCATE(rs_data(ista)%ind_intp_obs)
        NULLIFY(rs_data(ista)%ind_intp_obs)
        ALLOCATE(rs_data(ista)%ind_intp_obs(nobs_max))
        rs_data(ista)%ind_intp_obs = missval_int
        zind_intp_obs => rs_data(ista)%ind_intp_obs
      ELSE
        ALLOCATE(zind_intp_obs(nobs_max))
        zind_intp_obs = missval_int
      END IF

      ! .. 1D running index relative to the nominal ranges, azimuths and elevations as given by
      !    rs_meta(ista)%%az_start, %az_inc, %naz, %nra_obs, and %el_arr): 
      nobs  = 0
      DO o = start_rec,end_rec           ! Loop over Reports (one elevation per record)
        DO m = 1,rs_meta(ista)%naz_ncdf(itime,idata)  ! Loop over Azimut
!CDIR NODEP
          DO n = 1,rs_meta(ista)%nra_obs     ! Loop over Range
            nobs = nobs + 1
            IF ( azind(m,o) > missthr_int .AND. eleind(o) > missthr_int .AND. n <= rbins(o) ) THEN
              zind_intp_obs(nobs) = azind(m,o) + (n-1)*rs_meta(ista)%naz &
                   + (eleind(o)-1)*rs_meta(ista)%nra_obs*rs_meta(ista)%naz
            ELSE
!!$              rs_data(ista)%ind_intp_obs(nobs) = -9999
              zind_intp_obs(nobs) = -9999
            END IF
          END DO
        END DO
        ! .. Keep track of this elevation and add it to the list of present elevations:
        IF (rs_meta(ista)%nel_present > 0) THEN
          IF (ALL(rs_meta(ista)%ind_ele_present(1:rs_meta(ista)%nel_present) /= eleind(o) .AND. eleind(o) > missthr_int)) THEN
            rs_meta(ista)%nel_present = rs_meta(ista)%nel_present + 1
            rs_meta(ista)%ind_ele_present(rs_meta(ista)%nel_present) = eleind(o)
          END IF
        ELSE
          IF (eleind(o) > missthr_int) THEN
            rs_meta(ista)%nel_present = 1
            rs_meta(ista)%ind_ele_present(rs_meta(ista)%nel_present) = eleind(o)
          END IF
        END IF
      END DO

      IF (lequal_azi_alldatasets) THEN
!!$ UB, Alberto: set nobs only if there is a valid data file:
        DO n=1, ndatakind
          IF (rs_meta(ista)%obsfile(itime,n) /= obsfile_missingname) THEN
            rs_data(ista)%nobs_obs(n)     = nobs
          ELSE
            rs_data(ista)%nobs_obs(n)     = 0
          END IF
        END DO
      ELSE
        rs_data(ista)%nobs_obs(idata)     = nobs
      END IF

      DEALLOCATE(mabaz0,manel0,azind,eleind)

    ELSE
      ! .. Follow up call for an additional data set, therefore compare
      !    azimuts and elevations of this data set with %check:

      ! .. Some azimut- and/or elevation values in different data type files might deviate
      !    by 0.01 degree (rounding error???), so we need a little tolerance
      !    in elevation- and azimut checking:

      tol = nrep*1  ! tolerance in the sum of azimuts/elevations = nrep * 0.01 (nrep*1 because of factor 100 below!)

      zcheck(1) = SUM(time(start_rec:end_rec))
      zcheck(2) = SUM(NINT(manel(start_rec:end_rec)*100.))
      zcheck(3) = SUM(NINT(mabaz(start_rec:end_rec)*100.))
      zcheck(4) = SUM(rbins(start_rec:end_rec))

      IF ( rs_data(ista)%check(1)   /= zcheck(1)         .OR. &
           ABS( rs_data(ista)%check(2)-zcheck(2) ) > tol .OR. &
           ABS( rs_data(ista)%check(3)-zcheck(3) ) > tol .OR. &
           rs_data(ista)%check(4)   /= zcheck(4)              &
           ) THEN

        WRITE(*,'(a,i4,i8,4(/,a,2i12))') &
             TRIM(yzroutine)//': Error: Data set inconsistent '// &
             'from same radar and time: station: ',ista,rs_meta(ista)%station_id, &
             '  time checksum: ', zcheck(1), rs_data(ista)%check(1) , &
             '  elev checksum: ', zcheck(2), rs_data(ista)%check(2) , &
             '  azim checksum: ', zcheck(3), rs_data(ista)%check(3) , &
             '  rbin checksum: ', zcheck(4), rs_data(ista)%check(4)

        ldata_avail = .FALSE.
        IF (ASSOCIATED(datap)) DEALLOCATE(datap)
        NULLIFY(datap)
        DEALLOCATE(time,manel,rbins,mabaz)
        status = nf90_close(fid)
        IF (status /= NF90_NOERR) WRITE (*,'(a,/,a)') 'Error in closing file '//TRIM(ncdffile)//': '// &
                                                       TRIM(nf90_strerror(status))
        RETURN

      END IF

      IF (lequal_azi_alldatasets) THEN
        IF (ASSOCIATED(zind_intp_obs)) DEALLOCATE(zind_intp_obs)
        NULLIFY(zind_intp_obs)
        zind_intp_obs => rs_data(ista)%ind_intp_obs
      END IF

    END IF
    DEALLOCATE(time,manel,rbins,mabaz)

    ! allocate memory for the data:
    IF (ASSOCIATED(datap)) DEALLOCATE(datap)
    NULLIFY(datap)
    ALLOCATE(datap(nobs_max))

    ! get the data:   (with name exceptions for OSSEs)
    status1 = nf90_inq_varid (fid, ncdf_shortname, varid_data)
    IF (status1 /= NF90_NOERR) THEN
      SELECT CASE (TRIM(ncdf_shortname))
      CASE ('MDMVR')
        status1 = nf90_inq_varid (fid, 'radialwind', varid_data)
        IF (status1 /= NF90_NOERR) status1 = nf90_inq_varid (fid, 'vrsim', varid_data)
        IF (status1 /= NF90_NOERR) status1 = nf90_inq_varid (fid, 'vrobs', varid_data)
      CASE ('MHORRE0')
        status1 = nf90_inq_varid (fid, 'reflectivity', varid_data)
        IF (status1 /= NF90_NOERR) status1 = nf90_inq_varid (fid, 'zrsim', varid_data)
        IF (status1 /= NF90_NOERR) status1 = nf90_inq_varid (fid, 'zrobs', varid_data)
      END SELECT
    END IF
    status2 = nf90_get_att   (fid, varid_data, "_FillValue", missing_data)
    datap   = missing_data
    status  = nf90_get_var   (fid, varid_data, datap, nsta, ncount)

    IF (status == NF90_NOERR) THEN
      ldata_avail = .TRUE.
    ELSE
      WRITE (*,'(2a,6i6,/,2i12,3i6)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           ncdf_shortname//') from a NetCDF file! ', TRIM(nf90_strerror(status)), &
           fid, varid_data, status, status1, status2, nf90_noerr,  &
           start_rec, end_rec, rs_meta(ista)%nra_obs, rs_meta(ista)%naz_ncdf(itime,idata), nrep
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
    END IF

    status = nf90_close(fid)
    IF (status /= NF90_NOERR) WRITE (*,'(a,/,a)') 'Error in closing file '//TRIM(ncdffile)//': '// &
                                                  TRIM(nf90_strerror(status))

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_field_obs_dwd

  !==============================================================================

#ifdef HDF5_RADAR_INPUT

  SUBROUTINE read_field_obs_mch(ista, idata, itime, hdf5_shortname, &
       datap, zind_intp_obs, ldata_avail)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine reads observed data for each
    !              radar station for each time in a generic fashion.
    !              However, it expects the input format to be ODIM-HDF5
    !              and that there is one file per PPI-elevation per station per
    !              time, i.e., the way that MeteoSwiss stores its radar obs.
    !
    ! Method: The HDF5 files for one variable and one volume scan
    !         are opened within this routine.
    !         "datap" has to be a pointer to a 1D data vector in
    !         the rs_data(ista) structure, which points to the desired
    !         input field. "ncdf_shortname" denotes the NetCDF shortname
    !         of the input field as determined from ncdump.
    !         "ldata_avail" is a flag indicating failure or success of reading.
    !         "missing_data" is the file-internal missing value for the data.
    !
    ! IS CALLED FROM "READ_OBS_RAD".
    !
    ! Caveat: At the moment, this subroutine has some side-effects through
    !         some global variables / data vectors, which are only correctly
    !         integrated if this subroutine is called from "read_obs_rad"!
    !
    !------------------------------------------------------------------------------


    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(in):

    INTEGER, INTENT(IN)    :: &
         ista     , & ! Station nr.
         idata    , & ! data type index: 1=vr, 2=qv, 3=z, 4=qz
         itime        ! time index, from 1 ... rs_meta(ista)%nobs_times

    CHARACTER (len=*), INTENT(in) :: hdf5_shortname ! HDF5 short name of radar data field

    !.. INTENT(out):
    REAL (KIND=dp), POINTER  :: datap(:) ! storage vector for the data
    INTEGER,        POINTER  :: zind_intp_obs(:) ! storage vector for the data index
    LOGICAL, INTENT(out)     :: ldata_avail   ! flag for success of reading

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER   :: yzroutine = 'read_field_obs_mch'
    CHARACTER (LEN=80)             :: yzerrmsg

    INTEGER                        :: ii, kk, m, n, o, ierr, nel, nobs, nobs_max
    REAL    (KIND=dp), ALLOCATABLE :: manel(:), indata3d(:,:,:)

    CHARACTER(len=cobsflen)        :: infilename
    CHARACTER(len=cobsflen),ALLOCATABLE :: infilenames(:)
    INTEGER                        :: ninfiles

    INTEGER(kind=hid_t)            :: file_id

    INTEGER, ALLOCATABLE           :: rbins(:),nazis(:)
    LOGICAL, ALLOCATABLE           :: ldata_avail_dataset(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !    Check index idata:
    IF (idata < 1 .OR. idata > ndatakind) THEN
      WRITE (*,'(a,/,a,i2,a)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           hdf5_shortname//') from an ODIM-HDF5 file:', &
           '  idatakind = ', idata, ' is not between 1 and '//cndatakind//'!'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    ! .. Return silently if there is no obs file for this quantity:
    IF (TRIM(rs_meta(ista)%obsfile(itime,idata)) == obsfile_missingname) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    infilename(:) = ' '
    infilename = rs_meta(ista)%obsfile(itime,idata)

    CALL reconstruct_file_series_mch (yzroutine, infilename, infilenames, ninfiles, ierr)
    IF (ierr /= 0) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    
    !=====================================================================================================================
    ! .. Now all filenames of the respective timestep / volume scan are known, so read their metadata:
    !=====================================================================================================================

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE (infilenames)
      RETURN
    ENDIF


    !    Define three dimensional start and count arrays
    nel      = rs_meta(ista)%nel
    !    Define max. number of possible observations:
    nobs_max = rs_meta(ista)%nra_obs * rs_meta(ista)%naz * rs_meta(ista)%nel

    !    Allocate temporary storage array for the data fields from the HDF5-files:
    ALLOCATE( indata3d(rs_meta(ista)%nra_obs, rs_meta(ista)%naz, nel) )
    indata3d = miss_value

    ! .. Set %check - fields for checking of equivalence of
    !    available azimuts and elevations for different data sets
    !    from the same station from the same time:
    ALLOCATE( manel(nel) , &
              nazis(nel) , &
              rbins(nel) , &
              ldata_avail_dataset(ninfiles)  )

    manel(:) = missing_obs
    nazis(:) = 0
    rbins(:) = 0
    ldata_avail_dataset(:) = .TRUE.

    read_loop: DO ii=1, ninfiles

      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(ii)), H5F_ACC_RDONLY_F, file_id, ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5 file '//TRIM(ydirradarin)//TRIM(infilenames(ii))
        ldata_avail = .FALSE.
      ENDIF

      CALL get_dataset_odim_h5 (TRIM(yzroutine), TRIM(ydirradarin)//TRIM(infilenames(ii)), file_id, hdf5_shortname, &
           &                    idata, 'dataset1', 'data1', rs_meta(ista), ldata_avail_dataset(ii), &
           &                    indata3d, rbins, nazis, manel, &
           &                    lcheck_azi_sorting=.FALSE., lwarn_on_missing_startazA=.FALSE., lele_allowed_to_be_wrong=.FALSE.)

      CALL h5fclose_f(file_id, ierr)

    END DO read_loop

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. If any of the elevations could be sucessfully read, we set ldata_avail=.true.:
    ldata_avail = ANY(ldata_avail_dataset)

    IF (.NOT. ldata_avail) THEN
      DEALLOCATE(manel,rbins,nazis,indata3d,infilenames,ldata_avail_dataset)
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    CALL check_and_complete_dataset ( TRIM(yzroutine), rs_data(ista), rs_meta(ista), &
         &                            itime, idata, nel, nobs_max, manel, nazis, rbins, &
         &                            zind_intp_obs, ierr)
    IF (ierr /= 0) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,infilenames,ldata_avail_dataset)
      RETURN
    END IF
    
    DEALLOCATE(manel,rbins,nazis,infilenames,ldata_avail_dataset)

    ! .. Allocate memory for the 1D output vector:
    IF (ASSOCIATED(datap)) DEALLOCATE(datap)
    NULLIFY(datap)
    ALLOCATE(datap(nobs_max))

    ! .. Put the data into the 1D output vector:
    nobs = 0
    DO o = 1, rs_meta(ista)%nel        ! Loop over Reports (one elevation per record)
      DO m = 1, rs_meta(ista)%naz      ! Loop over Azimut
        DO n = 1,rs_meta(ista)%nra_obs ! Loop over Range
          nobs = nobs + 1
          datap(nobs) = indata3d(n,m,o)
        END DO
      END DO
    END DO

    DEALLOCATE(indata3d)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_field_obs_mch

  !==============================================================================

  SUBROUTINE read_field_obs_italy(ista, idata, itime, hdf5_shortname, &
       datap, zind_intp_obs, ldata_avail)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine reads observed data for each
    !              radar station for each time in a generic fashion.
    !              However, it expects the input format to be ODIM-HDF5
    !              and that there is one file per per station per
    !              time, i.e., the way that ARPA-SIMC in Italy stores its radar obs.
    !
    ! Method: The HDF5 files for one variable and one volume scan
    !         are opened within this routine.
    !         "datap" has to be a pointer to a 1D data vector in
    !         the rs_data(ista) structure, which points to the desired
    !         input field. "ncdf_shortname" denotes the NetCDF shortname
    !         of the input field as determined from ncdump.
    !         "ldata_avail" is a flag indicating failure or success of reading.
    !         "missing_data" is the file-internal missing value for the data.
    !
    ! IS CALLED FROM "READ_OBS_RAD".
    !
    ! Caveat: At the moment, this subroutine has some side-effects through
    !         some global variables / data vectors, which are only correctly
    !         integrated if this subroutine is called from "read_obs_rad"!
    !
    !------------------------------------------------------------------------------

    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(in):

    INTEGER, INTENT(IN)    :: &
         ista     , & ! Station nr.
         idata    , & ! data type index: 1=vr, 2=qv, 3=z, 4=qz
         itime        ! time index, from 1 ... rs_meta(ista)%nobs_times

    CHARACTER (len=*), INTENT(in) :: hdf5_shortname ! HDF5 short name of radar data field

    !.. INTENT(out):
    REAL (KIND=dp), POINTER      :: datap(:) ! storage vector for the data
    INTEGER, POINTER             :: zind_intp_obs(:) ! storage vector for the data index
    LOGICAL, INTENT(out)         :: ldata_avail   ! flag for success of reading

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER   :: yzroutine = 'read_field_obs_italy'
    CHARACTER (LEN=80)             :: yzerrmsg

    INTEGER                        :: ii, m, n, o, ierr, nel, nobs_max, nobs
    REAL    (KIND=dp), ALLOCATABLE :: manel(:), indata3d(:,:,:)

    CHARACTER(len=cobsflen)        :: infilename
    CHARACTER(len=60)              :: dsetname, vardataset
    CHARACTER(len=60)              :: tmpgroup

    INTEGER(kind=hid_t)            :: file_id

    INTEGER, ALLOCATABLE           :: rbins(:),nazis(:)
    LOGICAL, ALLOCATABLE            :: ldata_avail_dataset(:)

    CHARACTER(len=cobsflen),ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    INTEGER                        :: nwords

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !    Check index idata:
    IF (idata < 1 .OR. idata > ndatakind) THEN
      WRITE (*,'(a,/,a,i2,a)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           hdf5_shortname//') from an ODIM-HDF5 file:', &
           '  idatakind = ', idata, ' is not between 1 and '//cndatakind//'!'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    ! .. Return silently if there is no obs file for this quantity:
    IF (TRIM(rs_meta(ista)%obsfile(itime,idata)) == obsfile_missingname) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    CALL split_string (TRIM(rs_meta(ista)%obsfile(itime,idata)), '/', LEN(words), words, nwords)
    IF (nwords /= 2) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Not a valid obsfile name '//TRIM(rs_meta(ista)%obsfile(itime,idata))
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    infilename(:) = ' '
    infilename = TRIM(words(1))   ! This is the true filename
    vardataset = TRIM(words(2))   ! This is the name of the data-group within the /datasetXX groups for the different elevations
    IF (ALLOCATED(words)) DEALLOCATE(words)

    !=====================================================================================================================
    ! .. Now all filenames of the respective timestep / volume scan are known, so read their metadata:
    !=====================================================================================================================

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      DEALLOCATE(manel, rbins, nazis, indata3d)
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    ENDIF

    !    Define three dimensional start and count arrays
    nel      = rs_meta(ista)%nel
    !    Define max. number of possible observations:
    nobs_max = rs_meta(ista)%nra_obs * rs_meta(ista)%naz * rs_meta(ista)%nel

    !    Allocate temporary storage array for the data fields from the HDF5-files:
    ALLOCATE( indata3d(rs_meta(ista)%nra_obs, rs_meta(ista)%naz, nel) )
    indata3d = miss_value

    ! .. Set %check - fields for checking of equivalence of
    !    available azimuts and elevations for different data sets
    !    from the same station from the same time:
    ALLOCATE( manel(nel) , &
              nazis(nel) , &
              rbins(nel) , &
              ldata_avail_dataset(nel) )

    manel(:) = missing_obs
    nazis(:) = 0
    rbins(:) = 0
    ldata_avail_dataset(:) = .TRUE.

    CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilename), H5F_ACC_RDONLY_F, file_id, ierr)
    IF (ierr /= 0) THEN
      WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5 file '//TRIM(ydirradarin)//TRIM(infilename)
      ldata_avail = .FALSE.
    ENDIF

    read_loop: DO ii=1, nel

      dsetname(:) = ' '
      WRITE (dsetname, '(i3)') ii
      dsetname = 'dataset'//TRIM(ADJUSTL(dsetname))

      CALL get_dataset_odim_h5 (TRIM(yzroutine), TRIM(ydirradarin)//TRIM(infilename), file_id, hdf5_shortname, &
           &                    idata, dsetname, vardataset, rs_meta(ista), ldata_avail_dataset(ii), &
           &                    indata3d, rbins, nazis, manel, &
           &                    lcheck_azi_sorting=.FALSE., lwarn_on_missing_startazA=.FALSE., lele_allowed_to_be_wrong=.FALSE.)


    END DO read_loop

    ! .. Close the HDF5 file:
    CALL h5fclose_f(file_id, ierr)

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. If any of the elevations could be sucessfully read, we set ldata_avail=.true.:
    ldata_avail = ANY(ldata_avail_dataset)
    
    IF (.NOT. ldata_avail) THEN
      DEALLOCATE(manel,rbins,nazis,indata3d,ldata_avail_dataset)
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    CALL check_and_complete_dataset ( TRIM(yzroutine), rs_data(ista), rs_meta(ista), &
         &                            itime, idata, nel, nobs_max, manel, nazis, rbins, &
         &                            zind_intp_obs, ierr)
    IF (ierr /= 0) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,ldata_avail_dataset)
      RETURN
    END IF    
    
    DEALLOCATE(manel,rbins,nazis,ldata_avail_dataset)

    ! .. Allocate memory for the 1D output vector:
    IF (ASSOCIATED(datap)) DEALLOCATE(datap)
    NULLIFY(datap)
    ALLOCATE(datap(nobs_max))

    ! .. Put the data into the 1D output vector:
    nobs = 0
    DO o = 1, rs_meta(ista)%nel        ! Loop over Reports (one elevation per record)
      DO m = 1, rs_meta(ista)%naz      ! Loop over Azimut
        DO n = 1,rs_meta(ista)%nra_obs ! Loop over Range
          nobs = nobs + 1
          datap(nobs) = indata3d(n,m,o)
        END DO
      END DO
    END DO

    DEALLOCATE(indata3d)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_field_obs_italy

  SUBROUTINE read_field_obs_kit_h5(ista, idata, itime, hdf5_shortname, &
       datap, zind_intp_obs, ldata_avail)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine reads observed data from the KIT C-Band radar.
    !              However, it expects the input format to be KIT HDF5
    !              and that there is one file per per station per
    !              time, i.e., the way that KIT stores its radar obs.
    !
    !    These files are special in many ways:
    !
    !    - do not follow ODIM specifications
    !    - ODIM was used as a rough orientation, but names of datasets, attributes
    !      and general file organization is/are somewhat different
    !    - azimuts are not stored in startazA list, but along other informations
    !      in a derived type (data type H5_COMPOUND) for each ray, which makes it very complicated to read
    !    - strings are not fixed-length, but variable length, which requires
    !      to use the cvar_vlen=cvar argument to read_attr_odim() instead
    !      of cvar=cvar as for all other hdf5-providers
    !    - /scanX is a group, but /scanX/moment_Y is a dataset. In ODIM, /datasetX/dataY is also a group,
    !      and the actual data are in /datasetX/dataY/data
    !    - scan9 does not have a moment_12 (= ), but all other elevations have it
    !
    ! IS CALLED FROM "READ_OBS_RAD".
    !
    ! Caveat: At the moment, this subroutine has some side-effects through
    !         some global variables / data vectors, which are only correctly
    !         integrated if this subroutine is called from "read_obs_rad"!
    !
    !------------------------------------------------------------------------------

    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(in):

    INTEGER, INTENT(IN)    :: &
         ista     , & ! Station nr.
         idata    , & ! data type index: 1=vr, 2=qv, 3=z, 4=qz
         itime        ! time index, from 1 ... rs_meta(ista)%nobs_times

    CHARACTER (len=*), INTENT(in) :: hdf5_shortname ! HDF5 short name of radar data field

    !.. INTENT(out):
    REAL (KIND=dp), POINTER      :: datap(:) ! storage vector for the data
    INTEGER, POINTER             :: zind_intp_obs(:) ! storage vector for the data index
    LOGICAL, INTENT(out)         :: ldata_avail   ! flag for success of reading

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER   :: yzroutine = 'read_field_obs_kit_h5'
    CHARACTER (LEN=80)             :: yzerrmsg

    INTEGER                        :: ii, m, n, o, ierr, nel, nobs_max, nobs
    REAL    (KIND=dp), ALLOCATABLE :: manel(:), indata3d(:,:,:)

    CHARACTER(len=cobsflen)        :: infilename
    CHARACTER(len=60)              :: dsetname, vardataset
    CHARACTER(len=60)              :: tmpgroup

    INTEGER(kind=hid_t)            :: file_id

    INTEGER, ALLOCATABLE           :: rbins(:),nazis(:)
    LOGICAL, ALLOCATABLE            :: ldata_avail_dataset(:)

    CHARACTER(len=cobsflen),ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    INTEGER                        :: nwords

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !    Check index idata:
    IF (idata < 1 .OR. idata > ndatakind) THEN
      WRITE (*,'(a,/,a,i2,a)') 'WARNING '//TRIM(yzroutine)//': problem reading data (shortname = '// &
           hdf5_shortname//') from a KIT HDF5 file:', &
           '  idatakind = ', idata, ' is not between 1 and '//cndatakind//'!'
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    ! .. Return silently if there is no obs file for this quantity:
    IF (TRIM(rs_meta(ista)%obsfile(itime,idata)) == obsfile_missingname) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    
    CALL split_string (TRIM(rs_meta(ista)%obsfile(itime,idata)), '/', LEN(words), words, nwords)
    IF (nwords /= 2) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Not a valid obsfile name '//TRIM(rs_meta(ista)%obsfile(itime,idata))
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF
    infilename(:) = ' '
    infilename = TRIM(words(1))   ! This is the true filename
    vardataset = TRIM(words(2))   ! This is the name of the dataset within the /scanXX groups for the different elevations
    IF (ALLOCATED(words)) DEALLOCATE(words)

    !=====================================================================================================================
    ! .. Now all filenames of the respective timestep / volume scan are known, so read their metadata:
    !=====================================================================================================================

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      DEALLOCATE(manel, rbins, nazis, indata3d)
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    ENDIF

    !    Define three dimensional start and count arrays
    nel      = rs_meta(ista)%nel
    !    Define max. number of possible observations:
    nobs_max = rs_meta(ista)%nra_obs * rs_meta(ista)%naz * rs_meta(ista)%nel

    !    Allocate temporary storage array for the data fields from the HDF5-files:
    ALLOCATE( indata3d(rs_meta(ista)%nra_obs, rs_meta(ista)%naz, nel) )
    indata3d = miss_value

    ! .. Set %check - fields for checking of equivalence of
    !    available azimuts and elevations for different data sets
    !    from the same station from the same time:
    ALLOCATE( manel(nel) , &
              nazis(nel) , &
              rbins(nel) , &
              ldata_avail_dataset(nel) )

    manel(:) = missing_obs
    nazis(:) = 0
    rbins(:) = 0
    ldata_avail_dataset(:) = .TRUE.

    CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilename), H5F_ACC_RDONLY_F, file_id, ierr)
    IF (ierr /= 0) THEN
      WRITE(*,*) 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5 file '//TRIM(ydirradarin)//TRIM(infilename)
      ldata_avail = .FALSE.
    ENDIF

    read_loop: DO ii=1, nel

      dsetname(:) = ' '
      WRITE (dsetname, '(i3)') ii-1  ! index starts at scan0
      dsetname = 'scan'//TRIM(ADJUSTL(dsetname))

      CALL get_dataset_kit_h5 (TRIM(yzroutine), TRIM(ydirradarin)//TRIM(infilename), file_id, hdf5_shortname, &
           &                    idata, dsetname, vardataset, rs_meta(ista), ldata_avail_dataset(ii), &
           &                    indata3d, rbins, nazis, manel, &
           &                    lcheck_azi_sorting=.TRUE., lele_allowed_to_be_wrong=.FALSE.)


    END DO read_loop

    ! .. Close the HDF5 file:
    CALL h5fclose_f(file_id, ierr)

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. If any of the elevations could be sucessfully read, we set ldata_avail=.true.:
    ldata_avail = ANY(ldata_avail_dataset)
    
    IF (.NOT. ldata_avail) THEN
      DEALLOCATE(manel,rbins,nazis,indata3d,ldata_avail_dataset)
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      RETURN
    END IF

    CALL check_and_complete_dataset ( TRIM(yzroutine), rs_data(ista), rs_meta(ista), &
         &                            itime, idata, nel, nobs_max, manel, nazis, rbins, &
         &                            zind_intp_obs, ierr)
    IF (ierr /= 0) THEN
      ldata_avail = .FALSE.
      IF (ASSOCIATED(datap)) DEALLOCATE(datap)
      NULLIFY(datap)
      DEALLOCATE(manel,rbins,nazis,indata3d,ldata_avail_dataset)
      RETURN
    END IF    
    
    DEALLOCATE(manel,rbins,nazis,ldata_avail_dataset)

    ! .. Allocate memory for the 1D output vector:
    IF (ASSOCIATED(datap)) DEALLOCATE(datap)
    NULLIFY(datap)
    ALLOCATE(datap(nobs_max))

    ! .. Put the data into the 1D output vector:
    nobs = 0
    DO o = 1, rs_meta(ista)%nel        ! Loop over Reports (one elevation per record)
      DO m = 1, rs_meta(ista)%naz      ! Loop over Azimut
        DO n = 1,rs_meta(ista)%nra_obs ! Loop over Range
          nobs = nobs + 1
          datap(nobs) = indata3d(n,m,o)
        END DO
      END DO
    END DO

    DEALLOCATE(indata3d)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_field_obs_kit_h5

  !==============================================================================
  !
  ! Generic procedure to read one dataset from an ODIM hdf5 file
  !
  !==============================================================================

  SUBROUTINE get_dataset_odim_h5 (yzroutine, infilename, file_id, hdf5_shortname, &
       &                          idata, dsetname, vardataset, rsm, ldata_avail, &
       &                          indata3d, rbins, nazis, manel, lcheck_azi_sorting, lwarn_on_missing_startazA, &
       &                          lele_allowed_to_be_wrong)

    IMPLICIT NONE

    CHARACTER(len=*),      INTENT(in)    :: yzroutine        ! Name of caller for labeling error output
    CHARACTER(len=*),      INTENT(in)    :: infilename       ! Name and path of input file for labeling error output
    CHARACTER(len=*),      INTENT(in)    :: hdf5_shortname   ! Name of variable for labeling error output
    INTEGER,               INTENT(in)    :: idata            ! Flag for the kind of data (i_dbzh, i_vrad, ...)
    CHARACTER(len=*),      INTENT(in)    :: dsetname         ! Name of the datasetXX from datasetXX/dataYY/
    CHARACTER(len=*),      INTENT(in)    :: vardataset       ! Name of the dataYY    from datasetXX/dataYY/
    INTEGER(kind=hid_t),   INTENT(in)    :: file_id          ! File id of the pre-opened hdf5 input file
    TYPE(radar_meta_type), INTENT(inout) :: rsm              ! Meta data of the radar station
    LOGICAL,               INTENT(out)   :: ldata_avail      ! Flag to indicate whether data could be read. User has to preset .true. on input; will be set to .false. if some error occurs
    REAL(kind=dp),         INTENT(out)   :: indata3d(:,:,:)  ! 3D array to pass back the volume scan data
    REAL(kind=dp),         INTENT(inout) :: manel(:)         ! Vector of elevations sorted by value, will be filled when reading the data
    INTEGER,               INTENT(inout) :: rbins(:)         ! Vector of number of range bins for each elevation, will be filled when reading the data
    INTEGER,               INTENT(inout) :: nazis(:)         ! Vector of number of azimuts for each elevation, will be filled when reading the data
    LOGICAL,               INTENT(in)    :: lcheck_azi_sorting ! Whether to check the azimut sorting and if necessary, shift the start azi to North, based on attribute "startazA"
    LOGICAL,               INTENT(in)    :: lwarn_on_missing_startazA ! Whether or not to report if "startazA" is missing in the file and azimuth cannot be checked/shifted. In this case, we have to trust the azimut sorting in the file to be referenced to North.
    LOGICAL,               INTENT(in)    :: lele_allowed_to_be_wrong
    
    INTEGER                    :: m, n, nn, o, eleind, ierr, nra, naz, n_noazi, rotsign, a1gate, ndiff
    REAL(KIND=dp)              :: tol_ele, dval, gain, offset, nodata, undetected, tmpdiff, dazi, reldazi, azmid
    INTEGER                    :: ival
    CHARACTER(len=60)          :: cpath
    INTEGER, ALLOCATABLE       :: iarr2d(:,:), azindlist(:), azichecklist(:)
    REAL(kind=dp), ALLOCATABLE :: startazA(:)

    ldata_avail = .TRUE.
    
    CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/where', 'elangle', ierr, dval=dval, silent=.TRUE.)
    ! Search for the index of the nearest elevation in the nominal scan strategy:
    IF (ierr == 0) THEN
      eleind = 1
      tmpdiff = ABS(REAL(dval,KIND=dp)-rsm%el_arr(eleind))
      DO o = 1, rsm%nel
        IF ( ABS(REAL(dval,KIND=dp)-rsm%el_arr(o)) < tmpdiff ) THEN
          eleind = o
          tmpdiff = ABS(REAL(dval,KIND=dp)-rsm%el_arr(eleind))
        END IF
      END DO

      ! THOMAS modification
      ! Modified tolerance for Italian radars
      SELECT CASE (rsm%icountry)
      CASE (i_denmark, i_arpasim)
        tol_ele = 0.2_dp
      CASE default
        tol_ele = 0.1_dp
      END SELECT
      
      IF (tmpdiff <= tol_ele) THEN
        manel(eleind) = dval
      ELSE
        eleind = missval_int
        IF (.NOT.lele_allowed_to_be_wrong) THEN
          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Strange eleind problem in '//TRIM(infilename), TRIM(dsetname), dval, tmpdiff
        END IF
        ldata_avail = .FALSE.
      END IF
    ELSE
      eleind = missval_int
      ldata_avail = .FALSE.
    END IF

    IF (eleind >= missthr_int) THEN
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/where', 'nbins', ierr, ival=ival, silent=.TRUE.)
      rbins(eleind) = ival
      IF (rbins(eleind) <= 0 .OR. rbins(eleind) > rsm%nra_obs) THEN
        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Strange nrange problem in '//TRIM(infilename)
        ldata_avail = .FALSE.
      END IF

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/where', 'nrays', ierr, ival=ival, silent=.TRUE.)
      nazis(eleind) = ival
      
      IF (rsm%icountry == i_dwd .AND. rsm%station_name == 'BON') THEN
        ! Workaround for problematic files from KIT-Cube which may contain different numbers of azimuts
        ! for the same elevations of different variables, and/or different numbers of azimuts for different
        ! elevations. Mostly this is one more azimut (361) than nominal (360). As this happens
        ! very often and these data are not used operationally, those files have not been discarded right away,
        ! but in the data reader we use a1gate to separate out the superfluous datum.
        IF (ival > rsm%naz+1) THEN
          ldata_avail = .FALSE.
        ELSE IF (ival > rsm%naz) THEN
          ! We tolerate only 1 azimut more than nominal:
          nazis(eleind) = rsm%naz
        END IF

      END IF

    END IF

    IF (ldata_avail) THEN
     
      cpath(:) = ' '
      cpath = '/'//TRIM(dsetname)//'/'//TRIM(vardataset)

      CALL read_attr_odim (file_id, TRIM(cpath)//'/what', 'gain',     ierr, dval=gain)
      CALL read_attr_odim (file_id, TRIM(cpath)//'/what', 'offset',   ierr, dval=offset)
      CALL read_attr_odim (file_id, TRIM(cpath)//'/what', 'nodata',   ierr, dval=nodata)
      CALL read_attr_odim (file_id, TRIM(cpath)//'/what', 'undetect', ierr, dval=undetected)

      CALL read_dataset_odim (file_id, TRIM(cpath)//'/data', ierr, iarr2d=iarr2d)
      IF (ALLOCATED(iarr2d)) THEN

        nra = SIZE(iarr2d,dim=1)   ! nra might be smaller than rsm%nra for certain elevations!
        naz = SIZE(iarr2d,dim=2)   ! naz might be smaller than rsm%naz for certain elevations!

        IF (rsm%icountry == i_dwd .AND. rsm%station_name == 'BON') THEN
          IF (naz > nazis(eleind)) THEN
            ! Workaround for problematic datasets from KIT-Cube:
            ! If necessary, cut out the superfluous azimut:
            ndiff = naz - nazis(eleind)
            IF (ndiff > 1) THEN
              WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine)//': Cannot fix KITC data because too large naz=', naz, &
                   ' in '//TRIM(infilename)
              ! no data will be read below, because nazis(eleind) /= naz
            ELSE
              CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/where', 'a1gate', ierr, ival=ival, silent=.TRUE.)
              a1gate = ival + 1  ! The index a1gate in the file is relative to 0, but we need it w.r.t 1
              IF (a1gate < 1 .OR. a1gate > naz) THEN
                WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine)//': Wrong value of a1gate=', a1gate, &
                     ' in '//TRIM(infilename)
                ! no data will be read below, because nazis(eleind) /= naz
              ELSE
                WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine)//': deleting superfluous azimut with index=', a1gate-1, &
                     ' from '//TRIM(cpath)//' in '//TRIM(infilename)
                ! eliminate the superfluous azimut at the true end of the PPI-scan. This is the index right before a1gate:
                IF (a1gate > 1) THEN
                  iarr2d(:,a1gate-1:naz-1) = iarr2d(:,a1gate:naz)
                END IF
                iarr2d(:,naz) = miss_value
                naz = nazis(eleind)
              END IF
            END IF
          END IF
        END IF

        IF (rbins(eleind) == nra .AND. nazis(eleind) == naz) THEN

          ALLOCATE(azindlist(naz))

          IF (lcheck_azi_sorting) THEN

            CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/how', 'startazA', ierr, dvec=startaza)

            IF (ALLOCATED(startaza)) THEN

              ! Check azimut increment:
              dazi    = 0.0_dp  ! to check the absolute mean dazi
              reldazi = 0.0_dp  ! to determine the rotation direction of the antenna
              n = 0
              DO nn=2, naz
                tmpdiff = startaza(nn) - startaza(nn-1)
                reldazi = reldazi + tmpdiff / (ABS(tmpdiff)+eps)
                IF ( ABS(tmpdiff) < 360.0_dp - 5*rsm%az_inc ) THEN
                  dazi = dazi + tmpdiff
                  n = n + 1
                END IF
              END DO
              rotsign = NINT(SIGN(1.0_dp,reldazi))
              IF (n > 0) THEN
                dazi = dazi / n
              END IF

              IF (ABS(ABS(dazi)-rsm%az_inc) <= 0.02_dp) THEN

                ! The azimuth increment is correct. Next step is to check on the azimuth reference and sorting in ascending order.
                !  For this, we do the sorting ourselves:
                !  The first azimuth should start at about north (az_start) and the following should be sorted
                !  in ascending azimut direction:

                ALLOCATE(azichecklist(rsm%naz))
                  
                azindlist(:) = -1
                azichecklist(:) = -1
                DO nn=1, naz
                  azmid = startaza(nn) + 0.5_dp*rotsign*rsm%az_inc
!!$                  azmid = MODULO(azmid, 360.0_dp)
                  n = NINT((azmid-rsm%az_start)/rsm%az_inc) + 1
                  n = MODULO(n-1, rsm%naz) + 1
!!$                  IF (n > 0 .AND. n <= rsm%naz) THEN
                  azindlist(nn)   = n
                  azichecklist(n) = nn
!!$                  END IF
                END DO

                n_noazi = 0
                DO n=1, rsm%naz
                  IF (azichecklist(n) == -1) n_noazi = n_noazi + 1
                END DO
                IF (n_noazi > 0) THEN
                  WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Some azimuts missing in '// &
                       TRIM(infilename)
                  DO n=1, rsm%naz
                    IF (azichecklist(n) == -1) THEN
                      WRITE (*,'(a,i0,a,f0.2)') ' |--> Missing azimut no. ',n, ' for elevation ', rsm%el_arr(eleind)
                    END IF
                  END DO
                END IF

                DEALLOCATE (azichecklist)
                  
              ELSE

                WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Discard file, wrong azimut increment, '// &
                     TRIM(infilename)
                ldata_avail = .FALSE.

              END IF

              DEALLOCATE(startaza)

            ELSE

              IF (lwarn_on_missing_startazA) THEN
                WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Cannot check azimuth, because startazA missing in '// &
                     TRIM(infilename)
              END IF

              ! We have to trust the azimut sorting, namely that the first azimuth is north and that the following are sorted
              !  in ascending azimut direction:
              DO nn=1, naz
                azindlist(nn) = nn
              END DO

            END IF   ! startazA allocated

          ELSE

            ! We trust the azimut sorting, namely that the first azimuth is north and that the following are sorted
            !  in ascending azimut direction:
            DO nn=1, naz
              azindlist(nn) = nn
            END DO

          END IF  ! lcheck_azi_sorting

          IF (ldata_avail) THEN

            ! .. Keep track of this elevation and add it to the list of present elevations:
            IF (rsm%nel_present > 0) THEN
              IF (ALL(rsm%ind_ele_present(1:rsm%nel_present) /= eleind)) THEN
                rsm%nel_present = rsm%nel_present + 1
                rsm%ind_ele_present(rsm%nel_present) = eleind
              END IF
            ELSE
              rsm%nel_present = 1
              rsm%ind_ele_present(rsm%nel_present) = eleind         
            END IF
            
            SELECT CASE (idata)
            CASE (i_vrad)

              ! .. This is radial wind:

              ! NOTE: ras07 and ras11 files of DWD will lead to slight differences, because gain, offset and Nyquist velocity (dataset1/how/NI)
              !       are different among these file formats. Also, the wavelenght attributes are different and have slightly different values.

              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == NINT(undetected) .OR. iarr2d(m,nn) == NINT(nodata)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

              ! THOMAS modification: all Gattatico radial wind values are set to 
              ! nodata (missing_obs) due to a DPC coding problem whereby all values are
              ! undetect or equal to the maximum
              IF (rsm%station_id == 16199) THEN
                DO n = 1, naz
                  DO m = 1, nra
                      indata3d(m, n, eleind) = missing_obs
                  END DO
                END DO
              END IF

            CASE (i_dbzh)

              ! .. This is horizontal reflectivity:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == NINT(undetected)) THEN
                      indata3d(m, n, eleind) = zero_value
                    ELSE IF (iarr2d(m,nn) == NINT(nodata)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_zdr)

              ! .. This is differential reflectivity:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == NINT(undetected)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE IF (iarr2d(m,nn) == NINT(nodata)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                    ! Some crude quality control:
                    IF (ABS(indata3d(m, n, eleind)) > 0.95_dp*ABS(offset)) THEN
                      ! Inhibit noisy values around the extreme end near the undetected or nodata flag:
                      indata3d(m, n, eleind) = missing_obs
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_kdp, i_phidp)

              ! .. This is specific or total differential phase shift:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == NINT(undetected)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE IF (iarr2d(m,nn) == NINT(nodata)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_rhv)

              ! .. This is H/V correllation coefficient:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == NINT(undetected)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE IF (iarr2d(m,nn) == NINT(nodata)) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE default

              WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Discard file, wrong idata for hdf5_shortname '''// &
                   TRIM(hdf5_shortname)//''' when reading the 2D data from '// &
                   TRIM(infilename)
              ldata_avail = .FALSE.

            END SELECT

          END IF   ! ldata_avail

          DEALLOCATE(azindlist)

        ELSE

          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Discard file, wrong range or azimut dimension '// &
               'when reading the 2D data from '//TRIM(infilename)
          ldata_avail = .FALSE.

        END IF

        DEALLOCATE(iarr2d)

      ELSE

        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed reading the 2D data from '// &
             TRIM(infilename)
        ldata_avail = .FALSE.

      END IF
    END IF

  END SUBROUTINE get_dataset_odim_h5
  
  !==============================================================================
  !
  ! Procedure to read one dataset from a KIT hdf5 file
  !
  !==============================================================================

  SUBROUTINE get_dataset_kit_h5 (yzroutine, infilename, file_id, hdf5_shortname, &
       &                         idata, dsetname, vardataset, rsm, ldata_avail, &
       &                         indata3d, rbins, nazis, manel, lcheck_azi_sorting, &
       &                         lele_allowed_to_be_wrong)

    USE ISO_FORTRAN_ENV  ! for portable kind type identifiers int64, int32, int16, int8, REAL32, REAL64, REAL128
    
    IMPLICIT NONE

    CHARACTER(len=*),      INTENT(in)    :: yzroutine        ! Name of caller for labeling error output
    CHARACTER(len=*),      INTENT(in)    :: infilename       ! Name and path of input file for labeling error output
    CHARACTER(len=*),      INTENT(in)    :: hdf5_shortname   ! Name of variable for labeling error output
    INTEGER,               INTENT(in)    :: idata            ! Flag for the kind of data (i_dbzh, i_vrad, ...)
    CHARACTER(len=*),      INTENT(in)    :: dsetname         ! Name of the scanXX    from scanXX/moment_YY/
    CHARACTER(len=*),      INTENT(in)    :: vardataset       ! Name of the moment_YY from scanXX/moment_YY/
    INTEGER(kind=hid_t),   INTENT(in)    :: file_id          ! File id of the pre-opened hdf5 input file
    TYPE(radar_meta_type), INTENT(inout) :: rsm              ! Meta data of the radar station
    LOGICAL,               INTENT(out)   :: ldata_avail      ! Flag to indicate whether data could be read. User has to preset .true. on input; will be set to .false. if some error occurs
    REAL(kind=dp),         INTENT(out)   :: indata3d(:,:,:)  ! 3D array to pass back the volume scan data
    REAL(kind=dp),         INTENT(inout) :: manel(:)         ! Vector of elevations sorted by value, will be filled when reading the data
    INTEGER,               INTENT(inout) :: rbins(:)         ! Vector of number of range bins for each elevation, will be filled when reading the data
    INTEGER,               INTENT(inout) :: nazis(:)         ! Vector of number of azimuts for each elevation, will be filled when reading the data
    LOGICAL,               INTENT(in)    :: lcheck_azi_sorting ! Whether to check the azimut sorting and if necessary, shift the start azi to North, based on attribute "startazA"
    LOGICAL,               INTENT(in)    :: lele_allowed_to_be_wrong
    
    INTEGER                    :: m, n, nn, o, eleind, ierr, nra, naz, n_noazi, rotsign, a1gate, ndiff, nodata, undetected
    REAL(KIND=dp)              :: tol_ele, dval, gain, offset, tmpdiff, dazi, reldazi, azmid, dyn_range_min, dyn_range_max
    INTEGER                    :: ival
    CHARACTER(len=60)          :: cpath
    INTEGER, ALLOCATABLE       :: iarr2d(:,:), azindlist(:), azichecklist(:)
    REAL(kind=dp), ALLOCATABLE :: startazA(:)

    INTEGER, PARAMETER :: i8 = selected_int_KIND(15)  ! cover at least +-10^15
    INTEGER, PARAMETER :: i2 = selected_int_KIND(3)   ! cover at least +-10^3
    INTEGER, PARAMETER :: i1 = selected_int_kind(1)   ! cover at least +-10
    
    ! Data type for reading the compount H5T_COMPOUND type ray_header:
!!$ DATASET "/scan0/ray_header" {
!!$    DATATYPE  H5T_COMPOUND {
!!$      H5T_IEEE_F64LE "azimuth_start";
!!$      H5T_IEEE_F64LE "azimuth_stop";
!!$      H5T_IEEE_F64LE "elevation_start";
!!$      H5T_IEEE_F64LE "elevation_stop";
!!$      H5T_STD_I64LE  "timestamp";
!!$      H5T_IEEE_F64LE "az_speed";
!!$      H5T_IEEE_F64LE "el_speed";
!!$      H5T_STD_I64LE  "tm_stop";
!!$      H5T_STD_U32LE  "cur_prf";
!!$      H5T_STD_U16LE  "num_pulses";
!!$      H5T_IEEE_F64LE "burst_power";
!!$      H5T_IEEE_F64LE "burst_freq";
!!$      H5T_STD_U8LE   "afc_mode";
!!$      H5T_STD_U32LE  "afc_status";
!!$      H5T_STD_U32LE  "ifd_power_flags";
!!$      H5T_STD_U32LE  "adc_overflow_flags";
!!$      H5T_IEEE_F32LE "pci_dma_tranfer_rate";
!!$      H5T_IEEE_F32LE "pci_dma_fifo_fill";
!!$      H5T_IEEE_F32LE "tcp_transfer_rate";
!!$      H5T_IEEE_F32LE "host_throughput";
!!$      H5T_STD_U32LE  "shortpulse_bins";
!!$      H5T_IEEE_F32LE "pc_gain";             
!!$    }

    TYPE t_ray_header
      ! accuracy of data types has to be precise here, because they will be mapped to
      ! the corresponding C-types
      SEQUENCE
      REAL(real64)   :: azimuth_start
      REAL(real64)   :: azimuth_stop
      REAL(real64)   :: elevation_start
      REAL(real64)   :: elevation_stop
      INTEGER(int64) :: timestamp
      REAL(real64)   :: az_speed
      REAL(real64)   :: el_speed
      INTEGER(int64) :: tm_stop
      INTEGER(int32) :: cur_prf
      INTEGER(int16) :: num_pulses
      REAL(real64)   :: burst_power
      REAL(real64)   :: burst_freq
      INTEGER(int8)  :: afc_mode
      INTEGER(int32) :: afc_status
      INTEGER(int32) :: ifd_power_flags
      INTEGER(int32) :: adc_overflow_flags
      REAL(real32)   :: pci_dma_tranfer_rate
      REAL(real32)   :: pci_dma_fifo_fill 
      REAL(real32)   :: tcp_transfer_rate
      REAL(real32)   :: host_throughput
      INTEGER(int32) :: shortpulse_bins
      REAL(real32)   :: pc_gain             
    END TYPE t_ray_header

    TYPE(t_ray_header), ALLOCATABLE, TARGET :: ray_header(:)
    
    ldata_avail = .TRUE.
    
    CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/how', 'elevation', ierr, dval=dval, silent=.TRUE.)
    ! Search for the index of the nearest elevation in the nominal scan strategy:
    IF (ierr == 0) THEN
      eleind = 1
      tmpdiff = ABS(REAL(dval,KIND=dp)-rsm%el_arr(eleind))
      DO o = 1, rsm%nel
        IF ( ABS(REAL(dval,KIND=dp)-rsm%el_arr(o)) < tmpdiff ) THEN
          eleind = o
          tmpdiff = ABS(REAL(dval,KIND=dp)-rsm%el_arr(eleind))
        END IF
      END DO

      tol_ele = 0.1_dp
      
      IF (tmpdiff <= tol_ele) THEN
        manel(eleind) = dval
      ELSE
        eleind = missval_int
        IF (.NOT.lele_allowed_to_be_wrong) THEN
          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Strange eleind problem in '//TRIM(infilename), TRIM(dsetname), dval, tmpdiff
        END IF
        ldata_avail = .FALSE.
      END IF
    ELSE
      eleind = missval_int
      ldata_avail = .FALSE.
    END IF

    IF (eleind >= missthr_int) THEN
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/how', 'bin_count', ierr, ival=ival, silent=.TRUE.)
      rbins(eleind) = ival
      IF (rbins(eleind) <= 0 .OR. rbins(eleind) > rsm%nra_obs) THEN
        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Strange nrange problem in '//TRIM(infilename)
        ldata_avail = .FALSE.
      END IF

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname)//'/how', 'ray_count', ierr, ival=ival, silent=.TRUE.)
      nazis(eleind) = ival
      
      ! Workaround for problematic files from KIT which may contain different numbers of azimuts
      ! for the same elevations of different variables, and/or different numbers of azimuts for different
      ! elevations. Mostly this is one more azimut (361) than nominal (360). As this happens
      ! very often and these data are not used operationally, those files have not been discarded right away,
      ! but in the data reader we use a1gate to separate out the superfluous datum.
      IF (ival > rsm%naz+1) THEN
        ldata_avail = .FALSE.
      ELSE IF (ival > rsm%naz) THEN
        ! We tolerate only 1 azimut more than nominal:
        nazis(eleind) = rsm%naz
      END IF

    END IF

    IF (ldata_avail) THEN
     
      cpath(:) = ' '
      cpath = '/'//TRIM(dsetname)//'/'//TRIM(vardataset)

!!$ attribute  /scan0/moment_0/dyn_range_min    = offset
!!$ attribute  /scan0/moment_0/dyn_range_max    --> gain = (dyn_range_max-dyn_range_min) / 255

      CALL read_attr_odim (file_id, TRIM(cpath), 'dyn_range_max', ierr, rval=dyn_range_max) ! dtype = 32 bit!
      CALL read_attr_odim (file_id, TRIM(cpath), 'dyn_range_min', ierr, rval=dyn_range_min) ! dtype = 32 bit!
      
      gain   = (dyn_range_max - dyn_range_min) / 255.0_dp
      offset =  dyn_range_min

      nodata = 0
      undetected = 0
      
      CALL read_dataset_odim (file_id, TRIM(cpath), ierr, iarr2d=iarr2d)
      IF (ALLOCATED(iarr2d)) THEN

        nra = SIZE(iarr2d,dim=1)   ! nra might be smaller than rsm%nra for certain elevations!
        naz = SIZE(iarr2d,dim=2)   ! naz might be smaller than rsm%naz for certain elevations!

        IF (naz > nazis(eleind)) THEN
          ! Workaround for problematic datasets from KIT:
          ! If necessary, cut out the superfluous azimut:
          ndiff = naz - nazis(eleind)
          IF (ndiff > 1) THEN
            WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine)//': Cannot fix KITC data because too large naz=', naz, &
                 ' in '//TRIM(infilename)
            ! no data will be read below, because nazis(eleind) /= naz
          ELSE
            WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine)//': deleting all azimuts with index > ', nazis(eleind), &
                 ' from '//TRIM(cpath)//' in '//TRIM(infilename)
            ! eliminate the superfluous azimut at the true end of the PPI-scan:
            iarr2d(:,naz) = miss_value
            naz = nazis(eleind)
          END IF
        END IF

        IF (rbins(eleind) == nra .AND. nazis(eleind) == naz) THEN

          ALLOCATE(azindlist(naz))

          IF (lcheck_azi_sorting) THEN
            
            CALL read_ray_header_kit (file_id, '/'//TRIM(dsetname), ray_header, ierr)
            ! ray_header(:) should now be allocated with as many elements as the original ray_count in the file

            IF (ierr == 0 .AND. ALLOCATED(ray_header)) THEN

              ALLOCATE (startazA(naz))
              DO nn=1, naz
                startazA(nn) = ray_header(nn)%azimuth_start
              END DO

              ! Check azimut increment:
              dazi    = 0.0_dp  ! to check the absolute mean dazi
              reldazi = 0.0_dp  ! to determine the rotation direction of the antenna
              n = 0
              DO nn=2, naz
                tmpdiff = startaza(nn) - startaza(nn-1)
                reldazi = reldazi + tmpdiff / (ABS(tmpdiff)+eps)
                IF ( ABS(tmpdiff) < 360.0_dp - 5*rsm%az_inc ) THEN
                  dazi = dazi + tmpdiff
                  n = n + 1
                END IF
              END DO
              rotsign = NINT(SIGN(1.0_dp,reldazi))
              IF (n > 0) THEN
                dazi = dazi / n
              END IF

              IF (ABS(ABS(dazi)-rsm%az_inc) <= 0.02_dp) THEN

                ! The azimuth increment is correct. Next step is to check on the azimuth reference and sorting in ascending order.
                !  For this, we do the sorting ourselves:
                !  The first azimuth should start at about north (az_start) and the following should be sorted
                !  in ascending azimut direction:

                ALLOCATE(azichecklist(rsm%naz))
                  
                azindlist(:) = -1
                azichecklist(:) = -1
                DO nn=1, naz
                  azmid = startaza(nn) + 0.5_dp*rotsign*rsm%az_inc
!!$                  azmid = MODULO(azmid, 360.0_dp)
                  n = NINT((azmid-rsm%az_start)/rsm%az_inc) + 1
                  n = MODULO(n-1, rsm%naz) + 1
!!$                  IF (n > 0 .AND. n <= rsm%naz) THEN
                  azindlist(nn)   = n
                  azichecklist(n) = nn
!!$                  END IF
                END DO

                n_noazi = 0
                DO n=1, rsm%naz
                  IF (azichecklist(n) == -1) n_noazi = n_noazi + 1
                END DO
                IF (n_noazi > 0) THEN
                  WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Some azimuts missing in '// &
                       TRIM(infilename)
                  DO n=1, rsm%naz
                    IF (azichecklist(n) == -1) THEN
                      WRITE (*,'(a,i0,a,f0.2)') ' |--> Missing azimut no. ',n, ' for elevation ', rsm%el_arr(eleind)
                    END IF
                  END DO
                END IF

                DEALLOCATE (azichecklist)
                  
              ELSE

                WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Discard ', TRIM(cpath),', wrong azimut increment, '// &
                     TRIM(infilename)
                ldata_avail = .FALSE.

              END IF

              DEALLOCATE(startaza, ray_header)

            ELSE

              WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Cannot check azimuth, because /', &
                   TRIM(dsetname),'/ray_header could not be read from '//TRIM(infilename)

              ldata_avail = .FALSE.
                
            END IF   ! startazA allocated

          ELSE

            ! We trust the azimut sorting, namely that the first azimuth is north and that the following are sorted
            !  in ascending azimut direction:
            DO nn=1, naz
              azindlist(nn) = nn
            END DO

          END IF  ! lcheck_azi_sorting

          IF (ldata_avail) THEN

            ! .. Keep track of this elevation and add it to the list of present elevations:
            IF (rsm%nel_present > 0) THEN
              IF (ALL(rsm%ind_ele_present(1:rsm%nel_present) /= eleind)) THEN
                rsm%nel_present = rsm%nel_present + 1
                rsm%ind_ele_present(rsm%nel_present) = eleind
              END IF
            ELSE
              rsm%nel_present = 1
              rsm%ind_ele_present(rsm%nel_present) = eleind         
            END IF
            
            SELECT CASE (idata)
            CASE (i_vrad)

              ! .. This is radial wind:

              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == 0) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*255.0_dp/254.0_dp*(iarr2d(m,nn)-1) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_dbzh)

              ! .. This is horizontal reflectivity:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == 0) THEN
                      indata3d(m, n, eleind) = zero_value
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_zdr)

              ! .. This is differential reflectivity:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == 0) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                    ! Some crude quality control:
                    IF (ABS(indata3d(m, n, eleind)) > 0.95_dp*ABS(offset)) THEN
                      ! Inhibit noisy values around the extreme end near the undetected or nodata flag:
                      indata3d(m, n, eleind) = missing_obs
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_kdp, i_phidp)

              ! .. This is specific or total differential phase shift:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == 0) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE (i_rhv)

              ! .. This is H/V correllation coefficient:
              DO nn = 1, naz
                n = azindlist(nn)
                IF (n > 0) THEN
                  DO m = 1, nra
                    IF (iarr2d(m,nn) == 0) THEN
                      indata3d(m, n, eleind) = missing_obs
                    ELSE
                      indata3d(m, n, eleind) = gain*iarr2d(m,nn) + offset
                    END IF
                  END DO
                END IF
              END DO

            CASE default

              WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Discard file, wrong idata for hdf5_shortname '''// &
                   TRIM(hdf5_shortname)//''' when reading the 2D data from '// &
                   TRIM(infilename)
              ldata_avail = .FALSE.

            END SELECT

          END IF   ! ldata_avail

          DEALLOCATE(azindlist)

        ELSE

          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Discard file, wrong range or azimut dimension '// &
               'when reading the 2D data from '//TRIM(infilename)
          ldata_avail = .FALSE.

        END IF

        DEALLOCATE(iarr2d)

      ELSE

        WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': Failed reading the 2D data from '// &
             TRIM(infilename)
        ldata_avail = .FALSE.

      END IF
    END IF

  CONTAINS

    SUBROUTINE read_ray_header_kit (file_id, grpname, h, ierr)

      INTEGER(kind=hid_t), INTENT(in)  :: file_id
      CHARACTER(len=*),    INTENT(in)  :: grpname   ! path to ray_header
      TYPE(t_ray_header),  INTENT(inout), ALLOCATABLE, TARGET :: h(:)
      INTEGER,             INTENT(out) :: ierr


      INTEGER            :: error, ierror, ndims
      INTEGER(kind=hid_t):: dset_id, type_id, dsp_id, sample_type_id

      CHARACTER(len=60) :: dsetname

      INTEGER(hsize_t) :: datadims(1), maxdatadims(1)
      TYPE(C_PTR) :: f_ptr

      CHARACTER(len=*), PARAMETER :: yzroutine = 'get_dataset_kit_h5::read_ray_header_kit'

      ierr = 0
      IF (ALLOCATED(h)) DEALLOCATE(h)

      dsetname(:) = ' '
      dsetname    = TRIM(grpname)//'/ray_header'
      
      ! .. Get dimensions of the data set:
      error = 0
      CALL h5oopen_f(file_id, TRIM(dsetname), dset_id, ierror)
      error = error + ierror
      CALL h5dget_type_f        (dset_id, type_id, ierror)
      error = error + ierror
      CALL h5dget_space_f       (dset_id, dsp_id, ierror)
      error = error + ierror
      CALL h5sget_simple_extent_ndims_f (dsp_id, ndims, ierror)
      error = error + ierror
      IF (ndims /= 1) THEN
        WRITE (*, '(a,i0,a)') 'ERROR '//TRIM(yzroutine)//': ray_header dimension=', ndims, &
             ', but should be 1, '//TRIM(grpname)
        CALL h5sclose_f (dsp_id , ierror)
        CALL h5tclose_f (type_id, ierror)
        CALL h5oclose_f (dset_id, ierror)
        ierr = 1
        RETURN
      END IF
      CALL h5sget_simple_extent_dims_f (dsp_id, datadims, maxdatadims, ierror)
      ! this leads to ierror > 0, but datadims seems to be correctly determined
      
      CALL h5sclose_f (dsp_id , ierror)
      CALL h5tclose_f (type_id, ierror)
      CALL h5oclose_f (dset_id, ierror)

      IF (error > 0) THEN
        WRITE (*, '(a,i0,a)') 'ERROR '//TRIM(yzroutine)//': section 1 hdf5 problems reading '//TRIM(dsetname)
        ierr = error
        RETURN
      END IF

      
      ! .. Now get the ray_header data set:

      ALLOCATE(h(datadims(1)))

      ! NOTE: for the following to work, the hdf5 library version must be >= 1.8.8 and
      !       it has to be build with both configure options --enable-fortran --enable-fortran2003 !

      ! .. Open the dataset:
      CALL H5Dopen_f(file_id, TRIM(dsetname), dset_id, ierror)

      ! .. Create the H5T_COMPOUND datatype matching the storage size of type ray_header:
      CALL H5Tcreate_f(H5T_COMPOUND_F, h5offsetof(C_LOC(h(1)), C_LOC(h(2))), sample_type_id, ierror)

      ! .. Insert the memory offsets of the type components of type ray_header
      !     within the memory of the NEW compound TYPE :
      CALL H5Tinsert_f( sample_type_id, "azimuth_start", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%azimuth_start)), &
           h5kind_to_type(KIND(h(1)%azimuth_start),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "azimuth_stop", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%azimuth_stop)), &
           h5kind_to_type(KIND(h(1)%azimuth_stop),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "elevation_start", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%elevation_start)), &
           h5kind_to_type(KIND(h(1)%elevation_start),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "elevation_stop", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%elevation_stop)), &
           h5kind_to_type(KIND(h(1)%elevation_stop),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "timestamp", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%timestamp)), &
           h5kind_to_type(KIND(h(1)%timestamp),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "az_speed", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%az_speed)), &
           h5kind_to_type(KIND(h(1)%az_speed),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "el_speed", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%el_speed)), &
           h5kind_to_type(KIND(h(1)%el_speed),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "tm_stop", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%tm_stop)), &
           h5kind_to_type(KIND(h(1)%tm_stop),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "cur_prf", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%cur_prf)), &
           h5kind_to_type(KIND(h(1)%cur_prf),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "num_pulses", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%num_pulses)), &
           h5kind_to_type(KIND(h(1)%num_pulses),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "burst_power", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%burst_power)), &
           h5kind_to_type(KIND(h(1)%burst_power),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "burst_freq", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%burst_freq)), &
           h5kind_to_type(KIND(h(1)%burst_freq),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "afc_mode", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%afc_mode)), &
           h5kind_to_type(KIND(h(1)%afc_mode),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "afc_status", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%afc_status)), &
           h5kind_to_type(KIND(h(1)%afc_status),H5_INTEGER_KIND), ierror)    
      CALL H5Tinsert_f( sample_type_id, "ifd_power_flags", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%ifd_power_flags)), &
           h5kind_to_type(KIND(h(1)%ifd_power_flags),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "adc_overflow_flags", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%adc_overflow_flags)), &
           h5kind_to_type(KIND(h(1)%adc_overflow_flags),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "pci_dma_tranfer_rate", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%pci_dma_tranfer_rate)), &
           h5kind_to_type(KIND(h(1)%pci_dma_tranfer_rate),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "pci_dma_fifo_fill", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%pci_dma_fifo_fill)), &
           h5kind_to_type(KIND(h(1)%pci_dma_fifo_fill),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "tcp_transfer_rate", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%tcp_transfer_rate)), &
           h5kind_to_type(KIND(h(1)%tcp_transfer_rate),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "host_throughput", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%host_throughput)), &
           h5kind_to_type(KIND(h(1)%host_throughput),H5_REAL_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "shortpulse_bins", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%shortpulse_bins)), &
           h5kind_to_type(KIND(h(1)%shortpulse_bins),H5_INTEGER_KIND), ierror)
      CALL H5Tinsert_f( sample_type_id, "pc_gain", &
           h5offsetof(C_LOC(h(1)), C_LOC(h(1)%pc_gain)), &
           h5kind_to_type(KIND(h(1)%pc_gain),H5_REAL_KIND), ierror)

      ! .. Read the actual dataset ray_header into the type instance h:
      f_ptr = C_LOC(h(1))
      CALL H5Dread_f(dset_id, sample_type_id, f_ptr, ierr)
      IF (ierr > 0) THEN
        WRITE (*, '(a,i0,a)') 'ERROR '//TRIM(yzroutine)//': section 2 hdf5 problems reading '//TRIM(dsetname)
      END IF

      ! .. Clean up:
      CALL H5Tclose_f(sample_type_id, ierror)
      CALL h5dclose_f (dset_id, ierror)


    END SUBROUTINE read_ray_header_kit
    
  END SUBROUTINE get_dataset_kit_h5
  
#endif

  !==============================================================================

  SUBROUTINE check_and_complete_dataset ( yzroutine, rsd, rsm, itime, idata, nel, nobs_max, &
       &                                  manel, nazis, rbins, zind_intp_obs, ierr)

    CHARACTER(len=*),      INTENT(in)    :: yzroutine
    TYPE(radar_data_type), INTENT(inout) :: rsd
    TYPE(radar_meta_type), INTENT(in)    :: rsm
    INTEGER,               INTENT(in)    :: itime, idata, nel, nobs_max
    real(kind=dp),         INTENT(in)    :: manel(:)
    INTEGER,               INTENT(in)    :: nazis(:), rbins(:)
    INTEGER,               POINTER       :: zind_intp_obs(:)
    INTEGER,               INTENT(out)   :: ierr

    INTEGER :: m, n, o, nobs, tol, zcheck(3)

    ierr = 0
    
    IF (.NOT.lequal_azi_alldatasets .OR. MAXVAL(rsd%check) == -99) THEN

      ! .. Either the %check - fields have not been initialized so far, i.e., first call
      !    for the actual time, or lequal_azi_alldatasets = .false., i.e., which means
      !    that the following is called for every dataset and we do not assume
      !    that azimuts are exactly equal for all datasets.

      rsd%check(1) = SUM(NINT(manel(1:nel)*100.0_dp), mask=(manel > miss_threshold))
      rsd%check(2) = SUM(nazis(1:nel), mask=(manel > miss_threshold))
      rsd%check(3) = SUM(rbins(1:nel), mask=(manel > miss_threshold))

      ! .. Allocate and calculate the vector for storing the running index from
      !    which the indices for range, azimut and elevation can be restored:
      IF (ASSOCIATED(zind_intp_obs)) DEALLOCATE(zind_intp_obs)
      NULLIFY(zind_intp_obs)

      IF (lequal_azi_alldatasets) THEN
        IF (ASSOCIATED(rsd%ind_intp_obs)) DEALLOCATE(rsd%ind_intp_obs)
        NULLIFY(rsd%ind_intp_obs)
        ALLOCATE(rsd%ind_intp_obs(nobs_max))
        rsd%ind_intp_obs = missval_int
        zind_intp_obs => rsd%ind_intp_obs
      ELSE
        ALLOCATE(zind_intp_obs(nobs_max))
        zind_intp_obs = missval_int
      END IF

      nobs  = 0
      DO o = 1, rsm%nel        ! Loop over elevations
        DO m = 1, rsm%naz      ! Loop over Azimut
          DO n = 1,rsm%nra_obs ! Loop over Range
            nobs = nobs + 1
            zind_intp_obs(nobs) = m + (n-1)*rsm%naz + (o-1)*rsm%nra_obs*rsm%naz
          END DO
        END DO
      END DO

      IF (lequal_azi_alldatasets) THEN
        DO n=1, ndatakind
          IF (rsm%obsfile(itime,n) /= obsfile_missingname) THEN
            rsd%nobs_obs(n)     = nobs
          ELSE
            rsd%nobs_obs(n)     = 0
          END IF
        END DO
      ELSE
        rsd%nobs_obs(idata)     = nobs
      END IF

    ELSE

      ! .. This is a follow up call for an additional data set, therefore compare
      !    azimuts and elevations of this data set with %check:

      ! .. Some azimut- and/or elevation values in different data type files might deviate
      !    by 0.01 degree (rounding error???), so we need a little tolerance
      !    in elevation- and azimut checking:

      tol = nel*1  ! tolerance in the sum of azimuts/elevations = nel * 0.01 (nel*1 because of factor 100 below!)

      zcheck(1) = SUM(NINT(manel(1:nel)*100.0_dp), mask=(manel > miss_threshold))
      zcheck(2) = SUM(nazis(1:nel), mask=(manel > miss_threshold))
      zcheck(3) = SUM(rbins(1:nel), mask=(manel > miss_threshold))

      IF ( ABS( rsd%check(1)-zcheck(1) ) > tol .OR. &
                rsd%check(2)   /= zcheck(2)    .OR. &
                rsd%check(3)   /= zcheck(3)         &
         ) THEN

        WRITE(*,'(a,i4,i8,3(/,a,2i12))') &
             TRIM(yzroutine)//': WARNING '//TRIM(yzroutine)//': Data set inconsistent '// &
             'from same radar and time: station: ', rsm%station_id, &
             '  elev checksum: ', zcheck(1), rsd%check(1) , &
             '  azim checksum: ', zcheck(2), rsd%check(2) , &
             '  rbin checksum: ', zcheck(3), rsd%check(3)

        ierr = 1
        RETURN
        
      END IF

      IF (lequal_azi_alldatasets) THEN
        IF (ASSOCIATED(zind_intp_obs)) DEALLOCATE(zind_intp_obs)
        NULLIFY(zind_intp_obs)
        zind_intp_obs => rsd%ind_intp_obs
      END IF

    END IF

  END SUBROUTINE check_and_complete_dataset

  !==============================================================================
  
  SUBROUTINE read_obs_rad ( ista, time_mod, action )

    !------------------------------------------------------------------------------
    !
    ! Description: Read the radar data for
    !              the station with number "ista" and for model time
    !              "time_mod" (seconds since model start.
    !
    ! Method: - Open input files for vr, qv, z, qz
    !         - allocate necessary data vector components in rs_data-structure
    !         - read the data into this structure
    !         - round first azimut value to nearest multiple of dazi
    !
    ! MISSING: Interpolation to regular azimut grid!
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)          :: ista
    REAL (KIND=dp), INTENT(IN)   :: time_mod       ! seconds since model start
    CHARACTER(len=*), INTENT(in) :: action

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'read_obs_rad'
    INTEGER            :: i, j, n, ikind, start_rec, end_rec, irep, irep_obs, status, nwords, ierr
    REAL(kind=dp)      :: missing_vr, missing_qv, missing_z, missing_qz
    LOGICAL            :: file_identified
    INTEGER            :: nel, eleind, nobs_el, start_ind, end_ind ! THOMAS modification

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------


    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine),' ',action, ' on proc ', my_radar_id

    IF (action == 'read') THEN

      ! .. for safety, clean and deallocate leftover memory from pointers in rs_data(ista):
      CALL cleanup_read_obs_rad_ista ()

      ! .. search for the index of current time_mod in the obs_times:
      irep = missval_int
      IF (rs_meta(ista)%nobs_times > 0) THEN
        irep = get_obstime_ind_of_currtime ( rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) ) ! missval_int if not found
      END IF
      
      ! .. search for the index of current time_mod in the obs_times_obs:
      irep_obs = missval_int
      IF (rs_meta(ista)%nobs_times_obs > 0) THEN
        irep_obs = get_obstime_ind_of_currtime ( rs_meta(ista)%obs_times_obs(1:rs_meta(ista)%nobs_times_obs) ) ! missval_int if not found
      END IF
      
      ! .. Initialize data present flags, data pointers and nobs_obs:
      rs_meta(ista)%lobs_avail(irep,:) = .FALSE.

      ! .. for safety:
      rs_data(ista)%nobs_obs(:) = 0


      file_identified = .FALSE.

      format_choice: IF (irep == missval_int .OR. irep_obs == missval_int) THEN

        WRITE (*,'(a,i6.6,a,f12.1,a)') TRIM(yzroutine)//': Read nothing for ', rs_meta(ista)%station_id, &
             ' at model time ', time_mod, ' s because no radar obs data found!'

      ELSE

        ! .. Data found for the particular timestep, so read them:

        ! .. Determine file format from file name:
        !     Search for a file name /= obsfile_missingname and determine the file format as in read_meta_info_all():
        ikind = 0
        DO j=1, ndatakind
          IF (rs_meta(ista)%obsfile(irep,j) /= obsfile_missingname) THEN
            ikind = j
            EXIT
          END IF
        END DO

        IF (ikind == 0 .AND. ABS(rs_meta(ista)%obs_times(irep)-rs_meta(ista)%obs_times_obs(irep_obs)) < 1e-1_dp ) THEN
          WRITE (*,'(a,i6.6,a)') 'WARNING '//TRIM(yzroutine)//': no data files for station ', &
               rs_meta(ista)%station_id, ' at time '//rs_meta(ista)%obs_cdate(irep)//&
               '! Should not have happended at this point! No obs are read!'
          RETURN
        END IF

        ! .. Initialize the list of present elevations, will be set correctly by the below readers:
        rs_meta(ista)%nel_present = 0
        rs_meta(ista)%ind_ele_present(:) = missval_int

#ifdef HDF5_RADAR_INPUT
        IF ( rs_meta(ista)%icountry == i_dwd .AND. &
             ( TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_smss) .OR. &
               TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_mmss) ) ) THEN

          ! Either DWD single-sweep multi-moment hdf5 files of filename-convention
          !     "rasXX-pcpng01_sweeph5allm_any_00-2020020300053400-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5allm_any_00-2020020300055800-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5allm_any_01-2020020300062100-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5allm_any_02-2020020300064500-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5allm_any_03-2020020300070800-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5allm_any_04-2020020300073100-neu-10557-hd5"
          !
          ! or DWD single-sweep single-moment hdf5 files of filename-convention
          !     "rasXX-pcpng01_sweeph5onem_dbzh_00-2020020300053400-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5onem_dbzh_00-2020020300055800-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5onem_dbzh_01-2020020300062100-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5onem_dbzh_02-2020020300064400-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5onem_dbzh_03-2020020300070800-neu-10557-hd5"
          !     "rasXX-vol5minng01_sweeph5onem_dbzh_04-2020020300073100-neu-10557-hd5"

          rs_data(ista)%check(:) = -99

          ! .. Reading of data:
          !      ( read_field_obs_dwd () will allocate rs_data(ista)%ind_intp_obs(nobs_max) )
          IF (loutradwind) THEN
            ! .. will allocate rs_data(ista)%radwind_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs or rs_data(ista)%ind_intp_obs_vr(nobs_max)
            CALL read_field_obs_dwd_h5 (ista,i_vrad,irep, &
                 'VRAD',rs_data(ista)%radwind_obs, rs_data(ista)%ind_intp_obs_vr, &
                 rs_meta(ista)%lobs_avail(irep,i_vrad))
!!$ UB: quick fix for non-existent quality flags in DWD data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualvrad) = rs_meta(ista)%lobs_avail(irep,i_vrad)
          END IF
          IF (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) THEN
            ! .. will allocate rs_data(ista)%zh_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_dwd_h5 (ista,i_dbzh,irep, &
                 'DBZH',rs_data(ista)%zh_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_dbzh))
!!$ UB: quick fix for non-existent quality flags in DWD data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualdbzh) = rs_meta(ista)%lobs_avail(irep,i_dbzh)
          END IF
          IF (loutpolstd .OR. loutpolall) THEN
            ! .. will allocate rs_data(ista)%zdr_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_dwd_h5 (ista,i_zdr,irep, &
                 'ZDR',rs_data(ista)%zdr_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_zdr))
            CALL read_field_obs_dwd_h5 (ista,i_kdp,irep, &
                 'KDP',rs_data(ista)%kdp_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_kdp))
            CALL read_field_obs_dwd_h5 (ista,i_phidp,irep, &
                 'UPHIDP',rs_data(ista)%phidp_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_phidp))
            CALL read_field_obs_dwd_h5 (ista,i_rhv,irep, &
                 'RHOHV',rs_data(ista)%rhv_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_rhv))
          END IF
          
          file_identified = .TRUE.

        END IF
#endif

        IF ( TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_cdfin) .AND. &
             ANY( rs_meta(ista)%icountry == (/ i_dwd, i_meteoswiss, i_arpasim, i_belgium, i_denmark, &
                                               i_france, i_poland, i_czech, i_netherlands, i_slovakia/) ) ) THEN

          ! NOTE: if data are from an OSE / OSSE, rs_meta(ista)%icountry will not necessarily be i_dwd,
          !       because get_metadata_from_cdfin() has set it according to the background list for the respective
          !       station_id!

          ! .. This is either the old DWD cdfin-format "cdfin_<vr,qv,z,qz>_<ii>" or the
          !     new format "cdfin_z_id-010908_201307282200_201307282255_volscan"     or
          !                "cdfin_z_id-010908_201307282200_201307282255_precipscan"

          !.. allocate field for cross-checking whether the 4 data
          !   sets to be read in below are from the same
          !   scan:
          rs_data(ista)%check(:) = -99

          ! .. Reading of data:
          !      ( read_field_obs_dwd () will allocate rs_data(ista)%ind_intp_obs(nobs_max) )
          IF (loutradwind) THEN
            ! .. will allocate rs_data(ista)%radwind_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_vr(nobs_max)
            CALL read_field_obs_dwd (ista,i_vrad,irep, &
                 'MDMVR',rs_data(ista)%radwind_obs, rs_data(ista)%ind_intp_obs_vr, &
                 rs_meta(ista)%lobs_avail(irep,i_vrad),&
                 missing_vr)
            IF (lqc_flag) THEN
              ! .. will allocate rs_data(ista)%radwind_obs_q(nobs_max) and, if needed,
              !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_qv(nobs_max)
              CALL read_field_obs_dwd (ista,i_qualvrad,irep, &
                   '021218',rs_data(ista)%radwind_obs_q,  rs_data(ista)%ind_intp_obs_qv, &
                   rs_meta(ista)%lobs_avail(irep,i_qualvrad),&
                   missing_qv)
            ELSE
!!$ UB: quick fix, but not clean!
              rs_meta(ista)%lobs_avail(irep,i_qualvrad) = rs_meta(ista)%lobs_avail(irep,i_vrad)
            ENDIF
          END IF
          IF (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) THEN
            ! .. will allocate rs_data(ista)%zh_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_dwd (ista,i_dbzh,irep, &
                 'MHORRE0',rs_data(ista)%zh_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_dbzh),&
                 missing_z)
            IF (lqc_flag) THEN
              ! .. will allocate rs_data(ista)%zh_radar_obs_q(nobs_max) and, if needed,
              !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_qz(nobs_max)
              CALL read_field_obs_dwd (ista,i_qualdbzh,irep, &
                   '021218',rs_data(ista)%zh_radar_obs_q, rs_data(ista)%ind_intp_obs_qz, &
                   rs_meta(ista)%lobs_avail(irep,i_qualdbzh),&
                   missing_qz)
            ELSE
!!$ UB: quick fix, but not clean!
              rs_meta(ista)%lobs_avail(irep,i_qualdbzh) = rs_meta(ista)%lobs_avail(irep,i_dbzh)
            ENDIF
          END IF

          !.. Now map the missing values and "correct 0's" correctly:
          IF ((loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) .AND. rs_meta(ista)%lobs_avail(irep,i_dbzh)) THEN
            IF (rs_meta(ista)%obs_cdate(irep) < '20130625082500') THEN
              ! At that time, there was not really a distinction between
              !  "correct 0's" and "erroneous observations", both coded as "_fill_value". Because the
              !  former are assumed to happen much more than the latter,
              !  we assign "correct 0's" to them:
              WHERE (rs_data(ista)%zh_radar_obs == missing_z)
                rs_data(ista)%zh_radar_obs = zero_value
              END WHERE
              IF ( lqc_flag .AND. rs_meta(ista)%lobs_avail(irep,i_qualdbzh) ) THEN
                WHERE (rs_data(ista)%zh_radar_obs_q == missing_qz)
                  rs_data(ista)%zh_radar_obs_q = missing_obs
                END WHERE
              END IF
            ELSE
              ! now "_fill_values" are assigned to measurements in blocked sectors
              !  and/or values beyond the range (different elevations can have different ranges,
              !  and the data are padded to the largest range among the elevations).
              !  -32 dBZ are "correct 0's" or "invalid data recognized by the signal processor".
              ! This information is a little better than the above, but still not perfect.
              !   We assign "_fill_values" to the missing value -999.99 and "correct 0's"
              !   to the -32:
              WHERE (rs_data(ista)%zh_radar_obs == missing_z)
                rs_data(ista)%zh_radar_obs = missing_obs
              ELSEWHERE ( rs_data(ista)%zh_radar_obs < -31.99_dp .AND. rs_data(ista)%zh_radar_obs > miss_threshold )
                rs_data(ista)%zh_radar_obs = zero_value
              END WHERE
              IF ( lqc_flag .AND. rs_meta(ista)%lobs_avail(irep,i_qualdbzh) ) THEN
                WHERE (rs_data(ista)%zh_radar_obs_q == missing_qz)
                  rs_data(ista)%zh_radar_obs_q = missing_obs
                END WHERE
              END IF
            END IF
          END IF

          IF (loutradwind .AND. rs_meta(ista)%lobs_avail(irep,i_vrad)) THEN
            WHERE (rs_data(ista)%radwind_obs == missing_vr)
              rs_data(ista)%radwind_obs = missing_obs
            END WHERE
            IF ( lqc_flag .AND. rs_meta(ista)%lobs_avail(irep,i_qualvrad) ) THEN
              WHERE (rs_data(ista)%radwind_obs_q == missing_qv)
                rs_data(ista)%radwind_obs_q = missing_obs
              END WHERE
            END IF
          END IF

          file_identified = .TRUE.

        END IF

#ifdef HDF5_RADAR_INPUT

        ! .. OPERA files which are distributed to OPERA from NWS's and which have filename-convention
        !     "T_PAGA52_C_LSSW_20220104101500.hdf"
        !     "T_PAGB52_C_LSSW_20220104101500.hdf"
        !     "T_PAGC52_C_LSSW_20220104101500.hdf"
        !     "T_PAGD52_C_LSSW_20220104101500.hdf"
        !     "T_PAGE52_C_LSSW_20220104101500.hdf"
          
        !     "T_PAGA41_C_LSSW_20220104102500.hdf"
        !     "T_PAGB41_C_LSSW_20220104102500.hdf"
        !     "T_PAGC41_C_LSSW_20220104102500.hdf"
        !     "T_PAGD41_C_LSSW_20220104102500.hdf"
        !     "T_PAGE41_C_LSSW_20220104102500.hdf"

        ! .. Swiss HDF5 radar files (only their basenames):
        !     "MLA1520500007U.001.V.h5"
        !     "MLA1520500007U.001.Z.h5"
        !     "MLA1520500007U.002.V.h5"
        !     "MLA1520500007U.002.Z.h5"

        IF ( ( rs_meta(ista)%icountry == i_dwd        .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_smss) ) .OR. &
             ( rs_meta(ista)%icountry == i_meteoswiss .AND. &
             &     ( TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_mmss) .OR. & 
             &       TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_smss)  .OR. &
             &       TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_mmss) ) ) .OR. &
             ( rs_meta(ista)%icountry == i_france     .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_mmss) ) .OR. &
             ( rs_meta(ista)%icountry == i_belgium    .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
             ( rs_meta(ista)%icountry == i_poland     .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
             ( rs_meta(ista)%icountry == i_czech      .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
             ( rs_meta(ista)%icountry == i_slovakia   .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
             ( rs_meta(ista)%icountry == i_dwd   .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
             ( rs_meta(ista)%icountry == i_netherlands.AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_mmms) ) .OR. &
             ( rs_meta(ista)%icountry == i_denmark    .AND. &
             &     TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_2_opera_mmms) )      &
             ) THEN

          rs_data(ista)%check(:) = -99

          ! .. Reading of data:
          !      ( read_field_obs_opera_h5 () will allocate rs_data(ista)%ind_intp_obs(nobs_max) )
          IF (loutradwind) THEN
            ! .. will allocate rs_data(ista)%radwind_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_vr(nobs_max)
            CALL read_field_obs_opera_h5 (ista,i_vrad,irep, &
                 'VRAD',rs_data(ista)%radwind_obs, rs_data(ista)%ind_intp_obs_vr, &
                 rs_meta(ista)%lobs_avail(irep,i_vrad))
!!$ UB: quick fix for non-existent quality flags in OPERA data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualvrad) = rs_meta(ista)%lobs_avail(irep,i_vrad)
          END IF
          IF (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) THEN
            ! .. will allocate rs_data(ista)%zh_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_opera_h5 (ista,i_dbzh,irep, &
                 'DBZH',rs_data(ista)%zh_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_dbzh))
!!$ UB: quick fix for non-existent quality flags in OPERA data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualdbzh) = rs_meta(ista)%lobs_avail(irep,i_dbzh)
          END IF
          
          file_identified = .TRUE.

        END IF
#endif

!!$#ifdef HDF5_RADAR_INPUT
!!$
!!$        ! .. Swiss HDF5 radar files (only their basenames):
!!$        !     "PLA1520500007U.001.V.h5"
!!$        !     "PLA1520500007U.001.Z.h5"
!!$        !     "PLA1520500007U.002.V.h5"
!!$        !     "PLA1520500007U.002.Z.h5"
!!$
!!$        IF ( rs_meta(ista)%icountry == i_meteoswiss .AND. &
!!$             TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_smss) ) THEN
!!$
!!$          rs_data(ista)%check(:) = -99
!!$
!!$          ! .. Reading of data:
!!$          !      ( read_field_obs_mch () will allocate rs_data(ista)%ind_intp_obs(nobs_max) )
!!$          IF (loutradwind) THEN
!!$            ! .. will allocate rs_data(ista)%radwind_obs(nobs_max) and, if needed,
!!$            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_vr(nobs_max)
!!$            CALL read_field_obs_mch (ista,i_vrad,irep, &
!!$                 'VRAD',rs_data(ista)%radwind_obs, rs_data(ista)%ind_intp_obs_vr, &
!!$                 rs_meta(ista)%lobs_avail(irep,i_vrad))
!!$!!$ UB: quick fix for non-existent quality flags in MCH data, but not clean!
!!$            rs_meta(ista)%lobs_avail(irep,i_qualvrad) = rs_meta(ista)%lobs_avail(irep,i_vrad)
!!$          END IF
!!$          IF (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) THEN
!!$            ! .. will allocate rs_data(ista)%zh_radar_obs(nobs_max) and, if needed,
!!$            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_z(nobs_max)
!!$            CALL read_field_obs_mch (ista,i_dbzh,irep, &
!!$                 'DBZH',rs_data(ista)%zh_radar_obs, rs_data(ista)%ind_intp_obs_z, &
!!$                 rs_meta(ista)%lobs_avail(irep,i_dbzh))
!!$!!$ UB: quick fix for non-existent quality flags in MCH data, but not clean!
!!$            rs_meta(ista)%lobs_avail(irep,i_qualdbzh) = rs_meta(ista)%lobs_avail(irep,i_dbzh)
!!$          END IF
!!$          
!!$          file_identified = .TRUE.
!!$
!!$        END IF
!!$#endif

#ifdef HDF5_RADAR_INPUT

        
        IF ( rs_meta(ista)%icountry == i_arpasim .AND. &
             TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_mmms) ) THEN

          ! .. ODIM HDF5 of Italy:

          rs_data(ista)%check(:) = -99

          ! .. Reading of data:
          !      ( read_field_obs_italy () will allocate rs_data(ista)%ind_intp_obs(nobs_max) )
          IF (loutradwind) THEN
            ! .. will allocate rs_data(ista)%radwind_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_vr(nobs_max)
            CALL read_field_obs_italy (ista,i_vrad,irep, &
                 'VRAD',rs_data(ista)%radwind_obs, rs_data(ista)%ind_intp_obs_vr, &
                 rs_meta(ista)%lobs_avail(irep,i_vrad))
!!$ UB: quick fix for non-existent quality flags in MCH data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualvrad) = rs_meta(ista)%lobs_avail(irep,i_vrad)
          END IF
          IF (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) THEN
            ! .. will allocate rs_data(ista)%zh_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs) or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_italy (ista,i_dbzh,irep, &
                 'DBZH',rs_data(ista)%zh_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_dbzh))
!!$ UB: quick fix for non-existent quality flags in MCH data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualdbzh) = rs_meta(ista)%lobs_avail(irep,i_dbzh)
          END IF

          ! THOMAS modification:
          ! Radial wind observations associated to "nodata" or "undetect" reflectivities
          ! are set to "missing_obs"
          IF ((loutradwind) .AND. (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0))) THEN
              ! Correct missing_obs
              WHERE (rs_data(ista)%zh_radar_obs == missing_obs .OR. rs_data(ista)%zh_radar_obs == zero_value)
                rs_data(ista)%radwind_obs = missing_obs
              END WHERE

              ! THOMAS modification: Nyquist velocity is set as the maximum between the 
              ! maximum (in absolute value) of the read values (for the specific elevation)
              ! and a default value (7 m/s) considered as a reasonable minimum. This is 
              ! necessary since, at present, it is not possible to retrieve Nyquist 
              ! velocity for all the Italian radars.
              nel     = rs_meta(ista)%nel
              nobs_el = rs_meta(ista)%naz * rs_meta(ista)%nra
              DO eleind = 1, nel
                start_ind = 1 + (eleind-1) * nobs_el
                end_ind   =      eleind    * nobs_el
                rs_meta(ista)%ext_nyq(eleind, 1) = MAX( MAXVAL( ABS(rs_data(ista)%radwind_obs(start_ind:end_ind)), &
                     mask = ABS(rs_data(ista)%radwind_obs(start_ind:end_ind)) < ABS(missing_obs) ), 7.0_dp )
                WRITE (*,'(a,i6.6,a,f0.1,a,f0.2,a)') 'INFO '//TRIM(yzroutine)//': Nyquist velocity for radar ', &
                         rs_meta(ista)%station_id, ' at elevation ', rs_meta(ista)%el_arr(eleind), ' deg is set to ',&
                         rs_meta(ista)%ext_nyq(eleind, 1), ' m/s'
              END DO
          END IF

          file_identified = .TRUE.

        END IF
#endif

#ifdef HDF5_RADAR_INPUT
        IF ( rs_meta(ista)%icountry == i_dwd .AND. &
             ( TRIM(rs_meta(ista)%obsfile_format(irep)) == TRIM(c_format_h5_native_mmms) ) ) THEN

          ! .. KIT C-band files of filename-convention:
          !     "scan-sidpol-120km-14_20001_20230712123503_00.h5"
          !     "scan-sidpol-120km-14_20001_20230712124001_00.h5"
          !
          ! The station 20001 is listed under i_dwd and produces mmms files, which the original DWD-radars don't

          rs_data(ista)%check(:) = -99

          ! .. Reading of data:
          !      ( read_field_obs_dwd () will allocate rs_data(ista)%ind_intp_obs(nobs_max) )
          IF (loutradwind) THEN
            ! .. will allocate rs_data(ista)%radwind_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs or rs_data(ista)%ind_intp_obs_vr(nobs_max)
            CALL read_field_obs_kit_h5 (ista,i_vrad,irep, &
                 'VRAD',rs_data(ista)%radwind_obs, rs_data(ista)%ind_intp_obs_vr, &
                 rs_meta(ista)%lobs_avail(irep,i_vrad))
!!$ UB: quick fix for non-existent quality flags in DWD data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualvrad) = rs_meta(ista)%lobs_avail(irep,i_vrad)
          END IF
          IF (loutdbz .OR. (loutradwind .AND. itype_obserr_vr > 0)) THEN
            ! .. will allocate rs_data(ista)%zh_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_kit_h5 (ista,i_dbzh,irep, &
                 'DBZH',rs_data(ista)%zh_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_dbzh))
!!$ UB: quick fix for non-existent quality flags in DWD data, but not clean!
            rs_meta(ista)%lobs_avail(irep,i_qualdbzh) = rs_meta(ista)%lobs_avail(irep,i_dbzh)
          END IF
          IF (loutpolstd .OR. loutpolall) THEN
            ! .. will allocate rs_data(ista)%zdr_radar_obs(nobs_max) and, if needed,
            !    alternatively rs_data(ista)%ind_intp_obs or rs_data(ista)%ind_intp_obs_z(nobs_max)
            CALL read_field_obs_kit_h5 (ista,i_zdr,irep, &
                 'ZDR',rs_data(ista)%zdr_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_zdr))
            CALL read_field_obs_kit_h5 (ista,i_kdp,irep, &
                 'KDP',rs_data(ista)%kdp_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_kdp))
            CALL read_field_obs_kit_h5 (ista,i_phidp,irep, &
                 'PHIDP',rs_data(ista)%phidp_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_phidp))
            CALL read_field_obs_kit_h5 (ista,i_rhv,irep, &
                 'RHOHV',rs_data(ista)%rhv_radar_obs, rs_data(ista)%ind_intp_obs_z, &
                 rs_meta(ista)%lobs_avail(irep,i_rhv))
          END IF
          
          file_identified = .TRUE.

        END IF
#endif

        ! .. Add readers for other formats in the future!

        IF (.NOT. file_identified) THEN
          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': obs data file '//TRIM(rs_meta(ista)%obsfile(irep,ikind))// &
               ' is not of known type! Reading not possible!'
        ELSE

          ! .. Sort the list of present elevations in ascending order. This is not guaranteed by the above readers:
          IF (rs_meta(ista)%nel_present > 0) THEN
            CALL bubblesort ( rs_meta(ista)%ind_ele_present(1:rs_meta(ista)%nel_present) )
          END IF

          ! .. Optional range coarsening (re-defines the data vectors and rs_data%nobs_obs):
          IF (rs_meta(ista)%n_aggr_ra_obs > 1) THEN
            CALL range_coarsening ( ista, irep, ierr )
            ! .. Nothing is done with the error status ierr so far
          END IF

        END IF

      END IF format_choice

      IF (ALL(.NOT.rs_meta(ista)%lobs_avail(irep,:))) THEN
        ! No obs have been read. We re-set the list of present elevations to match the nominal strategy,
        !  in order to enable later output of simulated composits of all nominal elevations.
        rs_meta(ista)%nel_present = rs_meta(ista)%nel
        rs_meta(ista)%ind_ele_present(:) = missval_int
        rs_meta(ista)%ind_ele_present(1:rs_meta(ista)%nel_present) = (/ (i, i=1, rs_meta(ista)%nel_present) /)
      END IF

    ELSE IF (action == 'cleanup') THEN

      ! clean and deallocate memory from pointers in rs_data(ista):
      CALL cleanup_read_obs_rad_ista ()

    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine),' ',action,' on proc ', my_radar_id


  CONTAINS

    ! clean and deallocate memory from pointers in rs_data(ista):

    SUBROUTINE cleanup_read_obs_rad_ista

      IF (ASSOCIATED (rs_data(ista)%radwind_obs))     DEALLOCATE(rs_data(ista)%radwind_obs)
      IF (ASSOCIATED (rs_data(ista)%radwind_obs_q))   DEALLOCATE(rs_data(ista)%radwind_obs_q)
      IF (ASSOCIATED (rs_data(ista)%zh_radar_obs))    DEALLOCATE(rs_data(ista)%zh_radar_obs)
      IF (ASSOCIATED (rs_data(ista)%zh_radar_obs_q))  DEALLOCATE(rs_data(ista)%zh_radar_obs_q)
      IF (ASSOCIATED (rs_data(ista)%zdr_radar_obs))   DEALLOCATE(rs_data(ista)%zdr_radar_obs)
      IF (ASSOCIATED (rs_data(ista)%kdp_radar_obs))   DEALLOCATE(rs_data(ista)%kdp_radar_obs)
      IF (ASSOCIATED (rs_data(ista)%phidp_radar_obs)) DEALLOCATE(rs_data(ista)%phidp_radar_obs)
      IF (ASSOCIATED (rs_data(ista)%rhv_radar_obs))   DEALLOCATE(rs_data(ista)%rhv_radar_obs)
      IF (lequal_azi_alldatasets) THEN
        ! .. This is only allocated if lequal_azi_alldatasets = .true.:
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs))       DEALLOCATE(rs_data(ista)%ind_intp_obs)
      ELSE
        ! .. The next are allocated only if lequal_azi_alldatasets = .false.
        !     and only point to rs_data(ista)%ind_intp_obs otherwise:
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_vr))    DEALLOCATE(rs_data(ista)%ind_intp_obs_vr)
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_qv))    DEALLOCATE(rs_data(ista)%ind_intp_obs_qv)
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_z))     DEALLOCATE(rs_data(ista)%ind_intp_obs_z)
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_qz))    DEALLOCATE(rs_data(ista)%ind_intp_obs_qz)
      END IF
      NULLIFY(rs_data(ista)%radwind_obs)
      NULLIFY(rs_data(ista)%radwind_obs_q)
      NULLIFY(rs_data(ista)%zh_radar_obs)
      NULLIFY(rs_data(ista)%zh_radar_obs_q)
      NULLIFY(rs_data(ista)%zdr_radar_obs)
      NULLIFY(rs_data(ista)%kdp_radar_obs)
      NULLIFY(rs_data(ista)%phidp_radar_obs)
      NULLIFY(rs_data(ista)%rhv_radar_obs)
      NULLIFY(rs_data(ista)%ind_intp_obs)
      NULLIFY(rs_data(ista)%ind_intp_obs_vr)
      NULLIFY(rs_data(ista)%ind_intp_obs_qv)
      NULLIFY(rs_data(ista)%ind_intp_obs_z)
      NULLIFY(rs_data(ista)%ind_intp_obs_qz)

    END SUBROUTINE cleanup_read_obs_rad_ista

  END SUBROUTINE read_obs_rad

  SUBROUTINE get_posind_rounded_fullrev (info_in_case_of_error, az_volscan, el_volscan, naz_ncdf, nele_ncdf,  &
                                missing_az, az_start, az_inc, &
                                el_arr, nele,    &
                                azind, eleind, naz)

    !------------------------------------------------------------------------------
    !
    ! Description: For azimut- and elevation-data of a DWD volume scan,
    !              determine the usable data sector (constant elevation,
    !              continuous valid azimut range). Then, compute the
    !              "regular" azimut- and elevation indices valid for
    !              approximating the "true" ones by a regular azimut- grid
    !              which spans one full antenna revolution (360 degrees)
    !              with a constant increment of "az_inc".
    !
    ! Method:
    !
    !   Round "middle" azimut value of the sweeps to nearest multiple of az_inc and successively
    !   increment/decrement the azimut index by 1 from neighbour to neighbour.
    !   The direction of antenna rotation is correctly taken into account.
    !
    !   We take into account the fact that the azimut values given in
    !   the radar data files denote the end of the azimutal pulse averaging interval
    !   in the direction of antenna rotation. This is achieved by "cheating" with
    !   the nominal az_start, in that it is increased/decreased by 0.5*az_inc in the
    !   antenna rotation direction ("az_start_mod" below).
    !
    !   As a result, the subroutine returns regular azimut indices i in a nearest-neighbour
    !   sense for each ray of the sweep.
    !   From i, the nominal azimut value can be computed by:
    !
    !   az(i) = az_start + (i-1)*az_inc   ,  i = 1 ... naz (full antenna revolution)
    !
    !      where    naz = NINT(360 / az_inc)
    !
    !   and az_start denotes the center of the first azimut pulse averaging interval (sampling interval).
    !
    ! Outputs:
    !
    !   - naz : number of nominal azimuts for full antenna revolution
    !   - azind(1:naz_ncdf, 1:nele_ncdf): nominal azimut index (equals "i" above)
    !                                     for each radar bin in each elevation sweep
    !   - eleind(1:nele_ncdf): nominal elevation index for each elevation
    !
    !   - "missing" azimuts and elevations in the data file are flagged by azind(i) = missval_int
    !
    !
    ! POSSIBLE IMPROVEMENT: Interpolation to regular azimut grid!
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    CHARACTER (len=*), INTENT(in) ::  info_in_case_of_error ! string that should contain infos to identify time and radar station for the warning/error messages
    INTEGER, INTENT(in) ::  naz_ncdf, nele_ncdf  ! field dimensions azim, elev
    INTEGER, INTENT(in) ::  nele                 ! vector length of nominal elevations el_arr

    REAL (KIND=dp), INTENT(in) :: az_volscan(naz_ncdf,nele_ncdf) ! array of azimuts for each ray of a volume scan     [deg]
    REAL (KIND=dp), INTENT(in) :: el_volscan(naz_ncdf,nele_ncdf) ! array of elevations for each ray of a volume scan  [deg]
    REAL (KIND=dp), INTENT(in) ::  az_start, az_inc              ! nominal azimut start and increment [deg]
    REAL (KIND=dp), INTENT(in) ::  el_arr(nele)                  ! nominal elevation vector           [deg]
    REAL (KIND=dp), INTENT(in) ::  missing_az                    ! missing value in azimut array      [deg]

    ! .. INTENT(OUT):
    INTEGER, INTENT(out) ::  naz                 ! number of nominal azimuts for full antenna revolution
    INTEGER, INTENT(out) ::  azind(naz_ncdf,nele_ncdf) ! nominal azimut index of each ray
    INTEGER, INTENT(out) ::  eleind(nele_ncdf)   ! nominal elevation index of each elevation


    !------------------------------------------------------------------------------
    !
    ! .. Local variables:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_posind_rounded_fullrev'
    INTEGER              :: i, o, zazmin, zazmax, znaz_ncdf, rotsign, znaz, &
                            zazind_start, azmin_el, azmax_el, azmid
    REAL (KIND=dp)       :: reldiffazi(10), elemean, az_start_mod, &
                            azi_offset, azinomdiff(naz_ncdf), tmpdiff
    REAL (KIND=dp), PARAMETER :: eps   = 1e-20_dp
    REAL (KIND=dp), PARAMETER :: eps01 = 0.15_dp

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !.. number of azimuts for full antenna revolution:
    naz = NINT(360.0_dp / az_inc)

    !.. Initializations with missing values:
    azind  = missval_int
    eleind = missval_int

    !.. Determine the regular azimuts and their nominal indices for each elevation
    !   along with the nominal elevation index of the elevation:
    eleloop: DO o = 1, nele_ncdf

      !.. Clip the azimut range to its non-missing range (search from the back):
      DO i = naz_ncdf, 1, -1
        IF (az_volscan(i,o) /= missing_az) THEN
          znaz_ncdf = i
          EXIT
        END IF
      END DO

      !.. Antenna rotation direction +1 or -1, determined from the first few azimuts:
      zazmax = MIN(10,znaz_ncdf-1)
      reldiffazi = 0.0_dp
      reldiffazi(1:zazmax) = az_volscan(2:zazmax+1,o) - az_volscan(1:zazmax,o)
      reldiffazi(1:zazmax) = reldiffazi(1:zazmax) / ( ABS(reldiffazi(1:zazmax)) + eps )
      rotsign = NINT(SIGN(1.0_dp, SUM(reldiffazi)))

      !.. The azimut values in the radar files denote the start of the azimut averaging interval
      !   in antenna rotation direction. To ensure that the below determined
      !   nearest regular azimuts represent the center of the averaging intervals,
      !   cheat here by shifting the "regular" az_start accordingly by 0.5*az_inc:
      az_start_mod = az_start + rotsign*az_inc*0.5

      !.. Determine mean elevation of sweep from mid section of azimut interval
      !   (should be nearly constant over the sweep):
      zazmin = MAX( NINT(znaz_ncdf*0.5)-5, 1) ! Mid azi of ele sweep - 5
      zazmax = MIN( zazmin + 10, znaz_ncdf)                       ! Mid azi of ele sweep + 5
      elemean = SUM(el_volscan(zazmin:zazmax,o)) / (zazmax-zazmin+1.0)


      !.. Determine which of the nominal elevations is nearest to the actual elemean and
      !   store the corresponding elevation index in eleind(o):
      eleind(o) = 1
      tmpdiff = ABS(elemean-el_arr(eleind(o)))
!CDIR NODEP
      DO i = 1, nele
        IF ( ABS(elemean-el_arr(i)) < tmpdiff ) THEN
          eleind(o) = i
          tmpdiff = ABS(elemean-el_arr(i))
        END IF
      END DO

      IF ( tmpdiff > eps01 ) THEN
        eleind(o) = missval_int
        WRITE (*,'(a,i5,1x,f0.4,1x,a)') 'WARNING '//TRIM(yzroutine)//': an elevation from a' // &
             ' NetCDF file does not have a correct nominal elevation angle! eleind, elemean, data_info = ', &
             o, elemean, TRIM(info_in_case_of_error)
      END IF

      !.. Find the azimut part with constant elevation by searching from both sides
      !   of the azimut sector towards its center. This assumes that only
      !   the margins of the azimutal sweep deviate:
      azmin_el = missval_int
      azmax_el = missval_int
      DO i = 1, znaz_ncdf
        IF ( ABS(el_volscan(i,o)-elemean) <= eps01 ) THEN
          azmin_el = i
          EXIT
        END IF
      END DO
      DO i = znaz_ncdf, 1, -1
        IF ( ABS(el_volscan(i,o)-elemean) <= eps01 ) THEN
          azmax_el = i
          EXIT
        END IF
      END DO

      IF (azmin_el == missval_int .OR. azmax_el == missval_int) THEN
        WRITE (*,'(a,2i5,1x,f0.4,1x,a)') 'WARNING '//TRIM(yzroutine)//': a sweep from a NetCDF file' // &
             ' does not have azimuts at constant elevation!  naz, eleind, elemean, data_info = ', znaz_ncdf, o, elemean, &
             TRIM(info_in_case_of_error)
        CYCLE eleloop
      END IF

      !.. number of azimuts with constant elevation:
      znaz = azmax_el - azmin_el + 1

      IF ( znaz < 20 ) THEN
        WRITE (*,'(a,3i5,1x,a)') 'WARNING '//TRIM(yzroutine)//': a sweep from a NetCDF file' // &
             ' has less than 20 azimuts at constant elevation! ',azmin_el,azmax_el,o, TRIM(info_in_case_of_error)
      ELSE IF ( znaz < naz ) THEN
        WRITE (*,'(a,3i5,1x,a)') 'WARNING '//TRIM(yzroutine)//': a sweep from a NetCDF file' // &
             ' does not cover one complete antenna revolution! ',azmin_el,azmax_el,o, TRIM(info_in_case_of_error)
      END IF

      !.. Rounding of the nominal azimut index of the middle regular azimut
      !   and incrementing the indices for the next azimuts, starting from the
      !   first regular azimut. This assumes that
      !   the azimuts of the sweep are continuous. The direction of
      !   antenna rotation is correctly taken into account, so that
      !   the nominal azimut index is increasing with increasing azimut:

      azmid = NINT(0.5*(azmin_el+azmax_el))

      zazind_start = NINT( MODULO(az_volscan(azmid,o)-(azmid-azmin_el)*az_inc - &
                                  az_start_mod, 360.0_dp ) / az_inc ) + 1
      DO i = azmin_el, MIN(azmax_el, azmin_el+naz-1)
        azind(i,o) = MODULO( zazind_start+rotsign*i-azmin_el - 1, naz) + 1
      END DO

      !.. Check azimuts against nominal values, taking into account a rounding shift of
      !   the azimut determined from the mid point of the sweep:

      azi_offset = az_volscan(azmid,o) - ( az_start_mod + ( azind(azmid,o) - 1 ) * az_inc )

      IF ( azi_offset >= 360.0-2*az_inc ) THEN
        azi_offset = azi_offset - 360.0
      ELSE IF (azi_offset <= -360.0+2*az_inc) THEN
        azi_offset = azi_offset + 360.0
      END IF

!CDIR NODEP
      DO i = 1, MIN(znaz,naz)

        azinomdiff(i) = MODULO(az_start_mod + (azind(azmin_el+i-1,o)-1)*az_inc, 360.0_dp) - &
                        az_volscan(azmin_el+i-1,o) + azi_offset

        IF ( azinomdiff(i) >= 360.0-az_inc ) THEN
          azinomdiff(i) = azinomdiff(i) - 360.0
        ELSE IF (azinomdiff(i) <= -360.0+az_inc) THEN
          azinomdiff(i) = azinomdiff(i) + 360.0
        END IF
        azinomdiff(i) = ABS(azinomdiff(i))

      END DO

      IF ( MAXVAL(azinomdiff(1:MIN(znaz,naz))) > 0.15*ABS(az_inc) ) THEN
        WRITE (*,'(a,i5,1x,a)') 'WARNING '//TRIM(yzroutine)//': not all azimuts from NetCDF file' // &
             ' agree with their nominal value!  eleind, data_info = ',o, TRIM(info_in_case_of_error)
        DO i=1, MIN(znaz,naz)
          IF (azinomdiff(i) > 0.15*ABS(az_inc)) THEN
            WRITE (*,'(a,i4,4f10.4)') '  az_i = ', i, azinomdiff(i), az_volscan(azmin_el+i-1,o), &
                 az_volscan(azmin_el,o), azi_offset
          END IF
        END DO
      END IF

    END DO eleloop

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_posind_rounded_fullrev

  SUBROUTINE range_coarsening ( ista, irep, ierr )

    !------------------------------------------------------------------------------
    !
    ! Description: Range coarsening by aggregation of input observations to reduce the data
    !              amount. This goes together with a coarser range of the simulated
    !              synthetic data and can be useful if the model resolution is much
    !              coarser than the radar range resolution to save computing time
    !              and memory. Also the computation of composites and fdbk-files
    !              as well as the output of volume scan data will be faster without
    !              much loss of accuracy in the context of the NWP model.
    !
    ! Method:      Coarse range bins are formed by aggregating an integer multiple
    !              of original range bins. Only integer multiples are supported.
    !              The method of aggregation depens on the parameter:
    !
    !              ZH:     Average value in linear ZH space (not dBZ)
    !              VRAD:   Reflectivity-weighted average value (linear ZH is the weight)
    !              ZDR:    Compute ZV from ZH and ZDR, aggregate both separately,
    !                      derive ZDR_aggr as ZH_aggr / ZV_aggr, convert ZDR_aggr back
    !                      to dB space.
    !              RHOHV:  Reflectivity-weighted average value similar to VRAD
    !              KDP:    Simple average value
    !              PHIDP:  Outermost original value within the coarse range bin.
    !                      This should be consistent to KDP if both were consistent
    !                      in the original data.
    !
    !              Missing values (missing_obs) are excluded from the aggregation.
    !
    !              The new range center of each coarsened bin is defined similar
    !              to the original range bins, but with a coarser resolution:
    !
    !                  r_coarse(i) = i * ra_inc_coarse,  i = 1, ... , nra_coarse
    !
    !
    ! Prerequisites:    The meta data defining the range bins should have been
    !                   prepared properly (is done in radar_namelist_read.f90):
    !
    !       - rs_meta(ista)%nra_obs    should contain the original number of range bins
    !       - rs_meta(ista)%nra        should contain the coarsened number of range bins
    !       - rs_meta(ista)%ra_inc_obs should contain the original range resolution
    !       - rs_meta(ista)%ra_inc     should contain the coarse range resolution
    !       - rs_meta(ista)%n_aggr_ra_inc should contain the number of original
    !                                     range bins to aggregate
    !
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with INTENT(IN):

    INTEGER, INTENT(in)         :: ista, irep

    ! .. INTENT(OUT):
    INTEGER, INTENT(out)        :: ierr

    !------------------------------------------------------------------------------
    !
    ! .. Local variables:

    INTEGER                     :: iobs, i, ii, k, kk, m, mm, n, nn, o, nra, nra_obs, naz, nel, nobs_max, nobs_max_obs, &
                                   nobs_obs(ndatakind), n_aggr
    REAL(kind=dp), ALLOCATABLE  :: w_aggr(:), w_sum(:), var_aggr(:,:), var_aggr_wgt(:,:), var_sum(:), tmpzh(:,:), tmpzv(:,:)
    REAL(kind=dp)               :: zloc, indataloc
    INTEGER                     :: ira_offset
    
    CHARACTER(len=*), PARAMETER :: yzroutine = 'range_coarsening'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ierr = 0
    
    n_aggr = rs_meta(ista)%n_aggr_ra_obs

    IF (n_aggr > 1) THEN

      ! .. Get scan dimensions:
      nra_obs = rs_meta(ista)%nra_obs ! orig number of range bins
      nra     = rs_meta(ista)%nra     ! coarsened number of range bins
      naz     = rs_meta(ista)%naz
      nel     = rs_meta(ista)%nel
      nobs_max_obs = nra_obs * naz * nel
      nobs_max     = nra     * naz * nel
      nobs_obs(:) = rs_data(ista)%nobs_obs(:)

      ! .. Check consistency of nobs_obs among the different data sets. If not fullfilled,
      !    some coarsened data where the coarsening process depends on others
      !    (e.g. coarsened radial wind is reflectivity weighted) will have to be set to
      !    missing_value. Give a warning in that case:
      IF (ANY(nobs_obs /= MAXVAL(nobs_obs) .AND. nobs_obs /= 0)) THEN
        WRITE (*,'(a,i6.6,a,'//cndatakind//'(1x,i9),/,a,f0.2,a)') 'WARNING '//TRIM(yzroutine)// &
             ' station=', rs_meta(ista)%station_id, &
             ': differing nobs_obs among data sets:', nobs_obs(:), &
             '         Some coarsened data with cross-dependencies to others ' // &
             '(e.g. radial wind is reflectivity weighted) will have to be set to ',missing_obs,'!'
      END IF

      ! .. Allocate work arrays:
      ALLOCATE(var_aggr(nobs_max, n_aggr))
      ALLOCATE(var_aggr_wgt(nobs_max, n_aggr))
      ALLOCATE(var_sum(nobs_max))
      ALLOCATE(w_sum(nobs_max))

      ! .. Prepare weights for range coarsening:
      ALLOCATE (w_aggr(n_aggr))
      IF (MOD(n_aggr,2) == 0) THEN
        ! Even number of aggregated original range bins per coarse range bin.
        ! The summation weights at the borders are half of the other weights,
        !  because the borders of the coarse range bins coincide with the centers
        !  of the original range bins:
        w_aggr(:)      = 1.0_dp
        w_aggr(1)      = 0.5_dp
        w_aggr(n_aggr) = 0.5_dp
        ira_offset = n_aggr / 2
      ELSE
        ! Odd number of original range bins to aggregate.
        ! The summation weights are all equal, because the borders of the coarse
        !  range bins coincide with the borders of the original range bins:
        w_aggr(:) = 1.0_dp
        ira_offset = (n_aggr-1) / 2
      END IF

      ! .. Linearize reflectivity, because it is needed as such for aggregation and weighting of other parameters:
      IF (ASSOCIATED(rs_data(ista)%zh_radar_obs)) THEN
        CALL dbz_to_linear(rs_data(ista)%zh_radar_obs(:))
      END IF
      
      ! .. Radial wind, reflectivity weigthed:
      IF (ASSOCIATED(rs_data(ista)%radwind_obs)) THEN

        IF (.NOT.ASSOCIATED(rs_data(ista)%zh_radar_obs)) THEN
          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': range aggregation of radial wind needs reflectivity obs, '// &
               'but these are missing! Range-aggregated radial wind obs are set to missing!'
          var_sum(:) = missing_obs
          CALL redefine_var (rs_data(ista)%radwind_obs, var_sum, rs_data(ista)%nobs_obs(i_vrad))
          rs_meta(ista)%lobs_avail(irep,i_vrad) = .FALSE.
          ierr = 1
        ELSE
        
          ! Aggregation to a coarser range resolution, taking into accout
          !  missing values. We assume that the obs data are
          !  not specifically sorted and that the integer 1D address rs_data(ista)%ind_intp_obs(iobs)
          !  contains the info on original range, azi and ele position.
        
          ! .. Preparation: re-sort data into helper array var_aggr(:,:) for better vectorization of the following:
          CALL resort_to_var_aggr ( rs_data(ista)%radwind_obs(:) , rs_data(ista)%ind_intp_obs_vr, var_aggr    )
          ! .. Preparation: re-sort weights into helper array var_aggr_wgt(:,:):
          CALL resort_to_var_aggr ( rs_data(ista)%zh_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr_wgt )
          ! .. Compute aggregated values:
          CALL aggr_var_sum_wgt ( var_aggr, var_aggr_wgt, var_sum )       ! Computes var_sum from var_aggr(:,:) and var_aggr_wgt(:,:), takes missing_obs into account
          ! .. Re-define storage vector in rs_data(ista):
          CALL redefine_var ( rs_data(ista)%radwind_obs, var_sum, rs_data(ista)%nobs_obs(i_vrad) )
        END IF
      END IF

      ! .. Quality flags for radial wind:
      IF (ASSOCIATED(rs_data(ista)%radwind_obs_q)) THEN
        ! .. Preparation: re-sort data into helper array for better vectorization of the following:
        CALL resort_to_var_aggr ( rs_data(ista)%radwind_obs_q(:), rs_data(ista)%ind_intp_obs_qv, var_aggr )  ! Fills var_aggr(:,:)
        ! .. Compute aggregated values:
        CALL aggr_var_max ( var_aggr(:,:), var_sum(:) )       ! Computes var_sum(:) from var_aggr(:,:)
        ! .. Re-define storage vector in rs_data(ista):
        CALL redefine_var ( rs_data(ista)%radwind_obs_q, var_sum, rs_data(ista)%nobs_obs(i_qualvrad) )
      END IF

      ! .. ZDR, splitted into ZH and ZV:
      IF (ASSOCIATED(rs_data(ista)%zdr_radar_obs)) THEN

        IF (.NOT.ASSOCIATED(rs_data(ista)%zh_radar_obs)) THEN
          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': range aggregation of ZDR needs reflectivity obs, '// &
               'but these are missing! Range-aggregated ZDR obs are set to missing!'
          var_sum(:) = missing_obs
          CALL redefine_var (rs_data(ista)%zdr_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_zdr))
          rs_meta(ista)%lobs_avail(irep,i_zdr) = .FALSE.
          ierr = 2
        ELSE

          ALLOCATE (tmpzh(nobs_max,n_aggr), tmpzv(nobs_max,n_aggr))

          ! .. 1) ZH:
          ! .. Preparation: re-sort data into helper array for better vectorization of the following:
          CALL resort_to_var_aggr (rs_data(ista)%zh_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr)  ! Fills var_aggr(:,:)
          ! .. Store ZH temporarily:
          tmpzh(:,:) = var_aggr(:,:)

          ! .. 2) ZV:
          ! .. Preparation: re-sort data into helper array var_aggr(:,:) for better vectorization of the following:
          CALL dbz_to_linear(rs_data(ista)%zdr_radar_obs(:))
          CALL resort_to_var_aggr (rs_data(ista)%zdr_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr    )
          DO i = 1, n_aggr
            DO iobs = 1, nobs_max
              IF (tmpzh(iobs,i) > miss_threshold .AND. var_aggr(iobs,i) >= Z_crit_radar) THEN
                tmpzv(iobs,i) = tmpzh(iobs,i) / var_aggr(iobs,i)
              ELSE IF (tmpzh(iobs,i) > miss_threshold) THEN
                tmpzv(iobs,i) = tmpzh(iobs,i)
              ELSE
                tmpzv(iobs,i) = missing_obs
              END IF
            END DO
          END DO

          ! .. Compute aggregated values of ZH:
          CALL aggr_var_sum ( tmpzh(:,:), var_sum(:) )
          tmpzh(:,1) = var_sum(:)
          
          ! .. Compute aggregated values of ZV:
          CALL aggr_var_sum ( tmpzv(:,:), var_sum(:) )
          ! .. Re-define storage vector in rs_data(ista):

          DO iobs = 1, nobs_max
            IF (tmpzh(iobs,1) > miss_threshold .AND. var_sum(iobs) >= Z_crit_radar) THEN
              var_sum(iobs) = tmpzh(iobs,1) / var_sum(iobs)
            ELSE IF (tmpzh(iobs,1) > miss_threshold .AND. var_sum(iobs) > miss_threshold) THEN
              var_sum(iobs) = 1.0_dp
            ELSE
              var_sum(iobs) = missing_obs
            END IF
          END DO
          CALL linear_to_dbz (var_sum)
          CALL redefine_var ( rs_data(ista)%zdr_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_zdr) )
          ! .. 3) ZDR:
          CALL linear_to_dbz(var_sum)

          DEALLOCATE (tmpzh, tmpzv)
          
        END IF
      END IF

      ! .. RHOHV, reflectivity weighted:
      IF (ASSOCIATED(rs_data(ista)%rhv_radar_obs)) THEN
        IF (.NOT.ASSOCIATED(rs_data(ista)%zh_radar_obs)) THEN
          WRITE (*,*) 'WARNING '//TRIM(yzroutine)//': range aggregation of RHOHV needs reflectivity obs, '// &
               'but these are missing! Range-aggregated RHOHV obs are set to missing!'
          var_sum(:) = missing_obs
          CALL redefine_var (rs_data(ista)%rhv_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_rhv))
          rs_meta(ista)%lobs_avail(irep,i_rhv) = .FALSE.
          ierr = 2
        ELSE
          ! .. Preparation: re-sort data into helper array var_aggr(:,:) for better vectorization of the following:
          CALL resort_to_var_aggr ( rs_data(ista)%rhv_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr   )
          ! .. Preparation: re-sort weights into helper array var_aggr_wgt(:,:):
          CALL resort_to_var_aggr ( rs_data(ista)%zh_radar_obs(:) , rs_data(ista)%ind_intp_obs_z, var_aggr_wgt )
          ! .. Compute aggregated values:
          CALL aggr_var_sum_wgt ( var_aggr, var_aggr_wgt, var_sum )       ! Computes var_sum from var_aggr(:,:) and var_aggr_wgt(:,:), takes missing_obs into account
          CALL redefine_var ( rs_data(ista)%rhv_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_rhv) )
        END IF
      END IF

      ! .. KDP, aggregation by simple bin averaging:
      IF (ASSOCIATED(rs_data(ista)%kdp_radar_obs)) THEN
        ! .. Preparation: re-sort data into helper array var_aggr(:,:) for better vectorization of the following:
        CALL resort_to_var_aggr ( rs_data(ista)%kdp_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr )
        ! .. Compute aggregated values:
        CALL aggr_var_sum ( var_aggr, var_sum )       ! Computes var_sum from var_aggr(:,:), takes missing_obs into account
        CALL redefine_var ( rs_data(ista)%kdp_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_kdp) )
      END IF

      ! .. PHIDP, aggregation by simply taking the value of the outermost original range bin within the coarse range bin:
      IF (ASSOCIATED(rs_data(ista)%phidp_radar_obs)) THEN
        ! .. Preparation: re-sort data into helper array var_aggr(:,:) for better vectorization of the following:
        CALL resort_to_var_aggr ( rs_data(ista)%phidp_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr )
        ! .. Compute aggregated values:
        CALL aggr_var_phidp ( var_aggr, var_sum )       ! Computes var_sum from var_aggr(:,:), takes missing_obs into account
        CALL redefine_var ( rs_data(ista)%phidp_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_phidp) )
      END IF

      ! .. Quality flags for ZH:
      IF (ASSOCIATED(rs_data(ista)%zh_radar_obs_q)) THEN
        ! .. Preparation: re-sort data into helper array for better vectorization of the following:
        CALL resort_to_var_aggr ( rs_data(ista)%zh_radar_obs_q(:), rs_data(ista)%ind_intp_obs_qz, var_aggr )  ! Fills var_aggr(:,:)
        ! .. Compute aggregated values:
        CALL aggr_var_max ( var_aggr(:,:), var_sum(:) )       ! Computes var_sum(:) from var_aggr(:,:)
        ! .. Re-define storage vector in rs_data(ista):
        CALL redefine_var ( rs_data(ista)%zh_radar_obs_q, var_sum, rs_data(ista)%nobs_obs(i_qualdbzh) )
      END IF

      ! .. ZH, at last: horizontal reflectivity:
      IF (ASSOCIATED(rs_data(ista)%zh_radar_obs)) THEN
        ! .. Preparation: re-sort data into helper array for better vectorization of the following:
        CALL resort_to_var_aggr ( rs_data(ista)%zh_radar_obs(:), rs_data(ista)%ind_intp_obs_z, var_aggr )  ! Fills var_aggr(:,:)
        ! .. Compute aggregated values:
        CALL aggr_var_sum ( var_aggr(:,:), var_sum(:) )       ! Computes var_sum(:) from var_aggr(:,:)
        CALL linear_to_dbz( var_sum )
        ! .. Re-define storage vector in rs_data(ista):
        CALL redefine_var ( rs_data(ista)%zh_radar_obs, var_sum, rs_data(ista)%nobs_obs(i_dbzh) )
      END IF
      
      ! .. Clean up of memory:
      DEALLOCATE (w_sum, w_aggr)
      DEALLOCATE (var_aggr, var_aggr_wgt, var_sum)

      
      ! .. As a last stept, set the new coarsened 1D-adress of each bin,
      !     which is equal to the running index of the vector (see resort_to_var_aggr() and aggr_var_XXX()):
      IF (lequal_azi_alldatasets) THEN
        ! .. This is only allocated if lequal_azi_alldatasets = .true.:
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs))       CALL redefine_intp(rs_data(ista)%ind_intp_obs)
      ELSE
        ! .. The next are allocated only if lequal_azi_alldatasets = .false.
        !     and only point to rs_data(ista)%ind_intp_obs otherwise:
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_vr))    CALL redefine_intp(rs_data(ista)%ind_intp_obs_vr)
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_qv))    CALL redefine_intp(rs_data(ista)%ind_intp_obs_qv)
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_z))     CALL redefine_intp(rs_data(ista)%ind_intp_obs_z)
        IF (ASSOCIATED (rs_data(ista)%ind_intp_obs_qz))    CALL redefine_intp(rs_data(ista)%ind_intp_obs_qz)
      END IF

    END IF
                
    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    SUBROUTINE resort_to_var_aggr (input_field_obs, input_intp_obs, varaggr)

      REAL(kind=dp), INTENT(in)  :: input_field_obs(:)
      INTEGER,       INTENT(in)  :: input_intp_obs(:)
      REAL(kind=dp), INTENT(out) :: varaggr(:,:)

      REAL(kind=dp) :: zloc
      INTEGER       :: iobs, i, m, mm, n, o, tmpind
      
      varaggr(:,:) = missing_obs
!NEC$ ivdep
      DO iobs = 1, nobs_obs(i_dbzh)

        ! .. Get original nominal (ra, az, el) indices from the 1D address, note that 1D adress is for transposed (az, ra, el) sorting!
        CALL ind2sub3D (input_intp_obs(iobs), naz, nra_obs, n, m, o)

        mm = (m + ira_offset) / n_aggr  ! The coarsened range bin to where m contributes

        IF (mm <= nra .AND. mm > 0) THEN
          i = m + ira_offset - mm*n_aggr + 1  ! The rank of m in the list of contributors to mm
          
          zloc = input_field_obs(iobs)
          ! .. Sort into helper array where the first index is the new 1D address assuming (az, ra, el) sorting:
          CALL sub2ind3D (n, mm, o, naz, nra, tmpind)
          varaggr(tmpind, i) = zloc
          IF (MOD(n_aggr,2) == 0 .AND. i == 1 .AND. mm > 1) THEN
            ! If n_aggr is even and this is the value at the lower coarse bin border, it also contributes to
            !  the previous coarse range bin:
            CALL sub2ind3D (n, mm-1, o, naz, nra, tmpind)
            varaggr(tmpind, n_aggr) = zloc
          END IF
        END IF
        
      END DO

    END SUBROUTINE resort_to_var_aggr
      
    SUBROUTINE aggr_var_sum ( varaggr, varsum )

      REAL(kind=dp), INTENT(in)  :: varaggr(:,:)
      REAL(kind=dp), INTENT(out) :: varsum(:)
      INTEGER :: iobs, i
      
      varsum(:) = 0.0_dp
      w_sum(:) = 0.0_dp
      DO i = 1, n_aggr
        DO iobs = 1, nobs_max
          IF (varaggr(iobs,i) > miss_threshold) THEN
            varsum(iobs) = varsum(iobs) + w_aggr(i)*varaggr(iobs,i)
            w_sum(iobs) = w_sum(iobs) + w_aggr(i)
          END IF
        END DO
      END DO
      DO iobs = 1, nobs_max
        IF (w_sum(iobs) > eps) THEN
          varsum(iobs) = varsum(iobs) / w_sum(iobs)
        ELSE
          varsum(iobs) = missing_obs
        END IF
      END DO
    END SUBROUTINE aggr_var_sum

    SUBROUTINE aggr_var_sum_wgt ( varaggr, varaggrwgt, varsum )

      REAL(kind=dp), INTENT(in)  :: varaggr(:,:), varaggrwgt(:,:)
      REAL(kind=dp), INTENT(out) :: varsum(:)
      INTEGER :: iobs, i
      
      varsum(:) = 0.0_dp
      w_sum(:) = 0.0_dp
      DO i = 1, n_aggr
        DO iobs = 1, nobs_max
          IF (varaggr(iobs,i) > miss_threshold .AND. varaggrwgt(iobs,i) > miss_threshold) THEN
            varsum(iobs) = varsum(iobs) + w_aggr(i)*varaggrwgt(iobs,i)*varaggr(iobs,i)
            w_sum(iobs)   = w_sum(iobs)   + w_aggr(i)*varaggrwgt(iobs,i)
          END IF
        END DO
      END DO
      DO iobs = 1, nobs_max
        IF (w_sum(iobs) > eps) THEN
          varsum(iobs) = varsum(iobs) / w_sum(iobs)
        ELSE
          varsum(iobs) = missing_obs
        END IF
      END DO
    END SUBROUTINE aggr_var_sum_wgt

    SUBROUTINE aggr_var_phidp ( varaggr, varmax )

      REAL(kind=dp), INTENT(in)  :: varaggr(:,:)
      REAL(kind=dp), INTENT(out) :: varmax(:)
      INTEGER :: iobs, i, ic

      IF (MOD(n_aggr,2) ==0) THEN
        ic = n_aggr/2
      ELSE
        ic = n_aggr/2+1
      END IF
      varmax(:) = missing_obs
      DO i = 1, ic  ! loop from innermost original range bin to the center
        DO iobs = 1, nobs_max
          ! .. Aggregated PHIDP is the value at the center of the coarse range bin, or, if this value is missing,
          !     the last value before the center which is not missing:
          IF (varaggr(iobs,i) > miss_threshold) THEN
            varmax(iobs) = varaggr(iobs,i)
          END IF
        END DO
      END DO
    END SUBROUTINE aggr_var_phidp

    SUBROUTINE aggr_var_max ( varaggr, varmax )

      REAL(kind=dp), INTENT(in)  :: varaggr(:,:)
      REAL(kind=dp), INTENT(out) :: varmax(:)
      INTEGER :: iobs, i
      
      varmax(:) = missing_obs
      DO i = 1, n_aggr  ! loop from innermost original range bin to the outermost
        DO iobs = 1, nobs_max
          ! .. Aggregated PHIDP is the outermost value in the coarse range bin which is not missing:
          IF (varaggr(iobs,i) > miss_threshold) THEN
            varmax(iobs) = MAX(varmax(iobs), varaggr(iobs,i))
          END IF
        END DO
      END DO
    END SUBROUTINE aggr_var_max

    SUBROUTINE redefine_var (datap, varsum, nobs)
      REAL(kind=dp), POINTER    :: datap(:)
      REAL(kind=dp), INTENT(in) :: varsum(:)
      INTEGER, INTENT(out)      :: nobs
      DEALLOCATE(datap); NULLIFY(datap)
      ALLOCATE(datap(nobs_max))
      datap(:) = varsum(:)
      nobs = nobs_max
    END SUBROUTINE redefine_var
    
    SUBROUTINE redefine_intp (intp)
      INTEGER, POINTER :: intp(:)
      INTEGER          :: iobs
      DEALLOCATE(intp); NULLIFY(intp)
      ALLOCATE(intp(nobs_max))
      DO iobs = 1, nobs_max
        intp(iobs) = iobs  ! The index iobs is already the 1D address assuming (az, ra, el) sorting
      END DO
    END SUBROUTINE redefine_intp
    
  END SUBROUTINE range_coarsening
  
#endif   /* NETCDF */


END MODULE radar_obs_data_read
