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

MODULE radar_obs_meta_read

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

#ifdef NETCDF
  USE radar_kind, ONLY : dp, sp
  
  USE radar_dbzcalc_params_type, ONLY : t_dbzcalc_params, dbz_namlst_d

  USE radar_data, ONLY :  &
       cvarlen,           &
       unused_value, miss_value, missval_int, miss_threshold, obsfile_missingname, &
       i_fwo_ini,         & ! Timing flag for the initialization of the forward operator
       i_fwo_barrier,     & ! Timing flag for barrier waiting in MPI-communications (measure for load imbalance)
       icomm_radar,       & ! communicator for the group of radar-IO PEs + compute PEs
       num_radar,         & ! number of radar PEs (num_compute + num_radario)
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       radar_meta_type, radar_meta_type_onetime, rsm_init_strings_blanks, &
       nradsta_max, ndatakind, nobstimes_max, cobsflen, &
       nscanstrategies_max, &
       ydate_ini_mod, &
       rsm_multitime2onetime, rsm_onetime2multitime, &
       mpi_radar_meta_typ_alltimes, mpi_radar_meta_typ_onetime, &
       nradsta_all, &
       i_vrad, i_qualvrad, i_dbzh, i_qualdbzh, &
       i_zdr, i_kdp, i_phidp, i_rhv, i_ldr, i_cflags, &
       i_dwd, i_meteoswiss, i_arpasim, i_belgium, i_denmark, i_france, i_poland, i_czech, i_netherlands, i_slovakia, &
       c_format_cdfin, &
       c_format_h5_native_smss, c_format_h5_native_smms, c_format_h5_native_mmss, c_format_h5_native_mmms,    & 
       c_format_h5_2_opera_smss, c_format_h5_2_opera_smms, c_format_h5_2_opera_mmss, c_format_h5_2_opera_mmms
       
  USE radar_data_namelist, ONLY : &
       ldebug_radsim, &
       ydirradarin, ydirlistfile, &
       loutdbz, loutradwind, lqc_flag, loutpolstd, loutpolall, &
       itype_obserr_vr, itype_supobing

  USE radar_obs_meta_list, ONLY :  &
       &        set_scanname, get_meta_proto, &
       &        get_meta_network_all, get_elarr_precipscan

  USE radar_utilities, ONLY : get_filenames_in_directory, jul2yyyymmdd, split_string, &
       &                      diff_seconds, new_datetime, round_datestr_min, tolower, toupper
  USE radar_interface, ONLY : get_runtime_timings

  USE netcdf, ONLY :  &
       nf90_open , &
       nf90_noerr, &
       NF90_NOWRITE, &
       nf90_inquire, &
       nf90_inq_dimid, &
       nf90_inquire_dimension, &
       nf90_inq_varid , &
       nf90_get_var, &
       nf90_close, &
       nf90_strerror, &
       nf90_global

#ifdef HDF5_RADAR_INPUT
  ! requires HDF5 version 1.8.8 and higher!!!
  ! requires build of hdf5 with both configure options --enable-fortran --enable-fortran2003
    USE hdf5
    USE h5lt
    USE, INTRINSIC :: ISO_C_BINDING
#endif
#endif

#ifndef NOMPI
  USE mpi
#endif
  USE radar_parallel_utilities, ONLY: distribute_values_radar, distribute_path_radar

  !================================================================================
  !================================================================================

  IMPLICIT NONE

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

  PRIVATE

  !==============================================================================
  ! Public Subroutines:

#ifdef NETCDF
  PUBLIC :: read_meta_info_all
#ifdef HDF5_RADAR_INPUT
  PUBLIC :: read_attr_odim, read_dataset_odim, get_names_of_datasets_odim_h5, &
            reconstruct_file_series_opera, reconstruct_file_series_mch
#endif
#endif

  !==============================================================================
  ! Interface blocks for overloaded procedures:

#ifdef NETCDF
  INTERFACE get_scanstrategy_ID
    MODULE PROCEDURE     &
         get_scanstrategy_ID_m, &
         get_scanstrategy_ID_o
  END INTERFACE get_scanstrategy_ID

  INTERFACE merge_metadata
    MODULE PROCEDURE     &
         merge_metadata_m, &
         merge_metadata_o
  END INTERFACE merge_metadata
#endif

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

#ifdef NETCDF

  !==============================================================================
  !+ Module procedure in radar_src for reading radar meta data from netCDF files.
  !  Is called by input_radarnamelist() if lreadmeta_from_netcdf = .TRUE.
  !  It collects rs_meta_ncdf on PE0 (if icomm_radar) and distributes
  !    nradsta to all PEs (icomm_radar).
  !------------------------------------------------------------------------------

  SUBROUTINE read_meta_info_all ( rs_meta_ncdf, nradsta_ncdf, err_ista , miss_err )

    IMPLICIT NONE

    TYPE(radar_meta_type),    INTENT(inout)  :: rs_meta_ncdf(:) ! output meta structure, merged from prototype and NetCDF files
    INTEGER,                  INTENT(out)    :: nradsta_ncdf
    INTEGER,                  INTENT(out)    :: err_ista(nradsta_max)
    INTEGER,                  INTENT(out)    :: miss_err(nradsta_max)


    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER  :: yzroutine = 'read_meta_info_all'

    TYPE(radar_meta_type)         :: rsm_dum_multitimes, rsm_proto_dwd
    TYPE(radar_meta_type_onetime) :: rsm_dum_onetime
    TYPE(radar_meta_type), ALLOCATABLE :: rs_meta_buf_multitimes(:)
    TYPE(radar_meta_type_onetime), ALLOCATABLE :: rs_meta_buf_onetime(:)
    CHARACTER (LEN=80)            :: yzerrmsg
    INTEGER                       :: i, j, ii, ii_par, ii_end, ista, ista_multi, ista_one, n, &
                                     ipe_rad, ierr, ninfiles_good, ninfiles_onetime, ninfiles_multitimes, ninfiles, elecnt, &
                                     station_id, mpierr, mpistat(MPI_STATUS_SIZE), ista_end, nwords, fid, &
                                     i3u, i3o, nradsta_multi, nradsta_one

    CHARACTER(len=cobsflen),ALLOCATABLE :: infilenames(:)
    CHARACTER(len=15),      ALLOCATABLE :: infileformat(:)
    CHARACTER(len=15),      ALLOCATABLE :: infilevarname(:)
    INTEGER,                ALLOCATABLE :: infilecountry(:)
    LOGICAL,                ALLOCATABLE :: has_multiple_obstimes(:)
    CHARACTER(len=cobsflen),ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    CHARACTER(len=15)                   :: format_from_file, c_format_loc
    INTEGER                             :: icountry_loc

    TYPE(radar_meta_type)      :: rs_meta_network(nradsta_all)
    TYPE(t_dbzcalc_params)     :: dbz_meta_network(nradsta_all)
    LOGICAL                    :: fileexist(nradsta_max), checkflags(1:ndatakind), good, &
                                  have_read_a_file(0:num_radar-1), missing_obs_time, &
                                  dirlistfile_exist, first_appear, is_multisweep

    CHARACTER(len=300)         :: errstring
    CHARACTER(len=100)         :: checkflagtext(1:ndatakind)
    CHARACTER(len=30)          :: tmpword
    CHARACTER(len=14)          :: datestr2
    CHARACTER(len=14), ALLOCATABLE :: datestr1vec(:)
    CHARACTER(len=1)           :: elechar

    CHARACTER(len=cvarlen) :: obs_dwd_hdf5_varname_dbzh
    CHARACTER(len=cvarlen) :: obs_dwd_hdf5_varname_vrad
    CHARACTER(len=cvarlen) :: obs_dwd_hdf5_varname_zdr
    CHARACTER(len=cvarlen) :: obs_dwd_hdf5_varname_kdp
    CHARACTER(len=cvarlen) :: obs_dwd_hdf5_varname_phidp
    CHARACTER(len=cvarlen) :: obs_dwd_hdf5_varname_rhv


    !==============================================================================
    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE read_meta_info_all
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    err_ista(:) = 0
    miss_err(:) = 0

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_ini)
#endif

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    CALL rsm_init_strings_blanks(rs_meta_ncdf)
    
    ! .. To check for DWD native hdf5-files, we need its background prototype structure
    !    and for efficiency reasons the standard moment names in lowercase:
    rsm_proto_dwd = get_meta_proto ( icountry=i_dwd )
    obs_dwd_hdf5_varname_dbzh  = tolower(rsm_proto_dwd%obs_hdf5_varname_dbzh )
    obs_dwd_hdf5_varname_vrad  = tolower(rsm_proto_dwd%obs_hdf5_varname_vrad )
    obs_dwd_hdf5_varname_zdr   = tolower(rsm_proto_dwd%obs_hdf5_varname_zdr  )
    obs_dwd_hdf5_varname_kdp   = tolower(rsm_proto_dwd%obs_hdf5_varname_kdp  )
    obs_dwd_hdf5_varname_phidp = tolower(rsm_proto_dwd%obs_hdf5_varname_phidp)
    obs_dwd_hdf5_varname_rhv   = tolower(rsm_proto_dwd%obs_hdf5_varname_rhv  )
      
    IF (my_radar_id == 0) THEN

      IF (LEN_TRIM(ydirlistfile) > 0) THEN
        INQUIRE(file=TRIM(ydirlistfile), exist=dirlistfile_exist)
        IF ( .NOT.dirlistfile_exist ) THEN
          WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//'(): directory listing file '// &
                         TRIM(ydirlistfile)//' does not exist. Creating temporary file in rundir!'
        END IF
      ELSE
        dirlistfile_exist = .FALSE.
      END IF
      IF (dirlistfile_exist) THEN
        CALL get_filenames_in_directory (TRIM(ydirradarin), LEN(infilenames), infilenames, ninfiles, ierr, &
                                         my_radar_id, TRIM(ydirlistfile))
      ELSE
        CALL get_filenames_in_directory (TRIM(ydirradarin), LEN(infilenames), infilenames, ninfiles, ierr, &
                                         my_radar_id)
      END IF

      ALLOCATE (has_multiple_obstimes(ninfiles))
      has_multiple_obstimes = .FALSE.
      ALLOCATE (infileformat(ninfiles))
      infileformat(:)(:) = ' '
      ALLOCATE (infilecountry(ninfiles))
      infilecountry = -1
      ALLOCATE (infilevarname(ninfiles))
      infilevarname(:)(:) = ' '
      ALLOCATE (datestr1vec(ninfiles))
      datestr1vec(:)(:) = ' '

      IF (ninfiles >= 1000) THEN
        WRITE (*,'(a,i7,a,i7)') '  emvorado::'//TRIM(yzroutine)//': Checking radar obs filename no. ', 1, ' from ', ninfiles
      ENDIF

      ! .. Pre-filter the known files, so that ninfiles will be not too large
      !     when allocating rs_meta_buf(:) below:
      ninfiles_good       = 0
      ninfiles_multitimes = 0
      ninfiles_onetime    = 0
 
      filter_loop: DO ii=1, ninfiles

        IF (ninfiles >= 1000 .AND. MOD(ii,1000) == 0) THEN
          WRITE (*,'(a,i7,a,i7)') '  emvorado::'//TRIM(yzroutine)//': Checking radar obs filename no. ', &
                                  ii, ' from ', ninfiles
        ENDIF

        ! 1) Check for DWD hdf5 files:

#ifdef HDF5_RADAR_INPUT
        CALL split_string (TRIM(infilenames(ii)), '_', LEN(words), words, nwords)
        IF (nwords == 4 .AND. LEN_TRIM(words(1)) >= 13) THEN
          IF (words(1)(1:3) == 'ras' .AND. INDEX(TRIM(words(4)), 'hd5') > 0) THEN

            ! 1a) Check for DWD single-sweep multi-moment hdf5 files of filename-convention
            !     "rasXX-pcpng01_sweeph5allm_any_00-2020020300053400-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_00-2020020300055800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_01-2020020300062100-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_02-2020020300064500-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_03-2020020300070800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_04-2020020300073100-neu-10557-hd5"
            !
            ! and, in one go,
            !
            ! 1b) Check for DWD single-sweep single-moment hdf5 files of filename-convention
            !     "rasXX-pcpng01_sweeph5onem_dbzh_00-2020020300053400-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_00-2020020300055800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_01-2020020300062100-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_02-2020020300064400-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_03-2020020300070800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_04-2020020300073100-neu-10557-hd5"

            IF ( ( (TRIM(words(2)) == 'sweeph5allm' .AND. TRIM(words(3)) == 'any') .OR. &
                 TRIM(words(2)) == 'sweeph5onem' ) .AND. &
                 (words(1)(7:13) == 'vol5min' .OR. words(1)(7:9) == 'pcp') ) THEN

              IF (TRIM(words(2)) == 'sweeph5onem') THEN
                IF (        words(3)(1:4) == obs_dwd_hdf5_varname_vrad(1:4)       .OR. &
                     INDEX(TRIM(words(3)), TRIM(obs_dwd_hdf5_varname_zdr))   /= 0 .OR. &
                     INDEX(TRIM(words(3)), TRIM(obs_dwd_hdf5_varname_kdp))   /= 0 .OR. &
                     INDEX(TRIM(words(3)), TRIM(obs_dwd_hdf5_varname_phidp)) /= 0 .OR. &
                     INDEX(TRIM(words(3)), TRIM(obs_dwd_hdf5_varname_rhv))   /= 0 .OR. &
                           TRIM(words(3)) == TRIM(obs_dwd_hdf5_varname_dbzh)           &
                     ) THEN
                  ! Should be one of 'vradh','dbzh','zdr','kdp','uphidp','urhohv'
                  CONTINUE
                ELSE
                  IF (ldebug_radsim) THEN
                    WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // '(): datakind=', TRIM(words(3)), &
                         ' not implemented for DWD hdf5 single-moment files! File discarded! '//TRIM(infilenames(ii))
                  END IF
                  CYCLE filter_loop
                END IF
              END IF

              ! Only keep the filename if the time stamp in word(4) is unique after being rounded
              !  down to the next full 5 minute interval.
              ! Otherwise, another filename of the same volume scan already has been stored in the
              !  list of infilenames, and we only want to store unique volume scan indicators and
              !  variable/qantity/moment names (words(1)+words(4) and words(3)):
              ! For this, we assume that all elevations of a volume scan start within the same 5min
              !  interval, that is, between minutes 0 and 5, 5 and 10, 10 and 15 after the hour.
              !  Is that always the case???
              first_appear = .TRUE.
              datestr2 = round_datestr_min(words(4)(4:17), -5)  ! round time stamp of actual file down to the next full 5 min for comparison with previous files
              IF ( words(1)(7:13) == 'vol5min' ) THEN   ! actual file is a volscan
                DO i=1, ninfiles_good
                  IF ( infilecountry(i) == i_dwd .AND. &
                       ( TRIM(infileformat(i)) == TRIM(c_format_h5_native_smss) .OR. &
                         TRIM(infileformat(i)) == TRIM(c_format_h5_native_mmss) ) ) THEN
                    IF ( infilenames(i)(7:13) == 'vol5min' ) THEN   ! infilenames(i) is a volscan
                      i3o = INDEX(infilenames(i),'_',back=.TRUE.)   ! position of the last '_' in the string
                      IF ( datestr1vec(i) == datestr2                               .AND. &   ! same 5 min time interval
                           infilenames(i)(i3o+21:i3o+29) == words(4)(21:29)   .AND. &   ! same radar station
                           INDEX(infilenames(i),'_'//TRIM(words(3))//'_') > 0 .AND. &   ! same quantity/moment
                           INDEX(infilenames(i),'_'//TRIM(words(2))//'_') > 0  &        ! same file type (sweeph5onem, sweeph5allm)
                         ) THEN
                        first_appear = .FALSE.
                        ! Increment the elevation counter in the stored filename by 1:
                        READ  (infilenames(i)(i3o+1:i3o+2), '(i2.2)') elecnt
                        WRITE (infilenames(i)(i3o+1:i3o+2), '(i2.2)') elecnt + 1
                        ! Add eleindex and mmss part of the time stamp to the end of the filename:
                        infilenames(i) = TRIM(infilenames(i)) // '|' // words(4)(1:2) // '-' // words(4)(14:17)
                        EXIT
                      END IF
                    END IF
                  END IF
                END DO
              END IF
              IF (first_appear) THEN
                ninfiles_good                        = ninfiles_good + 1
                infilenames(ninfiles_good)           = infilenames(ii)
                i3o                                  = INDEX(infilenames(ii),'_',back=.TRUE.)   ! position of the last '_' in the string
                datestr1vec(ninfiles_good)           = round_datestr_min(infilenames(ii)(i3o+4:i3o+17), -5)  ! round time stamp of infilesnames(ii) down to the next full 5 min and store for later comparison
                infilecountry(ninfiles_good)         = i_dwd
                IF (TRIM(words(2)) == 'sweeph5onem') THEN
                  infileformat(ninfiles_good)        = c_format_h5_native_smss  ! smss: single-moment single-sweep
                  infilevarname(ninfiles_good)       = tolower(TRIM(words(3)))
                ELSE
                  infileformat(ninfiles_good)        = c_format_h5_native_mmss  ! mmss: multi-moment single-sweep
                  infilevarname(ninfiles_good)       = 'all'
                END IF
                IF (words(1)(7:13) == 'vol5min') THEN
                  ! Add eleindex and mmss part of the time stamp to the end of the filename:
                  infilenames(ninfiles_good) = TRIM(infilenames(ninfiles_good)) // '|' // words(4)(1:2) // '-' // words(4)(14:17)
                  ! Initialize the elevation counter in the stored filename with 1:
                  elecnt = 1
                  WRITE (infilenames(ninfiles_good)(i3o+1:i3o+2), '(i2.2)') elecnt
                END IF
                has_multiple_obstimes(ninfiles_good) = .FALSE.
                ninfiles_onetime                     = ninfiles_onetime + 1
              END IF
              CYCLE filter_loop

            END IF

          END IF
        END IF
#endif

        ! 2) Check for DWD CDFIN-files of filename-convention
        !     "cdfin_z_id-010908_201307282200_201307282255_volscan" or
        !     "cdfin_z_id-010908_201307282200_201307282255_precipscan"

        CALL split_string (TRIM(infilenames(ii)), '_', LEN(words), words, nwords)
        IF (INDEX(TRIM(infilenames(ii)), 'cdfin') > 0 .AND. nwords == 6) THEN
          IF (TRIM(words(6)) == 'volscan' .OR. TRIM(words(6)) == 'precipscan') THEN

            SELECT CASE (TRIM(words(2)))
            CASE ('vr','z','vrsim','zrsim')    ! 'vrsim' and 'zrsim' are for idealized OSSE-runs
              CONTINUE
            CASE ('qv','qz','qvobs','qzobs')
              IF ( .NOT. lqc_flag ) THEN
                ! In case of lqc_flag = .false. the files for the quality flags qv and qz are not needed,
                !  so just skip them:
                CYCLE filter_loop
              END IF
            CASE default
              WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // '(): datakind=', TRIM(words(2)), &
                   ' not implemented! File discarded! '//TRIM(infilenames(ii))
              CYCLE filter_loop
            END SELECT
            ninfiles_good                        = ninfiles_good + 1
            infilenames(ninfiles_good)           = infilenames(ii)
            infilecountry(ninfiles_good)         = i_dwd
            infileformat(ninfiles_good)          = c_format_cdfin  ! single-moment multi-sweep
            infilevarname(ninfiles_good)         = TRIM(words(2))
            has_multiple_obstimes(ninfiles_good) = .TRUE.
            ninfiles_multitimes                  = ninfiles_multitimes + 1
            CYCLE filter_loop

          END IF
        END IF

        ! 3) Check for Swiss files of filename-convention
        !     "PLA1520500007U.001.V.h5"
        !     "PLA1520500007U.001.Z.h5"
        !     "PLA1520500007U.002.V.h5"
        !     "PLA1520500007U.002.Z.h5"
        !    but eliminate files with doubled first word, because this is the same for all elevations.
        !    Later, read all files from one time and one station (=one volscan) on one processer,
        !     but save only the first word in the file list in rs_meta.
        !
        ! Email from Lorenzo Clementi (MCH): "ML" indicates  that it is a polar sweep file (followed by the letter
        !  indicating which radar site it belongs to). In the past, until autumn 2017, the first letter was "P",
        !  it has changed to "M" when we moved from Solaris to Linux and the internal DATA FORMAT has changed
        !  accordingly. "ML" is therefore fixed now, it will not change in the foreseeable future.
        ! There are various possibilities for the "0U" or "7U" part (0 could be either 0 or 1, U could be another letter).
        !  If you are writing a program, I would suggest that you use a wildcard to match those two letters.


#ifdef HDF5_RADAR_INPUT
        CALL split_string (TRIM(infilenames(ii)), '.', LEN(words), words, nwords)
        ! if there are 4 words and the file starts with an M or a P, followed by an L, it is a candidate:
        IF (nwords == 4 .AND. SCAN(TRIM(infilenames(ii)), 'MP') == 1 .AND. SCAN(TRIM(infilenames(ii)), 'L') == 2) THEN
          IF (TRIM(words(4)) == 'h5' .AND. LEN_TRIM(words(1)) == 14 .AND. LEN_TRIM(words(2)) == 3) THEN
            SELECT CASE (TRIM(words(3)))
            CASE ('V','Z','ZV','VZ')
              CONTINUE
            CASE default
              WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // '(): datakind=', TRIM(words(3)), &
                   ' not implemented! File discarded! '//TRIM(infilenames(ii))
              CYCLE filter_loop
            END SELECT

            ! Only keep the filename if the first word and third words are unique. Otherwise, another filename
            !  of the same volume scan already has been stored in the list of infilenames,
            !  and we only want to store the volume scan indicator and variable name (words(1) and words(3)):
            first_appear = .TRUE.
            ! work-around for a possible compiler bug on Cray:
            !  the expression "infilenames(i)(20:20+LEN_TRIM(words(3))-1)" fails, so we
            !  do it with "infilenames(i)(i3u:i3o)":
            i3u = 20
            i3o = i3u + LEN_TRIM(words(3)) - 1
            DO i=1, ninfiles_good
              IF ( infilecountry(i)== i_meteoswiss .AND. &
                   ( TRIM(infileformat(i)) == TRIM(c_format_h5_native_smss) .OR. &
                     TRIM(infileformat(i)) == TRIM(c_format_h5_native_mmss) ) ) THEN
                tmpword(:) = ' '
                !            tmpword = infilenames(i)(20:20+LEN_TRIM(words(3))-1)
                tmpword = infilenames(i)(i3u:i3o)
                IF (infilenames(i)(1:14) == TRIM(words(1)) .AND. &
                     TRIM(tmpword) == TRIM(words(3))) THEN
                  first_appear = .FALSE.
                  ! Increment the elevation counter in the stored filename by 1:
                  READ (infilenames(i)(16:18), '(i3.3)') elecnt
                  WRITE (infilenames(i)(16:18), '(i3.3)') elecnt + 1
                  ! Append the actual elevation identifier to the end of the filename:
                  infilenames(i) = TRIM(infilenames(i)) // '|' // words(2)(1:3)
                  EXIT
                END IF
              END IF
            END DO
            IF (first_appear) THEN
              ninfiles_good                        = ninfiles_good + 1
              infilenames(ninfiles_good)           = infilenames(ii)
              ! Initialize the elevation counter in the stored filename with 1: (just to be sure...)
              elecnt = 1
              WRITE (infilenames(ninfiles_good)(16:18), '(i3.3)') elecnt
              infilenames(ninfiles_good) = TRIM(infilenames(ninfiles_good)) // '|' // words(2)(1:3)
              infilecountry(ninfiles_good)         = i_meteoswiss
              infilevarname(ninfiles_good)         = TRIM(words(3))
              IF (LEN_TRIM(infilevarname(ninfiles_good)) > 1) THEN
                infileformat(ninfiles_good)        = c_format_h5_native_mmss  ! mmss: multi-moment single-sweep
              ELSE
                infileformat(ninfiles_good)        = c_format_h5_native_smss  ! smss: single-moment single-sweep
              END IF
              has_multiple_obstimes(ninfiles_good) = .FALSE.
              ninfiles_onetime                     = ninfiles_onetime + 1
            END IF
            CYCLE filter_loop

          END IF
        END IF
#endif


        ! 4) Check for OPERA files which are distributed to OPERA from NWS's and which have filename-convention
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

#ifdef HDF5_RADAR_INPUT
        CALL split_string (TRIM(infilenames(ii)), '_', LEN(words), words, nwords)
        ! if there are 5 words and the file starts with a T_PA then it is a candidate:
        IF ( nwords == 5 .AND. infilenames(ii)(1:4) == 'T_PA') THEN
          IF ( 17 <= LEN_TRIM(words(5)) .AND. LEN_TRIM(words(5)) <= 18 ) THEN

            ! This is an OPERA-file collected from the original countries at DWD:

!!$ From Belgium: EBUM  use T_PAGX... (DBZH of 9 elevations for corrected DBZH) and T_PAHZ... (VRAD of 9 elevations for VRAD) --> single-moment multi-sweep
!!$                2 stations: 41 (WID Wideumont) und 42 (JAB Jabekke)
!!$ From France: LFPW  --> multi-moment single-sweep, use the C-Bands only, by NOT putting the S-Bands into the background list.
!!$              For French radars, only the lowest 3 elevations are always present. There are not more than 6 actual Elevations in total in a volume scan.
!!$ From Denmark: EKMI use T_PAZZ...  (multi-moment multi-sweep)  Alternating scans optimized
!!$                for DBZH to 240 km and VRAD to 120 km, variable rsm%nra! about 240 (VRAD scan) or about 480 (DBZH scan).
!!$               To handle these data, we fix rsm%nra at 480
!!$ From Poland: SOWR make 10' volume scans, but within the 10' there are two types of scans on different elevation and resolution sets,
!!$               one for DBZH (GZ) starting at minute 0 and another for VRADH (HZ) starting at minute 3.
!!$              We take the ones for DBZH ONLY, because for the VRADH DATA we have no corresponding reflectivities for quality control.
!!$ From Czech Republic: OKPR, time in filename is endtime of scan! 5' scans, DBZH and VRADH
!!$
!!$ From Netherlands: 2 Radars, but only 3 low elevations where we have constant ra_inc. DBZH and VRADH
!!$
!!$ From Slovakia: LZIB, time in filename is nominal start time of scan. 5' scans, DBZH and VRADH
!!$

            IF ( ( TRIM(words(4)) == 'LSSW' .AND. words(5)(16:18) == 'hdf'  .AND. LEN_TRIM(infilenames(ii)) == 34 ) .OR. &
                 ( TRIM(words(4)) == 'LFPW' .AND. words(5)(16:18) == 'h5'   .AND. LEN_TRIM(infilenames(ii)) == 33  .AND. &
                 &                                             (words(2)(3:3) == 'Z') ) .OR. &
                 ( TRIM(words(4)) == 'EDZW' .AND. words(5)(16:18) == 'hdf'  .AND. LEN_TRIM(infilenames(ii)) == 34 .AND.  &
                 &                                             (words(2)(3:3) == 'G' .OR. words(2)(3:3) == 'H') ) .OR. &
                 ( TRIM(words(4)) == 'EBUM' .AND. words(5)(16:18) == 'hdf'  .AND. LEN_TRIM(infilenames(ii)) == 34 .AND.  &
                 &                                             (words(2)(3:4) == 'GX' .OR. words(2)(3:4) == 'HZ') ) .OR. &
                 ( TRIM(words(4)) == 'SOWR' .AND. words(5)(16:17) == 'h5'  .AND. LEN_TRIM(infilenames(ii)) == 33 .AND.  &
                 &                                             (words(2)(3:4) == 'GZ') ) .OR. &
                 ( TRIM(words(4)) == 'OKPR' .AND. words(5)(16:18) == 'hdf'  .AND. LEN_TRIM(infilenames(ii)) == 34 .AND.  &
                 &                                             (words(2)(3:4) == 'GZ' .OR. words(2)(3:4) == 'HZ') ) .OR. &
                 ( TRIM(words(4)) == 'LZIB' .AND. words(5)(16:18) == 'hdf'  .AND. LEN_TRIM(infilenames(ii)) == 34 .AND.  &
                 &                                             (words(2)(3:4) == 'GZ' .OR. words(2)(3:4) == 'HZ') ) .OR. &
                 ( TRIM(words(4)) == 'KITC' .AND. words(5)(16:18) == 'h5'  .AND. LEN_TRIM(infilenames(ii)) == 33 .AND.  &
                 &                                             (words(2)(3:4) == 'GZ' .OR. words(2)(3:4) == 'HZ') ) .OR. &
                 ( TRIM(words(4)) == 'EHDB' .AND. words(5)(16:17) == 'h5'  .AND. LEN_TRIM(infilenames(ii)) == 33 .AND.  &
                 &                                             (words(2)(3:4) == 'GZ') ) .OR. &
                 ( TRIM(words(4)) == 'EKMI' .AND. words(5)(16:17) == 'h5'   .AND. LEN_TRIM(infilenames(ii)) == 33 .AND.  &
                 &                                             (words(2)(3:4) == 'ZZ') ) &
                 ) THEN

              ! Discriminate ODIM single-/multi-moment single-/multi-sweep files:
              SELECT CASE (words(4))
              CASE ('LSSW')
                c_format_loc = c_format_h5_2_opera_mmss  ! mmss: multi-moment single-sweep
                icountry_loc = i_meteoswiss
                is_multisweep = .FALSE.
              CASE ('LFPW')
                c_format_loc = c_format_h5_2_opera_mmss  ! mmss: multi-moment single-sweep
                icountry_loc = i_france
                is_multisweep = .FALSE.
              CASE ('EDZW')
                c_format_loc = c_format_h5_2_opera_smss  ! smms: single-moment multi-sweep
                icountry_loc = i_dwd
                is_multisweep = .FALSE.
              CASE ('EBUM')
                c_format_loc = c_format_h5_2_opera_smms  ! smms: single-moment multi-sweep
                icountry_loc = i_belgium
                is_multisweep = .TRUE.
              CASE ('EKMI')
                c_format_loc = c_format_h5_2_opera_mmms  ! smms: multi-moment multi-sweep
                icountry_loc = i_denmark
                is_multisweep = .TRUE.
              CASE ('SOWR')
                c_format_loc = c_format_h5_2_opera_smms  ! smms: single-moment multi-sweep
                icountry_loc = i_poland
                is_multisweep = .TRUE.
              CASE ('OKPR')
                c_format_loc = c_format_h5_2_opera_smms  ! smms: single-moment multi-sweep
                icountry_loc = i_czech
                is_multisweep = .TRUE.
              CASE ('LZIB')
                c_format_loc = c_format_h5_2_opera_smms  ! smms: single-moment multi-sweep
                icountry_loc = i_slovakia
                is_multisweep = .TRUE.
              CASE ('KITC')   ! KIT-Cube
                c_format_loc = c_format_h5_2_opera_smms  ! smms: single-moment multi-sweep
                icountry_loc = i_dwd  ! we abuse i_dwd, because the file format is similar to DWD smms files
                is_multisweep = .TRUE.
              CASE ('EHDB')
                c_format_loc = c_format_h5_2_opera_mmms  ! mmms: single-moment multi-sweep
                icountry_loc = i_netherlands
                is_multisweep = .TRUE.
              END SELECT
              
              ! Only keep the filename if the station number, scantype identifier, distributing center and the
              !  datetimestring (rounded down to the next full 5 minutes) is unique. Otherwise, another filename
              !  of the same volume scan already has been stored in the list of infilenames,
              !  and we only want to store the filename for the first elevation and the number of elevations:
              first_appear = .TRUE.
              datestr2 = round_datestr_min(words(5)(1:14), -5)   ! round time stamp of actual file down to the next full 5 min for comparison with previous files
              DO i=1, ninfiles_good
                IF (infilecountry(i) == icountry_loc .AND. TRIM(infileformat(i)) == TRIM(c_format_loc)) THEN
                  IF ( infilenames(i)(12:15) == TRIM(words(4))      .AND. & ! check distributing center
                       datestr1vec(i)        == datestr2            .AND. & ! datetimestring rounded down to full 5'
                       infilenames(i)(7:8)   == TRIM(words(2)(5:6)) .AND. & ! station number
                       infilenames(i)(5:5)   == TRIM(words(2)(3:3)) ) THEN  ! scantype identifier
                    first_appear = .FALSE.
                    ! Increment the elevation counter in the stored filename by 1:
                    READ (infilenames(i)(6:6), '(a1)') elechar
                    WRITE (infilenames(i)(6:6), '(a1)') ACHAR( IACHAR(elechar) + 1 )
                    ! Append the actual elevation identifier and exact minutes/seconds to the end of the filename:
                    infilenames(i) = TRIM(infilenames(i)) // '|' // words(2)(4:4) // '-' // words(5)(11:14)
                    EXIT
                  END IF
                END IF
              END DO
              IF (first_appear) THEN
                ninfiles_good                        = ninfiles_good + 1
                infilenames(ninfiles_good)           = infilenames(ii)
                datestr1vec(ninfiles_good)           = round_datestr_min(infilenames(ii)(17:30), -5) ! round time stamp of infilesnames(ii) down to the next full 5 min and store for later comparison
                IF (.NOT.is_multisweep) THEN
                  ! Initialize the elevation counter in the stored filename with A:
                  elechar = 'A'
                  WRITE (infilenames(ninfiles_good)(6:6), '(a1)') elechar
                  ! Append informations on elevation identifier and exact minutes/seconds to the filename:
                  infilenames(ninfiles_good)         = TRIM(infilenames(ninfiles_good)) // &
                                                         '|' // words(2)(4:4) // '-' // words(5)(11:14)
                END IF
                infilecountry(ninfiles_good)         = icountry_loc
                infileformat(ninfiles_good)          = c_format_loc
                infilevarname(ninfiles_good)         = 'all'
                has_multiple_obstimes(ninfiles_good) = .FALSE.
                ninfiles_onetime                     = ninfiles_onetime + 1
              END IF
              CYCLE filter_loop
              
            END IF
            
          END IF
        END IF
#endif

        
        ! 5) Check for Italian files of filename-convention
        !     "odim_201410090010_01"
        !     "odim_201410090010_02"    (02 = station index)

#ifdef HDF5_RADAR_INPUT
        CALL split_string (TRIM(infilenames(ii)), '_', LEN(words), words, nwords)
        IF (nwords == 3 .AND. INDEX(TRIM(infilenames(ii)), 'odim') > 0) THEN
        ! THOMAS & VIRGI modification: filename
          IF (TRIM(words(1)) == 'odim' .AND. LEN_TRIM(words(2)) == 12 .AND. LEN_TRIM(words(3)) == 5) THEN

            ninfiles_good                        = ninfiles_good + 1
            infilenames(ninfiles_good)           = infilenames(ii)
            infilecountry(ninfiles_good)         = i_arpasim
            infileformat(ninfiles_good)          = c_format_h5_native_mmms  ! mmms: multi-moment multi-sweep
            infilevarname(ninfiles_good)         = 'all'
            has_multiple_obstimes(ninfiles_good) = .FALSE.
            ninfiles_onetime                     = ninfiles_onetime + 1
            CYCLE filter_loop

          END IF
        END IF
#endif

        ! 6) Check for KIT C-band files of filename-convention
        !     "scan-sidpol-120km-14_20001_20230712123503_00.h5"
        !     "scan-sidpol-120km-14_20001_20230712124001_00.h5"
        !
        ! The station is listed under i_dwd and produces mmms files, which the original DWD-radars don't

#ifdef HDF5_RADAR_INPUT
        CALL split_string (TRIM(infilenames(ii)), '_', LEN(words), words, nwords)
        IF (nwords == 4 .AND. INDEX(TRIM(infilenames(ii)), 'scan-sidpol') > 0) THEN
          IF (words(1)(1:11) == 'scan-sidpol' .AND. LEN_TRIM(words(2)) == 5 .AND. &
               LEN_TRIM(words(3)) == 14 .AND. words(4)(4:5) == 'h5') THEN

            ninfiles_good                        = ninfiles_good + 1
            infilenames(ninfiles_good)           = infilenames(ii)
            infilecountry(ninfiles_good)         = i_dwd  ! we abuse i_dwd here, because DWD itself does not provide mmms files
            infileformat(ninfiles_good)          = c_format_h5_native_mmms  ! mmms: multi-moment multi-sweep
            infilevarname(ninfiles_good)         = 'all'
            has_multiple_obstimes(ninfiles_good) = .FALSE.
            ninfiles_onetime                     = ninfiles_onetime + 1
            CYCLE filter_loop

          END IF
        END IF
#endif

        ! 7) Possible follow-on checks for other countries:
        !
        !     ...

      END DO filter_loop

      ! Clean up the helper vector datestr1vec:
      DEALLOCATE(datestr1vec)

    END IF   ! my_radar_id == 0


    ! .. Distribute valid infilenames and related informations to all PEs:

    IF (num_radar > 1) THEN
      CALL distribute_values_radar (ninfiles,            1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (ninfiles_good,       1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (ninfiles_onetime,    1, 0, icomm_radar, ierr)
      CALL distribute_values_radar (ninfiles_multitimes, 1, 0, icomm_radar, ierr)
    END IF

    IF (ninfiles < 1) THEN
      IF (my_radar_id == 0) THEN
        WRITE(*,'(2a)') TRIM(yzroutine) // ' WARNING: No radar input files found in directory ''' &
             //TRIM(ydirradarin)//''''
      END IF
      miss_err = 1
      RETURN
    END IF

    IF (ninfiles_good < 1) THEN
      IF (my_radar_id == 0) THEN
        WRITE(*,'(2a)') TRIM(yzroutine) // ' WARNING: No valid radar files found in directory ''' &
             //TRIM(ydirradarin)//''''
      END IF
      miss_err = 2
      RETURN
    END IF

    IF (num_radar > 1) THEN
      IF (my_radar_id /= 0) THEN
        ALLOCATE (infilenames(ninfiles_good))
        infilenames(:)(:) = ' '
        ALLOCATE (infileformat(ninfiles_good))
        infileformat(:)(:) = ' '
        ALLOCATE (infilecountry(ninfiles_good))
        infilecountry(:) = -1
        ALLOCATE (infilevarname(ninfiles_good))
        infilevarname(:)(:) = ' '
        ALLOCATE (has_multiple_obstimes(ninfiles_good))
        has_multiple_obstimes = .FALSE.
      END IF
      CALL distribute_values_radar (has_multiple_obstimes(1:ninfiles_good), ninfiles_good, 0, icomm_radar, ierr)
      CALL distribute_values_radar (infilecountry        (1:ninfiles_good), ninfiles_good, 0, icomm_radar, ierr)
      DO i=1, ninfiles_good
        CALL distribute_path_radar (infilenames(i),  icomm_radar)
        CALL distribute_path_radar (infileformat(i), icomm_radar)
        CALL distribute_path_radar (infilevarname(i),icomm_radar)
      END DO
    END IF


    ! .. Get the metadata background list of all known radar stations, to check / complement
    !    the actual meta data in the radar files against them:

    CALL get_meta_network_all ( dbz_namlst_d, rs_meta_network, dbz_meta_network )

    ! dbz_meta_network not used further!

    ! .. From all valid files in ydirradarin, determine whether they are radar files, and if yes, determine
    !    their scan strategy and their country-flag. Do this in parallel on all PEs to save
    !    one global communication for the meta data.
    !    A station is defined by having a certain station_id and a certain scan strategy.
    !    We expect that one input file for one station contains only scans with the exact same strategy.
    !     Different scan strategies for the same station require different input files and will lead
    !     to different "radars" internally.
    !    If there are more files from one station with the same scan strategy (e.g., files with hourly batches),
    !     add their obs_times to the existing rs_meta(ista) for this station.
    fileexist(:) = .FALSE.
    ista_multi = 0
    ista_one   = 0

    IF (my_radar_id == 0) THEN
      ALLOCATE (rs_meta_buf_multitimes(ninfiles_multitimes))
      ALLOCATE (rs_meta_buf_onetime   (ninfiles_onetime   ))
    ELSE
      ! dummy allocation for MPI-calls below:
      ALLOCATE (rs_meta_buf_multitimes(1))
      ALLOCATE (rs_meta_buf_onetime(1))
    END IF
    CALL rsm_init_strings_blanks (rs_meta_buf_multitimes)
    CALL rsm_init_strings_blanks (rs_meta_buf_onetime)

    IF (my_radar_id == 0 .AND. ninfiles_good >= 500) THEN
      WRITE (*,'(a,i7,a,i7)') '  emvorado::'//TRIM(yzroutine)//': Processing radar obs volume set no. ', &
                              1, ' from ', ninfiles_good
    ENDIF

    collect_loop: DO ii_par = 1, ninfiles_good, num_radar

      ii_end = MIN(ii_par+num_radar-1, ninfiles_good)

      ! Scan the files, filter for known filenames (radar files) and
      !  determine the station_id, scan strategy and from these,
      !  the station-identifier: Also rsm%icountry
      !------------------------------------------------------------

      have_read_a_file(:) = .FALSE.
      read_loop_0: DO ii = ii_par, ii_end

        IF (my_radar_id == 0 .AND. ninfiles_good >= 500 .AND. MOD(ii,500) == 0) THEN
          WRITE (*,'(a,i7,a,i7)') '  emvorado::'//TRIM(yzroutine)//': Processing radar obs volume set no. ', &
                                  ii, ' from ', ninfiles_good
        ENDIF

        ! determine the number of the PE for output of the radar with no. ista:
        ipe_rad = MOD(ii-1, num_radar)

        IF (my_radar_id == ipe_rad) THEN

            ! 0) Initialization of the dummy rsm_dum of type(radar_meta_type) to hold
            !     the configuration of the radar station from the actual input file:
            rsm_dum_onetime%nobs_times_obs    = 0
            rsm_dum_multitimes%nobs_times_obs = 0

#ifdef HDF5_RADAR_INPUT
          IF (infilecountry(ii) == i_dwd .AND. &
               ( TRIM(infileformat(ii)) == TRIM(c_format_h5_native_smss) .OR. &
                 TRIM(infileformat(ii)) == TRIM(c_format_h5_native_mmss) ) ) THEN

            ! 1a) Check for DWD single-sweep multi-moment hdf5 files of filename-convention
            !     "rasXX-pcpng01_sweeph5allm_any_00-2020020300053400-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_00-2020020300055800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_01-2020020300062100-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_02-2020020300064500-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_03-2020020300070800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5allm_any_04-2020020300073100-neu-10557-hd5"
            !
            ! and, in one go,
            !
            ! 1b) Check for DWD single-sweep single-moment hdf5 files of filename-convention
            !     "rasXX-pcpng01_sweeph5onem_dbzh_00-2020020300053400-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_00-2020020300055800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_01-2020020300062100-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_02-2020020300064400-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_03-2020020300070800-neu-10557-hd5"
            !     "rasXX-vol5minng01_sweeph5onem_dbzh_04-2020020300073100-neu-10557-hd5"

            IF ( TRIM(infileformat(ii)) == TRIM(c_format_h5_native_smss) ) THEN
              IF (      TRIM(infilevarname(ii)) == TRIM(obs_dwd_hdf5_varname_vrad)      ) THEN
                fid = i_vrad
              ELSE IF ( TRIM(infilevarname(ii)) == TRIM(obs_dwd_hdf5_varname_dbzh)      ) THEN
                fid = i_dbzh
              ! Using INDEX below (ie substring-containment instead of string equality as above)
              !  allows the user to select the actual version of a moment to be used by
              !  providing the desired ones in the yinradar-folder.
              !  (eg "zdr" or "uzdr" from raw data or "attcorrzdrcorr" from POLARA or even a mix
              !   of these - if multiple data are available for the same moment/time/elevation/...
              !   an error is thrown later on)
              ELSE IF ( INDEX(infilevarname(ii), TRIM(obs_dwd_hdf5_varname_zdr))   /= 0 ) THEN
                fid = i_zdr
              ELSE IF ( INDEX(infilevarname(ii), TRIM(obs_dwd_hdf5_varname_kdp))   /= 0 ) THEN
                fid = i_kdp
              ELSE IF ( INDEX(infilevarname(ii), TRIM(obs_dwd_hdf5_varname_phidp)) /= 0 ) THEN
                fid = i_phidp
              ELSE IF ( INDEX(infilevarname(ii), TRIM(obs_dwd_hdf5_varname_rhv))   /= 0 ) THEN
                fid = i_rhv
              ELSE
                WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // '(): datakind=', TRIM(infilevarname(ii)), &
                     ' not implemented! File discarded! '//TRIM(infilenames(ii))
                CYCLE read_loop_0
              END IF
            ELSE
              fid = -1
            END IF

            ! .. read meta data for the entire volscan or precip scan, for either single-moment files or multi-moment files:
            CALL get_metadata_from_h5_dwd ( TRIM(infilenames(ii)), fid, rsm_dum_onetime, rs_meta_network, ierr )

            IF (ierr == 0) THEN
              have_read_a_file(ipe_rad) = .TRUE.
            ELSE
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): ' // &
                   'file(s) "' // TRIM(infilenames(ii)) // &
                   '" has/have suspicious meta data and is/are discarded!'
            END IF

            rsm_dum_onetime%obsfile_format = TRIM(infileformat(ii))
            ! rsm_dum_onetime%obsfile was set above by get_metadata_from_h5_dwd ()
            
            CYCLE read_loop_0

          END IF

#endif
          ! 2) Check for DWD CDFIN-files of filename-convention
          !     "cdfin_z_id-010908_201307282200_201307282255_volscan" or
          !     "cdfin_z_id-010908_201307282200_201307282255_precipscan"

          IF ( infilecountry(ii) == i_dwd .AND. TRIM(infileformat(ii)) == TRIM(c_format_cdfin) ) THEN
          
            ! .. Read complete metadata from the file:

            SELECT CASE (TRIM(infilevarname(ii)))
            CASE ('vr', 'vrsim')
              fid = i_vrad
            CASE ('qv', 'qvobs')
              fid = i_qualvrad
            CASE ('z',  'zrsim')
              fid = i_dbzh
            CASE ('qz', 'qzobs')
              fid = i_qualdbzh
            CASE default
              WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // '(): datakind=', TRIM(infilevarname(ii)), &
                   ' not implemented! File discarded! '//TRIM(infilenames(ii))
              CYCLE read_loop_0
            END SELECT

            
            CALL get_metadata_from_cdfin ( TRIM(infilenames(ii)), fid, rsm_dum_multitimes, rs_meta_network, ierr )

            ! NOTE: if data are from an OSE / OSSE, rsm_dum_multitimes%icountry will not necessarily be i_dwd,
            !       because get_metadata_from_cdfin() has set it according to the background list for the respective
            !       station_id!
            
            IF (ierr == 0) THEN
              have_read_a_file(ipe_rad) = .TRUE.
            ELSE
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): ' // &
                   'file(s) "' // TRIM(infilenames(ii)) // &
                   '" has/have suspicious meta data and is/are discarded!'
            END IF

            rsm_dum_multitimes%obsfile_format = TRIM(infileformat(ii))
            ! rsm_dum_multitimes%obsfile was set above by get_metadata_from_cdfin ()

            CYCLE read_loop_0

          END IF


#ifdef HDF5_RADAR_INPUT

          ! 3) Check for Swiss files of filename-convention
          !     "PLA1520500007U.001.V.h5"
          !     "PLA1520500007U.001.Z.h5"
          !     "PLA1520500007U.002.V.h5"
          !     "PLA1520500007U.002.Z.h5"

          IF ( infilecountry(ii) == i_meteoswiss .AND. &
               ( TRIM(infileformat(ii)) == TRIM(c_format_h5_native_smss ) .OR. &
                 TRIM(infileformat(ii)) == TRIM(c_format_h5_native_mmss ) ) )  THEN
         
            ! .. Read meta data for the entire volscan using the generalized OPERA reader.
            !    This reader can work with the OPERA filenames "T_PA*" and the MCH filenames "PLA" / "MLA":
            !    NOTE: this reader is capable of multi/single-moment multi/sinlge-sweep files (any combination)
            !          but expects /datasetXX/DATA../ to reflect different elevations and /dataset../dataYY to reflect different moments.
            CALL get_metadata_from_h5_opera ( TRIM(infilenames(ii)), rsm_dum_onetime, rs_meta_network, format_from_file, ierr )

            IF (ierr /= 0) THEN
              
              ! .. the new generalized opera reader had a problem, so try again with the old native reader:

              IF ( TRIM(infileformat(ii)) == TRIM(c_format_h5_native_smss) ) THEN
                
                SELECT CASE (TRIM(infilevarname(ii)))
                CASE ('V')
                  fid = i_vrad
                CASE ('Z')
                  fid = i_dbzh
                CASE default
                  WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // '(): datakind=', TRIM(infilevarname(ii)), &
                       ' not implemented! File discarded! '//TRIM(infilenames(ii))
                  CYCLE read_loop_0
                END SELECT

                CALL get_metadata_from_h5_mch ( TRIM(infilenames(ii)), fid, rsm_dum_onetime, rs_meta_network, ierr )

                IF (ierr == 0) THEN
                  have_read_a_file(ipe_rad) = .TRUE.
                ELSE
                  WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // ' get_metadata_from_h5_mch(): ' // &
                       'file(s) "' // TRIM(infilenames(ii)) // &
                       '" has/have suspicious meta data and is/are discarded!'
                END IF

              ELSE

                WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // ' get_metadata_from_h5_opera(): ' // &
                     'file(s) "' // TRIM(infilenames(ii)) // &
                     '" has/have suspicious meta data and is/are discarded!'
                
              END IF
              
            ELSE IF (TRIM(infileformat(ii)) == TRIM(format_from_file)) THEN
              
              have_read_a_file(ipe_rad) = .TRUE.

            ELSE
              
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // ' get_metadata_from_h5_opera(): ' // &
                   'file(s) "' // TRIM(infilenames(ii)) // &
                   '" has/have suspicious meta data and is/are discarded!'
            END IF


            rsm_dum_onetime%obsfile_format = TRIM(infileformat(ii))
            ! rsm_dum_onetime%obsfile was set above by get_metadata_from_h5_mch ()

            CYCLE read_loop_0

          END IF
#endif

#ifdef HDF5_RADAR_INPUT

          ! 4) Check for OPERA files which are distributed to OPERA from NWS's and which have filename-convention
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

          IF ( (infilecountry(ii) == i_dwd        .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_smss) ) .OR. &
               (infilecountry(ii) == i_meteoswiss .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_mmss) ) .OR. &
               (infilecountry(ii) == i_france     .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_mmss) ) .OR. &
               (infilecountry(ii) == i_belgium    .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
               (infilecountry(ii) == i_poland     .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
               (infilecountry(ii) == i_czech      .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
               (infilecountry(ii) == i_slovakia   .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
               (infilecountry(ii) == i_dwd        .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_smms) ) .OR. &
               (infilecountry(ii) == i_netherlands.AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_mmms) ) .OR. &
               (infilecountry(ii) == i_denmark    .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_2_opera_mmms) )      &
               ) THEN

            ! .. read meta data for the entire volscan, use all files with same first and third word as infilenames(ii)
            !    NOTE: this reader is capable of multi/single-moment multi/sinlge-sweep files (any combination)
            !          but expects /datasetXX/DATA../ to reflect different elevations and /dataset../dataYY to reflect different moments:
            CALL get_metadata_from_h5_opera ( TRIM(infilenames(ii)), rsm_dum_onetime, rs_meta_network, format_from_file, ierr )

            IF (TRIM(infileformat(ii)) == TRIM(format_from_file)) THEN
              IF (ierr == 0) THEN
                have_read_a_file(ipe_rad) = .TRUE.
              ELSE
                WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): ' // &
                     'file(s) "' // TRIM(infilenames(ii)) // &
                     '" has/have suspicious meta data and is/are discarded!'
              END IF
            ELSE
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): ' // &
                   'file TYPE from filename "' // TRIM(infileformat(ii)) // &
                   '" does not match file TYPE from content "' // TRIM(format_from_file) // &
                   '"! Discarded file ' // TRIM(infilenames(ii))
            END IF
            
            CYCLE read_loop_0            

          END IF
#endif

#ifdef HDF5_RADAR_INPUT

          ! 5) Check for Italian files of filename-convention
          !     "odim_201410090010_16101"<
          !     "odim_201410090010_16998"

          IF ( infilecountry(ii) == i_arpasim .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_native_mmms) ) THEN

            ! .. read meta data for the entire volscan, use all files with same first and third word as infilenames(ii)
            CALL get_metadata_from_h5_italy ( TRIM(infilenames(ii)), rsm_dum_onetime, rs_meta_network, ierr )

            IF (ierr == 0) THEN
              have_read_a_file(ipe_rad) = .TRUE.
            ELSE
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): ' // &
                   'file(s) "' // TRIM(infilenames(ii)) // &
                   '" has/have suspicious meta data and is/are discarded!'
            END IF

            rsm_dum_onetime%obsfile_format = TRIM(infileformat(ii))
            ! rsm_dum_onetime%obsfile was set above by get_metadata_from_h5_italy()

            CYCLE read_loop_0

          END IF
#endif

#ifdef HDF5_RADAR_INPUT

          ! 6) Check for KIT C-band files of filename-convention:
          !     "scan-sidpol-120km-14_20001_20230712123503_00.h5"
          !     "scan-sidpol-120km-14_20001_20230712124001_00.h5"
          !
          ! The station is listed under i_dwd and produces mmms files, which the original DWD-radars don't
          
          IF ( infilecountry(ii) == i_dwd .AND. TRIM(infileformat(ii)) == TRIM(c_format_h5_native_mmms) ) THEN

            ! .. read meta data for the entire volscan, use all files with same first and third word as infilenames(ii)
            CALL get_metadata_from_h5_kitcband ( TRIM(infilenames(ii)), rsm_dum_onetime, rs_meta_network, ierr )

            IF (ierr == 0) THEN
              have_read_a_file(ipe_rad) = .TRUE.
            ELSE
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): ' // &
                   'file(s) "' // TRIM(infilenames(ii)) // &
                   '" has/have suspicious meta data and is/are discarded!'
            END IF

            rsm_dum_onetime%obsfile_format = TRIM(infileformat(ii))
            ! rsm_dum_onetime%obsfile was set above by get_metadata_from_h5_italy()

            CYCLE read_loop_0

          END IF
#endif

          ! 7) Possible follow-on checks for other countries:
          !
          !     ...


          ! .. If the loop reached here, we have an unknown file type:
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine) // '(): Unknown file tpye '//TRIM(infilenames(ii))

        END IF   ! my_radar_id == ipe_rad

      END DO read_loop_0

      IF (num_radar > 1) THEN

        distrib_loop_0: DO ii = ii_par, ii_end

          ! determine the number of the PE for output of the radar with no. ista:
          ipe_rad = MOD(ii-1, num_radar)

          ! collect all ista metadata on PE 0:
          CALL distribute_values_radar (have_read_a_file(ipe_rad), 1, ipe_rad, icomm_radar, mpierr)

          IF (have_read_a_file(ipe_rad)) THEN
            IF (has_multiple_obstimes(ii)) THEN
              ista_multi = ista_multi + 1
              IF (my_radar_id == ipe_rad .AND. my_radar_id == 0) THEN
                rs_meta_buf_multitimes(ista_multi) = rsm_dum_multitimes
              END IF
              IF (my_radar_id == ipe_rad .AND. my_radar_id /= 0) THEN
                CALL mpi_send ( rsm_dum_multitimes, 1, mpi_radar_meta_typ_alltimes, &
                     0, ii+5, icomm_radar, mpierr )
              END IF
              IF (my_radar_id /= ipe_rad .AND. my_radar_id == 0) THEN
                CALL mpi_recv ( rs_meta_buf_multitimes(ista_multi), 1, mpi_radar_meta_typ_alltimes, &
                     ipe_rad, ii+5, icomm_radar, mpistat, mpierr )
              END IF
            ELSE
              ista_one = ista_one + 1
              IF (my_radar_id == ipe_rad .AND. my_radar_id == 0) THEN
                rs_meta_buf_onetime(ista_one) = rsm_dum_onetime
              END IF
              IF (my_radar_id == ipe_rad .AND. my_radar_id /= 0) THEN
                CALL mpi_send ( rsm_dum_onetime, 1, mpi_radar_meta_typ_onetime, &
                     0, ii+15, icomm_radar, mpierr )
              END IF
              IF (my_radar_id /= ipe_rad .AND. my_radar_id == 0) THEN
                CALL mpi_recv ( rs_meta_buf_onetime(ista_one), 1, mpi_radar_meta_typ_onetime, &
                     ipe_rad, ii+15, icomm_radar, mpistat, mpierr )
              END IF
            END IF
          END IF

        END DO distrib_loop_0

      ELSE

        IF (have_read_a_file(0)) THEN
          IF (has_multiple_obstimes(ii_par)) THEN
            ista_multi = ista_multi + 1
            rs_meta_buf_multitimes(ista_multi) = rsm_dum_multitimes
          ELSE
            ista_one = ista_one + 1
            rs_meta_buf_onetime(ista_one) = rsm_dum_onetime
          END IF
        END IF
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_ini)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_barrier)
#endif

      ! .. A communication synchronisation might be necessary here:
      IF (num_radar > 1) THEN
        CALL mpi_barrier (icomm_radar, mpierr)
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_barrier)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_ini)
#endif

    END DO collect_loop

    ! =========================================================================
    !
    ! .. Now PE 0 has the full list of meta data from every file.
    !    There can be more than one file for each radarID, therefore merge
    !    and shrink the list now to the known radarIDs. Check also, if files
    !    for different datakinds have the same times and records.
    !
    ! .. Start by looking at all the files with multiple observation times:

    IF (my_radar_id == 0) THEN

      nradsta_multi = 0
      merge_loop_1: DO ii = 1, ista_multi

        ! .. Find out if the radarID from the actual station ii already exists in the previous radar list:
        ista = missval_int
        DO i=1, nradsta_multi
          IF ( rs_meta_buf_multitimes(i)%station_id     == rs_meta_buf_multitimes(ii)%station_id     .AND. &
               TRIM(rs_meta_buf_multitimes(i)%scanname) == TRIM(rs_meta_buf_multitimes(ii)%scanname) ) THEN
            ista = i
            EXIT
          END IF
        END DO

        IF ( ista < 0 ) THEN
          ! .. If the actual radarID does not yet exist, append it to the list of radar stations:
          nradsta_multi = nradsta_multi + 1
          IF (nradsta_multi <= nradsta_max) THEN
            rs_meta_buf_multitimes(nradsta_multi) = rs_meta_buf_multitimes(ii)
          ELSE
            WRITE (*,'(a,i0,a,i0,a)') 'ERROR '//TRIM(yzroutine) // ': nradsta_multi=', nradsta_multi, &
                 ' is larger than nradsta_max=', nradsta_max,'!'
            err_ista(:) = 2
            nradsta_multi = nradsta_multi - 1
          END IF
        ELSE
          ! .. If it does exist, merge the metadata from rs_meta_buf(ii) into rs_meta_buf(ista),
          !     including some checks:
          CALL merge_metadata( rs_meta_buf_multitimes(ii), rs_meta_buf_multitimes(ista), ierr )
          IF (ierr /= 0) THEN
            err_ista(ista) = ierr
          END IF
        END IF

      END DO merge_loop_1

      check_loop_1: DO ii = 1, nradsta_multi

        ! Check if all the needed files are present for each observation time,
        !  and keep only those observation times:
        !--------------------------------------------------------------------
        n = 0
        missing_obs_time = .FALSE.
        DO i=1, rs_meta_buf_multitimes(ii)%nobs_times_obs

          ! Logical mask corresponding to vr, qv, z, qz, etc., to check for internal consistency of the meta data in each present
          !  file, if the file is needed for the operator. The inverse check, namely if some obs data are missing for the actual
          !  configuration, is done at other places in the code. Each process which needs certain radar moments at the same time checks its
          !  pre-requisites itself.
          ! If an element is .TRUE., the corresponding datakind has to be checked:
          checkflags(:)          = .FALSE.
          checkflags(i_vrad)     = loutradwind
          checkflags(i_qualvrad) = lqc_flag
          checkflags(i_dbzh)     = ( loutdbz .OR. (loutradwind .AND. itype_supobing > 0 .AND. itype_obserr_vr > 0) )
          checkflags(i_qualdbzh) = lqc_flag
          checkflags(i_zdr)      = ( loutpolstd .OR. loutpolall )
          checkflags(i_kdp)      = ( loutpolstd .OR. loutpolall )
          checkflags(i_phidp)    = ( loutpolstd .OR. loutpolall )
          checkflags(i_rhv)      = ( loutpolstd .OR. loutpolall )
          checkflagtext(:)(:)       = ' '
          checkflagtext(i_vrad)     = 'loutradwind=.TRUE., but '
          checkflagtext(i_qualvrad) = 'lqc_falg=.TRUE., but '
          checkflagtext(i_dbzh)     = 'loutdbz=.TRUE. or (loutradwind=.TRUE. and itype_supobing > 0 and itype_obserr_vr > 0), but '
          checkflagtext(i_qualdbzh) = 'lqc_flag=.TRUE., but '
          checkflagtext(i_zdr)      = 'loutpol=.TRUE., but for ZDR '
          checkflagtext(i_kdp)      = 'loutpol=.TRUE., but for KDP '
          checkflagtext(i_phidp)    = 'loutpol=.TRUE., but for PHIDP '
          checkflagtext(i_rhv)      = 'loutpol=.TRUE., but for RHOHV '

          ! .. Check for each datakind and present obs file if
          !     there is obs_startrec, obs_endrec and naz_ncdf for all the obs_times_obs.
          !     If something is missing, do not store the obs_time:
          good = .TRUE.
          DO j=1, ndatakind
            ! In merge_metadata, lobs_avail is (temporarily) used to track time-parameter-pairs
            !   where multiple, hence ambiguous, obsdata files are found. Here, we switch off these
            !   data entries by re-setting obsfile as missing and set back lobs_avail to its initial
            !   state (to be used later on for its original purpose).
            ! JM: For multi-time data (as handled here), the ambiguity tracking in merge_metadata is
            !   not working properly (and not required so far), hence switched off for now.
            !   Therefore, so far no need to further handle tracked data here.
            !IF (rs_meta_buf_multitimes(ii)%lobs_avail(i,j)) THEN
            !  rs_meta_buf_multitimes(ii)%obsfile(i,j) = TRIM(obsfile_missingname)
            !  rs_meta_buf_multitimes(ii)%lobs_avail(i,j) = .FALSE.
            !END IF
            IF (checkflags(j)) THEN
              IF (TRIM(rs_meta_buf_multitimes(ii)%obsfile(i,j)) /= obsfile_missingname) THEN
                IF (rs_meta_buf_multitimes(ii)%obs_startrec(i,j) <= 0) THEN
                  good = .FALSE.
                  WRITE (*,'(a)') 'ERROR: '//TRIM(checkflagtext(j))//' wrong obs_startrec in obs data file ' // &
                       rs_meta_buf_multitimes(ii)%obsfile(i,j)
                ELSE IF (rs_meta_buf_multitimes(ii)%obs_endrec(i,j) <= 0) THEN
                  good = .FALSE.
                  WRITE (*,'(a)') 'ERROR: '//TRIM(checkflagtext(j))//' wrong obs_endrec in obs data file ' // &
                       rs_meta_buf_multitimes(ii)%obsfile(i,j)
                ELSE IF (rs_meta_buf_multitimes(ii)%naz_ncdf(i,j) <= 0) THEN
                  good = .FALSE.
                  WRITE (*,'(a)') 'ERROR: '//TRIM(checkflagtext(j))//' wrong naz_ncdf in obs data file ' // &
                       rs_meta_buf_multitimes(ii)%obsfile(i,j)
                END IF
              END IF
            END IF
          END DO
          IF (good) THEN
            n = n + 1
            rs_meta_buf_multitimes(ii)%obs_times_obs(n)  = rs_meta_buf_multitimes(ii)%obs_times_obs(i)
            rs_meta_buf_multitimes(ii)%obs_cdate(n)      = rs_meta_buf_multitimes(ii)%obs_cdate(i)
            rs_meta_buf_multitimes(ii)%obs_startrec(n,:) = rs_meta_buf_multitimes(ii)%obs_startrec(i,:)
            rs_meta_buf_multitimes(ii)%obs_endrec(n,:)   = rs_meta_buf_multitimes(ii)%obs_endrec(i,:)
            rs_meta_buf_multitimes(ii)%obsfile(n,:)      = rs_meta_buf_multitimes(ii)%obsfile(i,:)
            rs_meta_buf_multitimes(ii)%obsfile_format(n) = rs_meta_buf_multitimes(ii)%obsfile_format(i)
            rs_meta_buf_multitimes(ii)%naz_ncdf(n,:)     = rs_meta_buf_multitimes(ii)%naz_ncdf(i,:)
          ELSE
            WRITE (*,'(1x,a,i7,a,/,a)')  &
                 'WARNING '//TRIM(yzroutine)//'(): obs_time ' // &
                 TRIM(rs_meta_buf_multitimes(ii)%obs_cdate(i)) // &
                 ' for radar station ', rs_meta_buf_multitimes(ii)%station_id, &
                 ' discarded because of missing datafile or obs_startrec or obs_endrec or naz_ncdf!', &
                 'Please cross-check obs data files and loutradwind, loutdbz/itype_obserr_vr and lqc_flags, see previous ERROR(s)!'
            missing_obs_time = .TRUE.
          END IF

        END DO
        rs_meta_buf_multitimes(ii)%nobs_times_obs = n

        DO i=n+1, nobstimes_max
          ! Initialize metadata for not needed obs_times_obs with missing values in rs_meta:
          rs_meta_buf_multitimes(ii)%obs_times_obs(i)  = unused_value
          rs_meta_buf_multitimes(ii)%obs_cdate(i)      = 'YYYYMMDDHHMMSS'
          rs_meta_buf_multitimes(ii)%obs_startrec(i,:) = missval_int
          rs_meta_buf_multitimes(ii)%obs_endrec(i,:)   = missval_int
          rs_meta_buf_multitimes(ii)%obsfile(i,:)      = obsfile_missingname
          rs_meta_buf_multitimes(ii)%obsfile_format(i) = ' '
          rs_meta_buf_multitimes(ii)%naz_ncdf(i,:)     = 0
        END DO
        rs_meta_ncdf(ii) = rs_meta_buf_multitimes(ii)
        IF (missing_obs_time) THEN
          err_ista(ii) = 4
        END IF

      END DO check_loop_1

      ! Here nradsta_multi is the final number of radar stations with multitime-files.
      !
      ! Now append the radar stations with one-time files. For this,
      !  use the final multitime-type rs_meta_ncdf for intermediate and
      !  final storage, starting from index nradsta_multi+1:

      nradsta_one = 0
      merge_loop_2: DO ii = 1, ista_one

        ! .. Find out if the radarID from the actual station ii already exists in the previous radar list (onetime-files):
        ista = missval_int
        DO i = nradsta_multi+1, nradsta_multi+nradsta_one
          IF ( rs_meta_ncdf(i)%station_id     == rs_meta_buf_onetime(ii)%station_id     .AND. &
               TRIM(rs_meta_ncdf(i)%scanname) == TRIM(rs_meta_buf_onetime(ii)%scanname) ) THEN
            ista = i
            EXIT
          END IF
        END DO

        IF ( ista < 0 ) THEN
          ! .. If the actual radarID does not yet exist, append it to the list of radar stations:
          nradsta_one = nradsta_one + 1
          IF (nradsta_multi+nradsta_one <= nradsta_max) THEN
            rs_meta_ncdf(nradsta_multi+nradsta_one) = rsm_onetime2multitime( rs_meta_buf_onetime(ii) )
          ELSE
            IF (my_radar_id == 0) THEN
              WRITE (*,'(a,i0,a,i0,a)') 'ERROR '//TRIM(yzroutine) // ': nradsta_multi+nradsta_one = ', &
                   nradsta_multi+nradsta_one, ' is larger than nradsta_max = ', nradsta_max,'!'
            END IF
            err_ista(:) = 2
            nradsta_one = nradsta_one - 1
          END IF
        ELSE
          ! .. If it does exist, merge the metadata from rs_meta_buf(ii) into rs_meta_buf(ista),
          !     including some checks:
          CALL merge_metadata( rs_meta_buf_onetime(ii), rs_meta_ncdf(ista), ierr )
          IF (ierr /= 0) THEN
            err_ista(ista) = ierr
          END IF
        END IF

      END DO merge_loop_2

      check_loop_2: DO ii = nradsta_multi + 1, nradsta_multi + nradsta_one

        ! Check if all the needed files are present for each observation time,
        !  and keep only those observation times:
        !--------------------------------------------------------------------
        n = 0
        missing_obs_time = .FALSE.

        DO i=1, rs_meta_ncdf(ii)%nobs_times_obs

          ! Logical mask corresponding to vr, qv, z, qz, etc., to check for internal consistency of the meta data in each present
          !  file, if the file is needed for the operator. The inverse check, namely if some obs data are missing for the actual
          !  configuration, is done at other places in the code. Each process which needs certain radar moments at the same time checks its
          !  pre-requisites itself.
          ! If an element is .TRUE., the corresponding datakind has to be checked:
          checkflags(:)          = .FALSE.
          checkflags(i_vrad)     = loutradwind
          checkflags(i_qualvrad) = lqc_flag
          checkflags(i_dbzh)     = ( loutdbz .OR. (loutradwind .AND. itype_supobing > 0 .AND. itype_obserr_vr > 0) )
          checkflags(i_qualdbzh) = lqc_flag
          checkflags(i_zdr)      = ( loutpolstd .OR. loutpolall )
          checkflags(i_kdp)      = ( loutpolstd .OR. loutpolall )
          checkflags(i_phidp)    = ( loutpolstd .OR. loutpolall )
          checkflags(i_rhv)      = ( loutpolstd .OR. loutpolall )
          checkflagtext(:)(:)       = ' '
          checkflagtext(i_vrad)     = 'loutradwind=.TRUE., but '
          checkflagtext(i_qualvrad) = 'lqc_falg=.TRUE., but '
          checkflagtext(i_dbzh)     = 'loutdbz=.TRUE. or (loutradwind=.TRUE. and itype_supobing > 0 and itype_obserr_vr > 0), but '
          checkflagtext(i_qualdbzh) = 'lqc_flag=.TRUE., but '
          checkflagtext(i_zdr)      = 'loutpol=.TRUE., but for ZDR '
          checkflagtext(i_kdp)      = 'loutpol=.TRUE., but for KDP '
          checkflagtext(i_phidp)    = 'loutpol=.TRUE., but for PHIDP '
          checkflagtext(i_rhv)      = 'loutpol=.TRUE., but for RHOHV '

          ! .. Check for each datakind and present obs file if
          !     there is obs_startrec, obs_endrec and naz_ncdf for all the obs_times_obs.
          !     If something is missing, do not store the obs_time:
          good = .TRUE.
          DO j=1, ndatakind
            ! In merge_metadata, lobs_avail is (temporarily) used to track time-parameter-pairs
            !   where multiple, hence ambiguous, obsdata files are found. Here, we switch off these
            !   data entries by re-setting obsfile as missing and set back lobs_avail to its initial
            !   state (to be used later on for its original purpose).
            IF (rs_meta_ncdf(ii)%lobs_avail(i,j)) THEN
              rs_meta_ncdf(ii)%obsfile(i,j) = TRIM(obsfile_missingname)
              rs_meta_ncdf(ii)%lobs_avail(i,j) = .FALSE.
            END IF
            IF (checkflags(j)) THEN
              IF (TRIM(rs_meta_ncdf(ii)%obsfile(i,j)) /= obsfile_missingname) THEN
                IF (rs_meta_ncdf(ii)%obs_startrec(i,j) <= 0) THEN
                  good = .FALSE.
                  WRITE (*,'(a)') 'ERROR: '//TRIM(checkflagtext(j))//' wrong obs_startrec in obs data file ' // &
                       rs_meta_ncdf(ii)%obsfile(i,j)
                ELSE IF (rs_meta_ncdf(ii)%obs_endrec(i,j) <= 0) THEN
                  good = .FALSE.
                  WRITE (*,'(a)') 'ERROR: '//TRIM(checkflagtext(j))//' wrong obs_endrec in obs data file '
                ELSE IF (rs_meta_ncdf(ii)%naz_ncdf(i,j) <= 0) THEN
                  good = .FALSE.
                  WRITE (*,'(a)') 'ERROR: '//TRIM(checkflagtext(j))//' wrong naz_ncdf in obs data file '
                END IF
              END IF
            END IF
          END DO
          IF (good) THEN
            n = n + 1
            rs_meta_ncdf(ii)%obs_times_obs(n)  = rs_meta_ncdf(ii)%obs_times_obs(i)
            rs_meta_ncdf(ii)%obs_cdate(n)      = rs_meta_ncdf(ii)%obs_cdate(i)
            rs_meta_ncdf(ii)%obs_startrec(n,:) = rs_meta_ncdf(ii)%obs_startrec(i,:)
            rs_meta_ncdf(ii)%obs_endrec(n,:)   = rs_meta_ncdf(ii)%obs_endrec(i,:)
            rs_meta_ncdf(ii)%obsfile(n,:)      = rs_meta_ncdf(ii)%obsfile(i,:)
            rs_meta_ncdf(ii)%obsfile_format(n) = rs_meta_ncdf(ii)%obsfile_format(i)
            rs_meta_ncdf(ii)%naz_ncdf(n,:)     = rs_meta_ncdf(ii)%naz_ncdf(i,:)
          ELSE
            WRITE (*,'(1x,a,i7,a,/,a)')  &
                 'WARNING '//TRIM(yzroutine)//'(): obs_time ' // &
                 TRIM(rs_meta_ncdf(ii)%obs_cdate(i)) // &
                 ' for radar station ', rs_meta_ncdf(ii)%station_id, &
                 ' discarded because of missing datafile or obs_startrec or obs_endrec or naz_ncdf!', &
                 'Please cross-check obs data files and loutradwind, loutdbz/itype_obserr_vr and lqc_flags, see previous ERROR(s)!'
            missing_obs_time = .TRUE.
          END IF

        END DO
        rs_meta_ncdf(ii)%nobs_times_obs = n

        DO i=n+1, nobstimes_max
          ! Initialize metadata for not needed obs_times and obs_times_obs with missing values in rs_meta:
          rs_meta_ncdf(ii)%obs_times(i)      = unused_value
          rs_meta_ncdf(ii)%obs_times_obs(i)  = unused_value
          rs_meta_ncdf(ii)%obs_cdate(i)      = 'YYYYMMDDHHMMSS'
          rs_meta_ncdf(ii)%obs_startrec(i,:) = missval_int
          rs_meta_ncdf(ii)%obs_endrec(i,:)   = missval_int
          rs_meta_ncdf(ii)%obsfile(i,:)      = obsfile_missingname
          rs_meta_ncdf(ii)%obsfile_format(i) = ' '
          rs_meta_ncdf(ii)%naz_ncdf(i,:)     = 0
        END DO

        IF (missing_obs_time) THEN
          err_ista(ii) = 4
        END IF

      END DO check_loop_2
      
      nradsta_ncdf = nradsta_multi + nradsta_one

      ! .. Now nradsta is the actual number of stations and rs_meta_ncdf(1:nradsta) is the list
      !     of stations which have all the needed data kinds.
      !

      ! .. Store the obs times from the files in %obs_times for proper bookkeeping:
      DO ii=1, nradsta_ncdf
          rs_meta_ncdf(ii)%obs_times(:) = rs_meta_ncdf(ii)%obs_times_obs(:)
          rs_meta_ncdf(ii)%nobs_times   = rs_meta_ncdf(ii)%nobs_times_obs
      END DO

      !
      ! =========================================================================

      IF (nradsta_ncdf < 1) THEN
        err_ista = 1
        WRITE(*,'(2a)') TRIM(yzroutine) // ' ERROR: No consistent sets (vr,qv,z,qz,...) of required radar station files ' &
             //'could be found in directory ''', TRIM(ydirradarin)//''''
      END IF

    END IF   ! my_radar_id == 0

    IF (num_radar > 1) THEN
      CALL distribute_values_radar (nradsta_ncdf, 1, 0, icomm_radar, ierr)
    END IF

!!$    ! Distribute error-status:
!!$    errstring(:) = ' '
!!$    CALL global_values_radar(err_ista, nradsta_ncdf, "MAX", icomm_radar, 0, errstring, mpierr)

    ! Clean up:
    DEALLOCATE (rs_meta_buf_multitimes, rs_meta_buf_onetime)
    IF (ALLOCATED(infilenames)          ) DEALLOCATE(infilenames)
    IF (ALLOCATED(infileformat)         ) DEALLOCATE(infileformat)
    IF (ALLOCATED(infilevarname)        ) DEALLOCATE(infilevarname)
    IF (ALLOCATED(infilecountry)        ) DEALLOCATE(infilecountry)
    IF (ALLOCATED(has_multiple_obstimes)) DEALLOCATE(has_multiple_obstimes)
    IF (ALLOCATED(words)                ) DEALLOCATE(words)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE read_meta_info_all

  !==============================================================================

  SUBROUTINE get_scanstrategy_ID_m ( rsm, cfileident, manel, tol_ele, err )

    !------------------------------------------------------------------------------
    !
    ! Description: Determines the elevations given in the vector manel(:)
    !              match any nominal scan strategy of the radar station
    !              with the metadata given in rsm.
    !
    !              If manel matches any of the nominal/default scan
    !              strategies in rsm%el_arr_default(:,:), the scan strategy
    !              parameters rsm%el_arr(:) and rsm%nel are updated accordingly,
    !              and rsm%scanname is set.
    !
    !              cfileident is a string to identify the file from which
    !              manel was read, and is used for error messages.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! In, InOut, Outgoing:

    TYPE(radar_meta_type), INTENT(inout) :: rsm
    CHARACTER(LEN=*),      INTENT(in)    :: cfileident
    REAL(KIND=dp),         INTENT(in)    :: manel(:)
    REAL(KIND=dp),         INTENT(in)    :: tol_ele   ! tolerance for elevation angle matching [deg]
    INTEGER,               INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine= 'get_scanstrategy_ID_m'

    INTEGER                    :: i, j, ir, n, irepel, iscstr_def, iscanstrat
    INTEGER                    :: scanstrategy_score(nscanstrategies_max)
    REAL(KIND=dp), ALLOCATABLE :: elev(:)
    CHARACTER (LEN=80)         :: yzerrmsg
    CHARACTER (LEN=3)          :: cisc
    CHARACTER (LEN=3)          :: cirepel

    err = 0

    n = SIZE(manel)
    IF (n <= 0) THEN
      err = 6
      WRITE (*,'(a,i0,a,a,i6.6)') &
           TRIM(yzroutine)//' WARNING: no elevations found in NetCDF file ' // &
           TRIM(cfileident) // ' for station ID ', rsm%station_id
      RETURN
    END IF

    ALLOCATE(elev(n))

    irepel=1
    elev(1)    = manel(1)
    DO i = 2, n
      !.. Find and eliminate duplicated elevations for rsm%el_arr
      !   and sort elevations in ascending direction:
      DO ir = 1, irepel
        IF (manel(i) == elev(ir)) EXIT
        IF (manel(i) < elev(ir)) THEN
          ! manel(i) is smaller than the largest elevation so far,
          ! so put it inbetween the correct neighouring elevs:
          irepel   = irepel + 1
          elev(ir+1:irepel) = elev(ir:irepel-1)
          elev(ir) = manel(i)
          EXIT
        END IF
        IF (ir == irepel) THEN
          ! if we reach this point at the end of the loop,
          ! no equal or larger elevation has been found so far,
          ! so this is a new elevation to put at the end of the elev-vector:
          irepel   = irepel + 1
          elev(irepel) = manel(i)
        END IF
      END DO

    END DO

    ! determine which of the default scan strategies matches the
    ! radar data in the netcdf files best:
    scanstrategy_score(:) = 0
    DO ir = 1, nscanstrategies_max
      IF (rsm%nel_default(ir) >= irepel ) THEN

        DO j=1, rsm%nel_default(ir)
          DO i=1, irepel
            IF (ABS(rsm%el_arr_default(j,ir) - elev(i)) <= tol_ele) THEN
              scanstrategy_score(ir) = scanstrategy_score(ir) + 1
            END IF
          END DO
        END DO
        
      END IF
    END DO

    ! if there are candidates where all elevations from the obs file are contained, choose the one
    !  with the lowest number of nominal elevations:
    iscanstrat = -9999
    DO ir = 1, nscanstrategies_max
      IF (scanstrategy_score(ir) == irepel) THEN
        iscanstrat = ir
      END IF
    END DO
    IF (iscanstrat > 0) THEN
      DO ir = 1, nscanstrategies_max
        IF (scanstrategy_score(ir) == irepel .AND. rsm%nel_default(ir) < rsm%nel_default(iscanstrat)) THEN
          iscanstrat = ir
        END IF
      END DO
    END IF

    IF (iscanstrat < 0) THEN
      err = 6
      iscstr_def = 0
      DO ir = 1, nscanstrategies_max
        IF (rsm%nel_default(ir) > 0) THEN
          iscstr_def = iscstr_def + 1
        END IF
      END DO
      cisc(:) = ''
      WRITE(cisc, '(i3)') iscstr_def
      cirepel(:) = ''
      WRITE(cirepel, '(i3)') irepel
      WRITE (*,'(a,i6.6,a,i0,a,'//TRIM(ADJUSTL(cirepel))//'(1x,f7.2),2(a,'//TRIM(ADJUSTL(cisc))//'(i0,1x)))') &
           TRIM(yzroutine)//' WARNING: no suitably matching registered default scan strategy found for file ' // &
           TRIM(cfileident) // ' (station ID=', rsm%station_id,', nele=', irepel, ', elevs=', elev(1:irepel), &
           '). The registered scan strategies have nele=', rsm%nel_default(1:iscstr_def), &
           ' nmatch=', scanstrategy_score(1:iscstr_def)
    ELSE

      rsm%nel                  = rsm%nel_default(iscanstrat)
      rsm%el_arr(1:rsm%nel)    = rsm%el_arr_default(1:rsm%nel, iscanstrat)

      CALL set_scanname ( rsm )

    END IF

    !.. round obs elevations to next multiple of 0.05  degrees:
!!$    DO ir = 1, irepel
!!$      elev(ir) = NINT(elev(ir)/0.05) * 0.05
!!$    END DO

    DEALLOCATE(elev)

  END SUBROUTINE get_scanstrategy_ID_m

  SUBROUTINE get_scanstrategy_ID_o ( rsm, cfileident, manel, tol_ele, err )

    !------------------------------------------------------------------------------
    !
    ! Description: Determines the elevations given in the vector manel(:)
    !              match any nominal scan strategy of the radar station
    !              with the metadata given in rsm.
    !
    !              If manel matches any of the nominal/default scan
    !              strategies in rsm%el_arr_default(:,:), the scan strategy
    !              parameters rsm%el_arr(:) and rsm%nel are updated accordingly,
    !              and rsm%scanname is set.
    !
    !              cfileident is a string to identify the file from which
    !              manel was read, and is used for error messages.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! In, InOut, Outgoing:

    TYPE(radar_meta_type_onetime), INTENT(inout) :: rsm
    CHARACTER(LEN=*),              INTENT(in)    :: cfileident
    REAL(KIND=dp),                 INTENT(in)    :: manel(:)
    REAL(KIND=dp),                 INTENT(in)    :: tol_ele   ! tolerance for elevation angle matching [deg]
    INTEGER,                       INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_scanstrategy_ID_o'

    INTEGER                    :: i, j, ir, n, irepel, iscstr_def, iscanstrat
    INTEGER                    :: scanstrategy_score(nscanstrategies_max)
    REAL(KIND=dp), ALLOCATABLE :: elev(:)
    CHARACTER (LEN=80)         :: yzerrmsg
    CHARACTER (LEN=3)          :: cisc
    CHARACTER (LEN=3)          :: cirepel

    err = 0

    n = SIZE(manel)
    IF (n <= 0) THEN
      err = 6
      WRITE (*,'(a,i0,a,a,i6.6)') &
           TRIM(yzroutine)//' WARNING: no elevations found in NetCDF file ' // &
           TRIM(cfileident) // ' for station ID ', rsm%station_id
      RETURN
    END IF

    ALLOCATE(elev(n))

    irepel=1
    elev(1)    = manel(1)
    DO i = 2, n
      !.. Find and eliminate duplicated elevations for rsm%el_arr
      !   and sort elevations in ascending direction:
      DO ir = 1, irepel
        IF (manel(i) == elev(ir)) EXIT
        IF (manel(i) < elev(ir)) THEN
          ! manel(i) is smaller than the largest elevation so far,
          ! so put it inbetween the correct neighouring elevs:
          irepel   = irepel + 1
          elev(ir+1:irepel) = elev(ir:irepel-1)
          elev(ir) = manel(i)
          EXIT
        END IF
        IF (ir == irepel) THEN
          ! if we reach this point at the end of the loop,
          ! no equal or larger elevation has been found so far,
          ! so this is a new elevation to put at the end of the elev-vector:
          irepel   = irepel + 1
          elev(irepel) = manel(i)
        END IF
      END DO

    END DO

    ! determine which of the default scan strategies matches the
    ! radar data in the netcdf files best:
    scanstrategy_score(:) = 0
    DO ir = 1, nscanstrategies_max
      IF (rsm%nel_default(ir) >= irepel ) THEN

        DO j=1, rsm%nel_default(ir)
          DO i=1, irepel
            IF (ABS(rsm%el_arr_default(j,ir) - elev(i)) <= tol_ele) THEN
              scanstrategy_score(ir) = scanstrategy_score(ir) + 1
            END IF
          END DO
        END DO
        
      END IF
    END DO

    ! if there are candidates where all elevations from the obs file are contained, choose the one
    !  with the lowest number of nominal elevations:
    iscanstrat = -9999
    DO ir = 1, nscanstrategies_max
      IF (scanstrategy_score(ir) == irepel) THEN
        iscanstrat = ir
      END IF
    END DO
    IF (iscanstrat > 0) THEN
      DO ir = 1, nscanstrategies_max
        IF (scanstrategy_score(ir) == irepel .AND. rsm%nel_default(ir) < rsm%nel_default(iscanstrat)) THEN
          iscanstrat = ir
        END IF
      END DO
    END IF

    IF (iscanstrat < 0) THEN
      err = 6
      iscstr_def = 0
      DO ir = 1, nscanstrategies_max
        IF (rsm%nel_default(ir) > 0) THEN
          iscstr_def = iscstr_def + 1
        END IF
      END DO
      cisc(:) = ''
      WRITE(cisc, '(i3)') iscstr_def
      cirepel(:) = ''
      WRITE(cirepel, '(i3)') irepel
      WRITE (*,'(a,i6.6,a,i0,a,'//TRIM(ADJUSTL(cirepel))//'(1x,f7.2),2(a,'//TRIM(ADJUSTL(cisc))//'(i0,1x)))') &
           TRIM(yzroutine)//' WARNING: no suitably matching registered default scan strategy found for file ' // &
           TRIM(cfileident) // ' (station ID=', rsm%station_id,', nele=', irepel, ', elevs=', elev(1:irepel), &
           '). The registered scan strategies have nele=', rsm%nel_default(1:iscstr_def), &
           ' nmatch=', scanstrategy_score(1:iscstr_def)
    ELSE

      rsm%nel                  = rsm%nel_default(iscanstrat)
      rsm%el_arr(1:rsm%nel)    = rsm%el_arr_default(1:rsm%nel, iscanstrat)

      CALL set_scanname ( rsm )

    END IF

    !.. round obs elevations to next multiple of 0.05  degrees:
!!$    DO ir = 1, irepel
!!$      elev(ir) = NINT(elev(ir)/0.05) * 0.05
!!$    END DO

    DEALLOCATE(elev)

  END SUBROUTINE get_scanstrategy_ID_o

  !==============================================================================

  SUBROUTINE get_metadata_from_cdfin ( infilename, fid, rsm, rsm_bg, err )

    TYPE(radar_meta_type),   INTENT(inout) :: rsm
    TYPE(radar_meta_type),   INTENT(in)    :: rsm_bg(:)  ! Background-list for cross-checking
    CHARACTER(LEN=*),        INTENT(in)    :: infilename
    INTEGER,                 INTENT(in)    :: fid
    INTEGER,                 INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_metadata_from_cdfin'

    INTEGER                  :: i, ncid, ierr
    INTEGER                  :: status, RecordDimID, AziDim, RBDim, EleDim
    INTEGER                  :: mii, niii, len_report, len_azimut,len_range,len_ele, &
                                mhp, azres, rngbs, station_id
    REAL (KIND=dp)           :: mlah, mloh, enyq, hnyq, dprf, mrglen, azoff, rgoff
    INTEGER                  :: mngave, mnipul
    ! NetCDF variable id for
    INTEGER                  :: varid_yssosn ! short name of station
    INTEGER                  :: varid_mii    ! WMO-block-number
    INTEGER                  :: varid_niii   ! WMO-station number
    INTEGER                  :: varid_mlah   ! latitude
    INTEGER                  :: varid_mloh   ! longitude
    INTEGER                  :: varid_mhp    ! height of station
    INTEGER                  :: varid_rgoff  ! offset of radial (021203)
    INTEGER                  :: varid_azoff  ! offset of azimuths (021204)
    INTEGER                  :: varid_rngbs  ! range bin size (021201)
    INTEGER                  :: varid_azres  ! azimutal resolution (021202)
    INTEGER                  :: varid_manel  ! antenna elevation
    INTEGER                  :: varid_time   ! reference time
    INTEGER                  :: varid_date   ! reference date
    INTEGER                  :: varid_enyq   ! extended NYQUIST Velocity (021236)
    INTEGER                  :: varid_hnyq   ! high NYQUIST Velocity (021237)
    INTEGER                  :: varid_dprf   ! DualPRF Ratio (002194)
    INTEGER                  :: varid_mrglen ! RANGE-GATE LENGTH
    INTEGER                  :: varid_mngave ! NUMBER OF GATES AVERAGED
    INTEGER                  :: varid_mnipul ! NUMBER OF INTEGRATED PULSES
    INTEGER                  :: varid_mjjj
    INTEGER                  :: varid_mmm
    INTEGER                  :: varid_myy
    INTEGER                  :: varid_mgg
    INTEGER                  :: varid_ngg
    INTEGER                  :: varid_nsecm
    REAL(KIND=dp), ALLOCATABLE :: manel(:)
    INTEGER,       ALLOCATABLE :: date(:), time(:)
    INTEGER,       ALLOCATABLE :: nc_mjjj(:), nc_mmm(:), nc_myy(:)  , &
                                  nc_mgg(:),  nc_ngg(:), nc_nsecm(:)
    CHARACTER (LEN=3)          :: yssosn
    CHARACTER (LEN=80)         :: yzerrmsg
    CHARACTER (LEN=14), ALLOCATABLE       :: cdate(:)
    LOGICAL                    :: station_found

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    err = 0

    status = nf90_open (TRIM(ydirradarin)//TRIM(infilename), NF90_NOWRITE, ncid)

    IF (status /= nf90_noerr) THEN
      WRITE(*,'(i5,1x,3a)') my_radar_id, 'ERROR opening NetCDF file ', TRIM(ydirradarin)//TRIM(infilename)//' ', &
           TRIM(nf90_strerror(status))
      err = 1
      RETURN
    END IF

    status = nf90_inquire          (ncid, unlimitedDimId = RecordDimID)
    status = nf90_inquire_dimension(ncid, RecordDimID,len = len_report) !Number
    status = nf90_inq_dimid        (ncid, "Loop_000_maxlen", AziDim)
      IF (status /= NF90_NOERR) status = nf90_inq_dimid(ncid, "n_azimuth", AziDim)
    status = nf90_inquire_dimension(ncid, AziDim, len = len_azimut) ! Number of azimuth
    status = nf90_inq_dimid        (ncid, "Loop_001_maxlen", RBDim)
      IF (status /= NF90_NOERR) status = nf90_inq_dimid(ncid, "n_range", RBDim)
    status = nf90_inquire_dimension(ncid, RBDim, len = len_range)   ! Number of range bins
    status = nf90_inq_dimid        (ncid, "section1_length", EleDim)
      IF (status /= NF90_NOERR) status = nf90_inq_dimid(ncid, "n_elevation", EleDim)
    status = nf90_inquire_dimension(ncid, EleDim, len = len_ele)   ! Max. number of elevations (can be wrong, so do not use!)

    status = nf90_inq_varid        (ncid, 'MII',   varid_mii)            ! station_id(1:3)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'country_ID',varid_mii)
    status = nf90_get_var          (ncid, varid_mii,  mii,  start=(/1/)) ! Station ID country
    status = nf90_inq_varid        (ncid, 'NIII',  varid_niii)           ! station_id(4:6)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'station_ID_national',varid_niii)
    status = nf90_get_var          (ncid, varid_niii, niii, start=(/1/)) ! Local station ID

    status = nf90_inq_varid        (ncid, 'YSSOSN', varid_yssosn)   ! m%station_name
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'station_name',varid_yssosn)
    ! YSSOSN is missing in the quality flag files, so don't read it if it is not there!
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var          (ncid, varid_yssosn, yssosn, start=(/1/)) ! short name of station
    ELSE
      yssosn = 'XXX'
    ENDIF

    station_id = mii*1000 + niii

    IF ( len_report < 1 ) THEN
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // ' ERROR: NetCDF file ' // &
           TRIM(infilename) // ' for radar station ', station_id, &
           ' do not contain any data records!'
      err = 5
      status = nf90_close(ncid)
      RETURN
    END IF

    ! .. Search the actual station in the background list and use it's meta data
    !     as a background list:
    station_found = .FALSE.
    DO i = 1, SIZE(rsm_bg)
      IF ( station_id == rsm_bg(i)%station_id ) THEN
        ! .. Copy the base values into rs_meta_ncdf(ista).
        rsm = rsm_bg(i)
        station_found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT.station_found) THEN
      err = 3
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // ' WARNING: radar station ', station_id, &
           ' not found in background list rsm_bg!'
      status = nf90_close(ncid)
      RETURN
    END IF

    ! .. Now partially override with the meta data from the file:
    IF (yssosn /= 'XXX') THEN
      rsm%station_name     = yssosn
    END IF
    rsm%station_id       = station_id
    rsm%naz_ncdf(:,fid)  = len_azimut
    rsm%nra              = len_range    ! max. number of range bins occuring in a volume scan
    rsm%nra_obs          = len_range
    rsm%nrep_ncdf(fid)   = len_report

    status = nf90_inq_varid (ncid, 'MLAH',  varid_mlah)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'station_latitude',varid_mlah)
    status = nf90_get_var   (ncid, varid_mlah, mlah, start=(/1/)) ! Latitude of station
    rsm%lat         = mlah
    status = nf90_inq_varid (ncid, 'MLOH',  varid_mloh)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'station_longitude',varid_mloh)
    status = nf90_get_var   (ncid, varid_mloh, mloh, start=(/1/)) ! Longitude of station
    rsm%lon         = mloh
    status = nf90_inq_varid (ncid, 'MHP' ,  varid_mhp)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'station_height',varid_mhp)
    status = nf90_get_var   (ncid, varid_mhp,  mhp,  start=(/1/)) ! Altitude of station
    rsm%alt_msl_true     = mhp
    status = nf90_inq_varid (ncid, '021201',varid_rngbs)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'range_resolution',varid_rngbs)
    status = nf90_get_var   (ncid, varid_rngbs,rngbs,start=(/1/)) ! Range bin resolution in m
    rsm%ra_inc      = rngbs
    rsm%ra_inc_obs  = rngbs
    status = nf90_inq_varid (ncid, '021202',varid_azres)
      IF (status .NE. NF90_NOERR ) status = nf90_inq_varid (ncid, 'azimuthal_resolution',varid_azres)
    status = nf90_get_var   (ncid, varid_azres,azres,start=(/1/)) ! Azimuthal resolution
    rsm%az_inc      = azres

! ra_start has the wrong value of 0.0 in the data files, so do not overtake it!
!   value is taken from the background list!
!!$    status = nf90_inq_varid (ncid, '021203',varid_rgoff)
!!$      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'range_start',varid_rgoff)
!!$    IF (status == NF90_NOERR) THEN
!!$      status = nf90_get_var   (ncid, varid_rgoff,rgoff,start=(/1/)) ! Range bin resolution in m
!!$      rsm%ra_start    = rgoff  ! wrong value of 0.0 in cdfin-files!
!!$    ELSE
!!$      rsm%ra_start    = miss_value
!!$    END IF

! az_start is not necessarily the center of the first output azimut sampling interval, so do not overtake it!
!   value is taken from the background list!
!!$    status = nf90_inq_varid (ncid, '021204',varid_azoff)
!!$      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'azimuth_start',varid_azoff)
!!$    IF (status == NF90_NOERR) THEN
!!$      status = nf90_get_var   (ncid, varid_azoff,azoff,start=(/1/)) ! Starting azimuth
!!$      rsm%az_start    = azoff
!!$    ELSE
!!$      rsm%az_start    = miss_value
!!$    END IF

    status = nf90_inq_varid (ncid, '021236',varid_enyq)    !extended Nyquist Velocity
    IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'extended_nyquist',varid_enyq)
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var   (ncid, varid_enyq, enyq, start=(/1/)) ! extended Nyquist Velocity
      rsm%ext_nyq(:,:)     = enyq
    END IF
    status = nf90_inq_varid (ncid, '021237',varid_hnyq)    !high Nyquist Velocity
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'high_nyquist',varid_hnyq)
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var   (ncid, varid_hnyq, hnyq, start=(/1/)) ! high Nyquist Velocity
      rsm%high_nyq(:)    = hnyq
    END IF
    status = nf90_inq_varid (ncid, '002194',varid_dprf)    !DualPRF ratio
    IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'dualPRF_ratio',varid_dprf)
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var   (ncid, varid_dprf, dprf, start=(/1/)) ! DualPRF ratio
      rsm%dualprf_ratio(:) = dprf
    END IF
    status = nf90_inq_varid (ncid, 'MRGLEN',varid_mrglen)  !RANGE-GATE LENGTH
    IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'range_gate_length',varid_mrglen)
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var   (ncid, varid_mrglen, mrglen, start=(/1/)) ! RANGE-GATE LENGTH
      rsm%rngate_len  = mrglen
    END IF
    status = nf90_inq_varid (ncid, 'MNGAVE',varid_mngave)  !NUMBER OF GATES AVERAGED
    IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'n_ranges_averaged',varid_mngave)
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var   (ncid, varid_mngave, mngave, start=(/1/)) ! NUMBER OF GATES AVERAGED
      rsm%num_gates   = mngave
    END IF
    status = nf90_inq_varid (ncid, 'MNIPUL',varid_mnipul)  !NUMBER OF INTEGRATED PULSES
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'n_pulses_averaged',varid_mnipul)
    IF (status == NF90_NOERR) THEN
      status = nf90_get_var   (ncid, varid_mnipul, mnipul, start=(/1/)) ! NUMBER OF INTEGRATED PULSES
      rsm%num_pulses  = mnipul
    END IF


    ALLOCATE(time    (len_report), &
         date    (len_report), &
         cdate   (len_report), &
         manel   (len_report), &
         nc_mjjj (len_report), &
         nc_mmm  (len_report), &
         nc_myy  (len_report), &
         nc_mgg  (len_report), &
         nc_ngg  (len_report), &
         nc_nsecm(len_report)   )

    status = nf90_inq_varid   (ncid, 'MANEL', varid_manel)
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'ppi_elevation',varid_manel)
    status = nf90_get_var     (ncid, varid_manel, manel ,  (/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'TIME', varid_time)
    status = nf90_get_var     (ncid, varid_time , time    ,(/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'DATE', varid_date)
    status = nf90_get_var     (ncid, varid_date , date    ,(/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'MJJJ',  varid_mjjj)    !Year
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'year',varid_mjjj)
    status = nf90_get_var     (ncid, varid_mjjj , nc_mjjj ,(/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'MMM',   varid_mmm)     !Month
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'month',varid_mmm)
    status = nf90_get_var     (ncid, varid_mmm  , nc_mmm  ,(/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'MYY',   varid_myy)     !Day
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'day',varid_myy)
    status = nf90_get_var     (ncid, varid_myy  , nc_myy  ,(/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'MGG',   varid_mgg)     !Hour
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'hour',varid_mgg)
    status = nf90_get_var     (ncid, varid_mgg  , nc_mgg  ,(/ 1 /), (/ len_report /))
    status = nf90_inq_varid   (ncid, 'NGG',   varid_ngg)     !Minute
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'minute',varid_ngg)
    status = nf90_get_var     (ncid, varid_ngg  , nc_ngg  ,(/ 1 /) , (/ len_report /))
    status = nf90_inq_varid   (ncid, 'NSECM', varid_nsecm)   !Second within a minute
      IF (status /= NF90_NOERR ) status = nf90_inq_varid (ncid, 'second',varid_nsecm)
    status = nf90_get_var     (ncid, varid_nsecm, nc_nsecm,(/ 1 /), (/ len_report /))

    status = nf90_close(ncid)

    ! The time information is the start time of the volume scan. Store it in cdate:
    DO i = 1, len_report
      IF (nc_mjjj(i) < 100) THEN
        nc_mjjj (i) = 2000 + nc_mjjj(i)
        IF (nc_mjjj(i) > 2059)  nc_mjjj (i) = nc_mjjj(i) - 100
      END IF
!!$           WRITE (cdate(i),'(i4.4,5i2.2)') &
!!$                nc_mjjj(i), nc_mmm(i), nc_myy(i), nc_mgg(i), nc_ngg(i), nc_nsecm(i)
!!$!    nsecm is not 0 (as in the nominal time) but the true seconds part of the end time of
!!$!          the scan, so disregard nsecm and take time-field (= "hhmm00") instead:
      WRITE (cdate(i),'(i4.4,2i2.2,i6.6)') &
           nc_mjjj(i), nc_mmm(i), nc_myy(i), time(i)
    END DO

    ! Extract blocks of obstimes from vector cdate(:); one possible error
    !  status (9) of this subroutine (size(cdate)<0)  equals err_ista(ista) = 5
    !  from above and has already been checked. Another error (10) is if the obs times are not sorted
    !  in monotonic ascending order. If this should be the case, the data are most likely corrupt.
    CALL append_obstimes_to_rsm ( rsm, TRIM(infilename), (/ fid /), ydate_ini_mod, cdate, ierr )
    IF (ierr /= 0) THEN
      err = ierr
    END IF

    ! .. recognizes also DWD precip scans:
    CALL get_scanstrategy_ID ( rsm, TRIM(infilename), manel, 0.1_dp, ierr )
    IF (ierr /= 0) THEN
      err = ierr
    END IF

    ! .. Store also the infilename of this file in the rsm structure:
    rsm%obsfile(:,fid) = TRIM(infilename)

    DEALLOCATE(time, date, cdate,      &
         manel, nc_mjjj, nc_mmm, &
         nc_myy, nc_mgg, nc_ngg, &
         nc_nsecm )

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_metadata_from_cdfin

  !==============================================================================

#ifdef HDF5_RADAR_INPUT

  !==============================================================================

  SUBROUTINE get_metadata_from_h5_dwd ( infilename, fid, rsmo, rsmm_bg, err )

    !==============================================================================
    !
    ! .. HDF5-files from DWD: two different file types:
    !    1) one file per station, per variable, per obs_time and per PPI-elevation (10 elevations)
    !    2) one file per station, per obs_time and per PPI-elevation (10 elevations) for all variables
    !
    !    and two different scan strategies:
    !    a) volume scans
    !    b) precipitation scans (1 elevation, horizon-following variable elevation)
    !
    !    This routine reads all the files for all sweeps of one volume or
    !    precipitation scan.
    !    The basename of the files is expected in the string "infilename".
    !    Depending on options a) or b) above:
    !    - for a), the infilename has to provide an appended "|" separated list
    !      of elevation indicators and minute/second strings ("mmss")
    !      needed to reconstruct the individual file names for each sweep.
    !    - for b), no such list is expected, just the plain filename.
    !
    !    The expected value of parameter fid depends on options 1) and 2) above:
    !    - for 1), it has to contain the unique integer index which represents
    !      the quantity in the file (e.g., 1=vr, 3=hor. reflectivity)
    !    - for 2), it has to be an integer < 1
    !
    !==============================================================================

    TYPE(radar_meta_type_onetime), INTENT(inout) :: rsmo
    TYPE(radar_meta_type),         INTENT(in)    :: rsmm_bg(:)  ! Background-list for cross-checking
    CHARACTER(LEN=*),              INTENT(in)    :: infilename
    INTEGER,                       INTENT(in)    :: fid
    INTEGER,                       INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_metadata_from_h5_dwd'

    CHARACTER(len=cobsflen), ALLOCATABLE   :: infilenames(:)
    INTEGER                                :: i, k, kk, ii, ierr, ninfiles, nwords
    INTEGER(kind=hid_t) :: file_id      ! File identifier
    INTEGER             :: fids(ndatakind)  ! file IDs
                                            ! Index 1 = i_vrad
                                            ! Index 2 = i_qualvrad
                                            ! Index 3 = i_dbzh
                                            ! Index 4 = i_qual_dbzh
                                            ! Index 5 = i_zdr
                                            ! Index 6 = i_kdp
                                            ! Index 7 = i_phidp
                                            ! Index 8 = i_rhv
                                            ! Index 9 = i_ldr
                                            ! Index 10 = i_cflags

    INTEGER             :: ival, year, fidindlist(ndatakind)
    INTEGER             :: mii, station_id, station_id_last, naz, naz_last, nra
    REAL(KIND=dp)       :: dval, magn
    REAL(KIND=dp)       :: alt_msl_true, lat, lon, az_inc, az_start, ra_inc, az_start_last
    CHARACTER(len=500)  :: cval
    CHARACTER(len=10), ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    CHARACTER(len=20)   :: quantity, quantity_last, dsetlist(ndatakind)
    CHARACTER(len=60)   :: cpath
    CHARACTER(len=7)    :: juldate   ! YYYYMMDD
    CHARACTER (LEN=14)  :: cdate, cdate_last
    CHARACTER (len=5)   :: statid_char
    LOGICAL             :: station_found

    REAL(KIND=dp), ALLOCATABLE :: manel(:), v_nyq(:), high_prf(:), dualprf_ratio(:)

    err = 0

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' for '//TRIM(infilename)//' on proc ', my_radar_id

    ! Re-construct the filename series of the PPI-files from the given input filename:
    !  Examples:
    !  - ras11-vol5minng01_sweeph5allm_any_10-2020020300205700-tur-10832-hd5|00-2057|01-2120|02-2143|03-2207|04-2230|05-2300|06-2322|07-2335|08-2348|09-2402
    !  - rasXX-pcpng01_sweeph5allm_any_00-2020020300053400-neu-10557-hd5

    k = INDEX(infilename,'_',back=.TRUE.)   ! position of the last '_' in the string
    IF (infilename(7:13) == 'vol5min') THEN
      ! .. This is a volume scan distributed over several files:
      READ (infilename(k+1:k+2), '(i2.2)') ninfiles
      kk = INDEX(infilename,'|') - 1   ! length of the actual filename (anything before the first '|'
      CALL split_string (TRIM(infilename), '|', LEN(words), words, nwords)
      IF (nwords /= ninfiles+1) THEN
        WRITE(*,'(a)') 'ERROR '//TRIM(yzroutine)//': Inconsistent infilename '//TRIM(infilename)
        IF (ALLOCATED(words)) DEALLOCATE(words)
        err = 50
        RETURN
      END IF
    ELSE IF (infilename(7:9) == 'pcp') THEN
      ! .. This is a precipitation scan in one file:
      ninfiles = 1
      kk       = LEN_TRIM(infilename)
    END IF
    ALLOCATE(infilenames(ninfiles))
    DO i=1, ninfiles
      infilenames(i) = infilename(1:kk)
      IF (infilename(7:13) == 'vol5min') THEN
        infilenames(i)(k+1:k+2)   = words(i+1)(1:2)
        infilenames(i)(k+14:k+17) = words(i+1)(4:7)
      END IF
    END DO
    IF (ALLOCATED(words)) DEALLOCATE(words)

    ! .. Now all filenames of the respective timestep / volume scan are known, so read their metadata:

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      err = 1
      RETURN
    ENDIF

    ! .. Initializations for the below loop over all files of the volume scan:
    ALLOCATE (manel(ninfiles));         manel(:)         = miss_value
    ALLOCATE (v_nyq(ninfiles));         v_nyq(:)         = miss_value
    ALLOCATE (high_prf(ninfiles));      high_prf(:)      = miss_value
    ALLOCATE (dualprf_ratio(ninfiles)); dualprf_ratio(:) = miss_value
    cdate(:) = ' '
    cdate_last = 'YYYYMMDDhhmmss'
    quantity_last(:) = 'X'
    station_id = 999999
    station_id_last = 999999
    mii = missval_int
    naz = missval_int
    nra = missval_int
    ra_inc = miss_value
    alt_msl_true = miss_value
    lat = miss_value
    lon = miss_value
    az_inc = miss_value
    az_start = miss_value
    az_start_last = miss_value
    ! Initialize the type holding the radar meta data:
    rsmo = rsm_multitime2onetime ( get_meta_proto   ( icountry=i_dwd ) )

    DO ii=1, ninfiles

      ! Open the file for the particular timestep
      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(ii)), H5F_ACC_RDONLY_F, file_id, ierr)
      IF (ierr /= 0) THEN
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '// &
                        TRIM(ydirradarin)//TRIM(infilenames(ii))
        err = ierr
      ENDIF


!!$ File content as seen by "h5dump -n 1 ras11-vol5minng01_sweeph5allm_any_09-2020020301040300-hnr-10339-hd5":
!!$HDF5 "ras11-vol5minng01_sweeph5allm_any_09-2020020301040300-hnr-10339-hd5" {
!!$FILE_CONTENTS {
!!$ group      /
!!$ attribute  /Conventions
!!$ group      /dataset1
!!$ group      /dataset1/data1
!!$ dataset    /dataset1/data1/data
!!$ group      /dataset1/data1/what
!!$ attribute  /dataset1/data1/what/gain
!!$ attribute  /dataset1/data1/what/nodata
!!$ attribute  /dataset1/data1/what/offset
!!$ attribute  /dataset1/data1/what/quantity
!!$ attribute  /dataset1/data1/what/undetect
!!$ group      /dataset1/data10
!!$ dataset    /dataset1/data10/data
!!$ group      /dataset1/data10/what
!!$ attribute  /dataset1/data10/what/gain
!!$ attribute  /dataset1/data10/what/nodata
!!$ attribute  /dataset1/data10/what/offset
!!$ attribute  /dataset1/data10/what/quantity
!!$ attribute  /dataset1/data10/what/undetect
!!$ ...
!!$ group      /dataset1/how
!!$ attribute  /dataset1/how/NI
!!$ attribute  /dataset1/how/highprf
!!$ attribute  /dataset1/how/lowprf
!!$ attribute  /dataset1/how/polmode
!!$ attribute  /dataset1/how/pulsewidth
!!$ attribute  /dataset1/how/rpm
!!$ attribute  /dataset1/how/scan_index
!!$ attribute  /dataset1/how/startazA
!!$ attribute  /dataset1/how/startazT
!!$ attribute  /dataset1/how/startelA
!!$ attribute  /dataset1/how/stopazA
!!$ attribute  /dataset1/how/stopazT
!!$ attribute  /dataset1/how/stopelA
!!$ ...
!!$ group      /dataset1/what
!!$ attribute  /dataset1/what/enddate
!!$ attribute  /dataset1/what/endtime
!!$ attribute  /dataset1/what/product
!!$ attribute  /dataset1/what/startdate
!!$ attribute  /dataset1/what/starttime
!!$ group      /dataset1/where
!!$ attribute  /dataset1/where/a1gate
!!$ attribute  /dataset1/where/elangle
!!$ attribute  /dataset1/where/nbins
!!$ attribute  /dataset1/where/nrays
!!$ attribute  /dataset1/where/rscale
!!$ attribute  /dataset1/where/rstart
!!$ attribute  /dataset1/where/startaz
!!$ attribute  /dataset1/where/stopaz
!!$ ...
!!$ group      /what
!!$ attribute  /what/date
!!$ attribute  /what/object
!!$ attribute  /what/source
!!$ attribute  /what/time
!!$ attribute  /what/version
!!$ attribute  /what/version-sub
!!$ group      /where
!!$ attribute  /where/height
!!$ attribute  /where/lat
!!$ attribute  /where/lon
!!$ }

      ! .. read some attributes from different groups/datasets. The choice of the
      !    optional arguments has to match the data type of the attribute in the
      !    file:  dval = double/real/float   ,  ival = integer  ,  cval = character string
      !    A type mismatch might lead to a runtime crash!

      CALL read_attr_odim (file_id, '/where', 'height', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-alt_msl_true) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': alt_msl_true differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', alt_msl_true
        err = 6
      END IF
      alt_msl_true = dval

      CALL read_attr_odim (file_id, '/where', 'lat', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-lat) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': lat differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', lat
        err = 7
      END IF
      lat         = dval

      CALL read_attr_odim (file_id, '/where', 'lon', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-lon) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': lon differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', lon
        err = 8
      END IF
      lon         = dval

      CALL read_attr_odim (file_id, '/dataset1/where', 'nbins', ierr, ival=ival)
      nra = MAX(ival, nra)

      ! The nrays in the file might not always be the nominal number, which we
      !  want to use from the file. For the actual reading of the file, we
      !  need the actual number, which might be larger.
      CALL read_attr_odim (file_id, '/dataset1/where', 'nrays', ierr, ival=ival)
      IF (ii == 1) THEN
        naz_last = ival
      ELSE IF (ii > 1 .AND. ival /= naz_last) THEN
        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
             ': naz differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', ival, ' /= ', naz_last
        err = 9
        naz_last = ival
      END IF
      naz         = MAX(ival, naz)

      CALL read_attr_odim (file_id, '/dataset1/where', 'elangle', ierr, dval=manel(ii))

      ! Check/find out which variables are in the file:
      IF (fid > 0) THEN
        
 !!$ FIXME: ras07-files: here we use urhohv and uphidp, but these are problematic according to M. Werner. Check! Make namelist parameters instead of hardcoding names!
       ! This is a sweeph5onem file with one moment only, so we can check the consistency
        !  of the quantity across infilenames:
        quantity(:) = ' '
        CALL read_attr_odim (file_id, '/dataset1/data1/what', 'quantity', ierr, cval=quantity)
        IF (ii > 1 .AND. TRIM(quantity) /= TRIM(quantity_last)) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity differs from previous PPI file in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(quantity_last)
          err = 10
        END IF
        IF (fid == i_vrad .AND. TRIM(quantity) /= rsmo%obs_hdf5_varname_vrad(1:4) .AND. &
                                TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_vrad)) THEN
          ! VRAD and VRADH should be allowed, but not VRADV
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the filename in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_vrad)
          err = 11
        END IF
        IF (fid == i_dbzh .AND. TRIM(quantity) == rsmo%obs_hdf5_varname_dbzh(1:3) .AND. &
                                TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
          ! DBZ and DBZH should be allowed, but not DBZV
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the filename in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_dbzh)
          err = 12
        END IF
        !IF (fid == i_zdr .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_zdr)) THEN
        IF (fid == i_zdr .AND. INDEX( tolower(TRIM(quantity)), tolower(TRIM(rsmo%obs_hdf5_varname_zdr)) ) == 0 ) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the filename in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_zdr)
          err = 13
        END IF
        !IF (fid == i_kdp .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_kdp)) THEN
        IF (fid == i_kdp .AND. INDEX( tolower(TRIM(quantity)), tolower(TRIM(rsmo%obs_hdf5_varname_kdp)) ) == 0 ) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the filename in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_kdp)
          err = 14
        END IF
        !IF (fid == i_phidp .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_phidp)) THEN
        IF (fid == i_phidp .AND. INDEX( tolower(TRIM(quantity)), tolower(TRIM(rsmo%obs_hdf5_varname_phidp)) ) == 0 ) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the filename in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_phidp)
          err = 14
        END IF
        !IF (fid == i_rhv .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_rhv)) THEN
        IF (fid == i_rhv .AND. INDEX( tolower(TRIM(quantity)), tolower(TRIM(rsmo%obs_hdf5_varname_rhv)) ) == 0 ) THEN
          WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the filename in ' // &
               TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_rhv)
          err = 15
        END IF
        quantity_last = quantity
        
      ELSE
        
!!$ FIXME: ras11-files: here we use rhohv and uphidp. Check! Make namelist parameters instead of hardcoding names!
        ! This is a sweeph5allm file with many moments (all in /dataset1 for DWD), so we have to
        !  see which quantities are in the file which match the implemented fid identifiers:
        dsetlist(:)(:) = ' '
        fids(:)        = -1
        probe_loop: DO kk=1, 100
          cpath(:) = ' '
          WRITE (cpath, '("/dataset1/data",i0,"/what")') kk
          quantity(:) = ' '
          CALL read_attr_odim (file_id, TRIM(cpath), 'quantity', ierr, cval=quantity, silent=.TRUE.)
          IF (ierr == 0) THEN
            IF ( TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_vrad) .OR. &
                 TRIM(quantity) == rsmo%obs_hdf5_varname_vrad(1:4)) THEN
              ! This assumes that VRAD and VRADH are not in a file at the same time
              dsetlist(i_vrad)(:) = ' '
              WRITE (dsetlist(i_vrad), '("data",i0)') kk
              fids(i_vrad) = i_vrad
            ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_dbzh) .OR. &
                     TRIM(quantity) == rsmo%obs_hdf5_varname_dbzh(1:3)) THEN
              ! This assumes that DBZ and DBZH are not in a file at the same time
              dsetlist(i_dbzh)(:) = ' '
              WRITE (dsetlist(i_dbzh), '("data",i0)') kk
              fids(i_dbzh) = i_dbzh
            ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_zdr)) THEN
              dsetlist(i_zdr)(:) = ' '
              WRITE (dsetlist(i_zdr), '("data",i0)') kk
              fids(i_zdr) = i_zdr
            ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_kdp)) THEN
              dsetlist(i_kdp)(:) = ' '
              WRITE (dsetlist(i_kdp), '("data",i0)') kk
              fids(i_kdp) = i_kdp
            ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_phidp)) THEN
              dsetlist(i_phidp)(:) = ' '
              WRITE (dsetlist(i_phidp), '("data",i0)') kk
              fids(i_phidp) = i_phidp
            ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_rhv)) THEN
              dsetlist(i_rhv)(:) = ' '
              WRITE (dsetlist(i_rhv), '("data",i0)') kk
              fids(i_rhv) = i_rhv
            END IF
          ELSE
            ! There are no more datasets in the file, so exit the loop:
            EXIT probe_loop
          END IF
        END DO probe_loop
      END IF


      ! Find out the Nyquist velocity for radial wind data sets:
      !  In newer files there is an attribute '/dataset1/how/NI' (Nyquist Interval)
      !  If that is not found, we resort to the negative offset of the VRADH dataset.
      !  Here, the Nyquist velocity is (nearly) the negative offset, but we have to do some
      !  correction and rounding here to compensate for differences between
      !  the offset and the true Nyquist interval. Remember that the coded values
      !  of vr classes start at v_lower = offset + 1*gain and the last value is
      !  v_upper = offset + N*gain, so that v_nyq = offset + 0.5*gain is the best estimate.

      ! First try: Nyquist-velocity from /dataset1/how:
      CALL read_attr_odim (file_id, '/dataset1/how', 'NI', ierr, dval=dval, silent=.TRUE.)
      IF (ierr == 0) THEN
        v_nyq(ii) = dval
      ELSE
        ! Second try: use the data set for VRADH to look at offset:
        IF (fid < 1) THEN
          ! This is a sweeph5allm file with many moments. Use the data set for VRADH to look at -offset:
          IF (fids(i_vrad) > 0) THEN
            cpath(:) = ' '
            cpath = '/dataset1/'//TRIM(dsetlist(i_vrad))//'/what'
            CALL read_attr_odim (file_id, TRIM(cpath), 'offset', ierr, dval=dval)
            v_nyq(ii) = -dval
            CALL read_attr_odim (file_id, TRIM(cpath), 'gain', ierr, dval=dval)
            v_nyq(ii) = v_nyq(ii) - 0.5_dp*dval ! "true" Nyquist veloc. from DWD hdf5 files
!            ! round to 2.5 significant digits:
!            magn = 10.0_dp**FLOOR(LOG10(v_nyq(ii)))
!            v_nyq(ii) = NINT(v_nyq(ii)/magn*20.0_dp)*0.05_dp * magn
! instead, round differently for precip and volscans:
!!$            IF (infilename(7:9) == 'pcp') THEN
!!$              ! round to 0.25 m/s:
!!$              v_nyq(ii) = NINT(v_nyq(ii)/0.25_dp)*0.25_dp
!!$            ELSE
!!$              ! round to 0.5 m/s:
!!$              v_nyq(ii) = NINT(v_nyq(ii)/0.5_dp)*0.5_dp
!!$            END IF
          END IF
        ELSE
          ! This is a sweeph5onem file with one moment only:
          IF (fid == i_vrad) THEN
            CALL read_attr_odim (file_id, '/dataset1/data1/what', 'offset', ierr, dval=dval)
            v_nyq(ii) = -dval  ! Nyquist-velocity is the negative offset
            CALL read_attr_odim (file_id, '/dataset1/data1/what', 'gain', ierr, dval=dval)
            v_nyq(ii) = v_nyq(ii) - 0.5_dp*dval ! "true" Nyquist veloc. from DWD hdf5 files
!            ! round to 2.5 significant digits:
!            magn = 10.0_dp**FLOOR(LOG10(v_nyq(ii)))
!            v_nyq(ii) = NINT(v_nyq(ii)/magn*20.0_dp)*0.05_dp * magn
! instead, round differently for precip and volscans:
!!$            IF (infilename(7:9) == 'pcp') THEN
!!$              ! round to 0.25 m/s:
!!$              v_nyq(ii) = NINT(v_nyq(ii)/0.25_dp)*0.25_dp
!!$            ELSE
!!$              ! round to 0.5 m/s:
!!$              v_nyq(ii) = NINT(v_nyq(ii)/0.5_dp)*0.5_dp
!!$            END IF
          END IF
        END IF
      END IF

      CALL read_attr_odim (file_id, '/dataset1/how', 'highprf', ierr, dval=dval)
      high_prf(ii) = dval
      CALL read_attr_odim (file_id, '/dataset1/how', 'lowprf', ierr, dval=dval)
      dualprf_ratio(ii) = high_prf(ii) / dval

      CALL read_attr_odim (file_id, '/dataset1/where', 'rscale', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-ra_inc) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': ra_inc differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', ra_inc
        err = 13
      END IF
      ra_inc         = dval

!      CALL read_attr_odim (file_id, '/dataset1/where', 'rstart', ierr, dval=dval)  ! rstart = ra_inc!

!!$ Because naz from the file is not reliable, we also take az_inc from the background list!
!!$      IF (naz > 0) THEN
!!$        az_inc = 360.0_dp / naz
!!$      ELSE
!!$        az_inc = 1.0_dp
!!$      END IF
!!$
!!$ Also, startaz is not in every DWD hdf5-file. Therefore we take it also from the background list!
!!$      CALL read_attr_odim (file_id, '/dataset1/where', 'startaz', ierr, dval=dval)
!!$      IF (ii > 1 .AND. ABS(dval-az_start_last) > 1E-6_dp) THEN
!!$        WRITE (*,'(a,f10.4,a,f10.4)') 'WARNING '//TRIM(yzroutine)// &
!!$             ': az_start differs from previous PPI file in ' // &
!!$             TRIM(infilenames(ii))//': ', dval, ' /= ', az_start_last
!!$      END IF
!!$      az_start_last = dval
!!$      az_start = dval + 0.5_dp*rsmo%az_inc   ! Use the nominal value of az_inc to compute the center of the azimut bin

      CALL read_attr_odim (file_id, '/what', 'source', ierr, cval=cval)
      CALL get_stationID_from_HDF5_source (cval, station_id, mii, statid_char, ierr)
      IF (ierr /= 0) THEN
        err = 14
      END IF
      IF (ii > 1 .AND. station_id /= station_id_last) THEN
        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
             ': station_id differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', station_id, ' /= ', station_id_last
        err = 15
      END IF
      station_id_last = station_id

      ! station_name: overtake later from background list

      ! date and time (start time of the volume scan in this case):
      cdate(:) = ' '
      CALL read_attr_odim (file_id, '/what', 'date', ierr, cval=cval)
      IF ( ierr /= 0) THEN
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
             ': date could not be read from file ' // &
             TRIM(infilename)
        err = 16
      END IF
      cdate = TRIM(cval)

      CALL read_attr_odim (file_id, '/what', 'time', ierr, cval=cval)
      IF ( ierr /= 0) THEN
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
             ': time could not be read from file ' // &
             TRIM(infilenames(ii))
        err = 17
      END IF
      cdate = TRIM(cdate)//TRIM(cval)
      ! The nominal start time of the volume scan is the next lower full 5 min interval, so we round down to it:
      IF (rsmo%dt_obs(3) > 0.0_dp) THEN
        ! assuming new notation, where the time increment is in dt_obs(3):
        cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(3)/60))
      ELSE IF (rsmo%dt_obs(1) > 0.0_dp) THEN
        ! assuming old notation, where the time increment is in dt_obs(1) and (2), (3) are miss_value:
        cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(1)/60))        
      ELSE
        ! This should not happen, except that the user give explicitly a wrong dt_obs-triplet in the namelist:
        WRITE (*,'(a,3(1x,f0.1))') 'WARNING '//TRIM(yzroutine)// &
             ': specification of dt_obs(1:3) is wrong for station in file ' // &
             TRIM(infilenames(ii))//': dt_obs =', rsmo%dt_obs
        err = 18
      END IF
      
      ! If in future we want to use the nominal endtime, round cdate up to the next full 5 min to denote the nominal end time of the scan
      ! NOTE: if cdate is exactly a full 5 min time, it will be rounded up to the next upper 5 min interval!
!!$      cdate = round_datestr_min(cdate, NINT(rsmo%dt_obs(1)/60))  ! care for dt_obs(1) and dt_obs(3)!!!                                   
      IF (ii > 1 .AND. cdate /= cdate_last) THEN
        WRITE (*,'(a,a,a,a)') 'WARNING '//TRIM(yzroutine)// &
             ': date/time differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', cdate, ' /= ', cdate_last
        err = 19
      END IF
      cdate_last = cdate

      CALL h5fclose_f(file_id, ierr)

    END DO

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)


    ! .. Search the actual station in the background list and use it's meta data
    !     as a background list:
    station_found = .FALSE.
    DO i = 1, SIZE(rsmm_bg)
      IF ( station_id == rsmm_bg(i)%station_id ) THEN
        ! .. Copy the base values into rs_meta_ncdf(ista).
        rsmo = rsm_multitime2onetime ( rsmm_bg(i) )
        station_found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT.station_found) THEN
      err = 3
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // 'WARNING '//TRIM(yzroutine)//': radar station ', station_id, &
           ' not found in background list rsmm_bg!'
      DEALLOCATE(manel, v_nyq, high_prf, dualprf_ratio)
      IF (ALLOCATED(infilenames)) DEALLOCATE(infilenames)
      RETURN
    END IF

    ! .. Now partially override with the meta data from the file:
    !     (rsmo%station_name is taken from rsmm_bg)
    rsmo%station_id       = station_id
    rsmo%nel              = ninfiles
    rsmo%alt_msl_true     = alt_msl_true
    rsmo%ra_inc           = ra_inc
    rsmo%ra_inc_obs       = ra_inc
!!$    rsmo%az_start         = az_start ! not available in every file, use background val
!!$    rsmo%az_inc           = az_inc   ! not reliable from the file, use background val
!!$    rsmo%naz              = naz      ! not reliable from the file, use background val
    rsmo%nra              = nra
    rsmo%nra_obs          = nra
    rsmo%lon              = lon
    rsmo%lat              = lat

    IF (fid < 1) THEN
      ! all variables are in one file:
      rsmo%nrep_ncdf(:)     = ninfiles
      rsmo%naz_ncdf(:,:)    = naz  ! max. value of all infiles of this volume scan, can be > 360!
      DO i=1, ndatakind
        IF (fids(i) > 0) THEN
          rsmo%obsfile(:,i) = TRIM(infilename)//'/'//TRIM(dsetlist(i))
        END IF
      END DO
    ELSE
      ! each variable is in a different file:
      rsmo%nrep_ncdf(fid)   = ninfiles
      rsmo%naz_ncdf(:,fid)  = naz  ! max. value of all infiles of this volume scan, can be > 360!
      ! .. Store also the infilename in the rsm structure:
      rsmo%obsfile(:,fid)   = TRIM(infilename)
    END IF

    ! Extract obstime from vector cdate; the only possible error
    !  status of this subroutine (cdate not valid)  equals err_ista(ista) = 5
    !  from above, which has already been checked:
    IF (fid < 1) THEN
      CALL put_obstimes_to_rsm_onetime ( rsmo, TRIM(infilename), fids, ydate_ini_mod, cdate, ierr )
    ELSE
      CALL put_obstimes_to_rsm_onetime ( rsmo, TRIM(infilename), (/ fid /), ydate_ini_mod, cdate, ierr )
    END IF
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq, high_prf, dualprf_ratio)
      IF (ALLOCATED(infilenames)) DEALLOCATE(infilenames)
      RETURN
    END IF

    ! .. recognizes also DWD precip scans:
    CALL get_scanstrategy_ID ( rsmo, TRIM(infilename), manel, 0.1_dp, err )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq, high_prf, dualprf_ratio)
      RETURN
    ELSE
      rsmo%ext_nyq(:,:)     = unused_value
      rsmo%high_nyq(:)      = unused_value
      rsmo%prf(:)           = unused_value
      rsmo%dualprf_ratio(:) = unused_value
      IF (fid == i_vrad .OR. (fid < 1 .AND. ANY(v_nyq(:) >= 0.0))) THEN
        DO ii = 1, ninfiles
          DO i = 1, rsmo%nel
            IF (ABS(manel(ii)-rsmo%el_arr(i)) <= 0.1) THEN
              rsmo%ext_nyq(i,:)     = v_nyq(ii)
              rsmo%high_nyq(i)      = v_nyq(ii)
              rsmo%prf(i)           = high_prf(ii)
              rsmo%dualprf_ratio(i) = dualprf_ratio(ii)
            END IF
          END DO
        END DO
      END IF
    END IF

    ! .. Clean up:
    DEALLOCATE(manel, v_nyq, high_prf, dualprf_ratio)
    IF (ALLOCATED(infilenames)) DEALLOCATE(infilenames)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_metadata_from_h5_dwd

  !==============================================================================

  SUBROUTINE get_metadata_from_h5_mch ( infilename, fid, rsmo, rsmm_bg, err )

    !==============================================================================
    !
    ! .. HDF5-files from MCH: one file per station, per variable, per obs_time
    !    and per PPI-elevation (20 elevations)
    !
    !    This routine reads all the 20 files for one volume scan of one variable.
    !    The basename of the files for the volume scan is expected in the string
    !    "infilename".
    !
    !==============================================================================

    TYPE(radar_meta_type_onetime), INTENT(inout) :: rsmo
    TYPE(radar_meta_type),   INTENT(in)    :: rsmm_bg(:)  ! Background-list for cross-checking
    CHARACTER(LEN=*),        INTENT(in)    :: infilename
    INTEGER,                 INTENT(in)    :: fid
    INTEGER,                 INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_metadata_from_h5_mch'

    CHARACTER(len=cobsflen), ALLOCATABLE :: infilenames(:)
    INTEGER                              :: i, k, ii, kk, ierr, ninfiles, nwords
    INTEGER(kind=hid_t) :: file_id      ! File identifier

    INTEGER             :: ival, year
    INTEGER             :: mii, station_id, station_id_last, naz, naz_last, nra
    REAL(KIND=dp)       :: dval
    REAL(KIND=dp)       :: alt_msl_true, lat, lon, az_inc, az_start, ra_inc
    CHARACTER(len=500)  :: cval
    CHARACTER(len=10), ALLOCATABLE :: words(:)    ! for results of subroutine split_string()
    CHARACTER(len=20)   :: quantity, quantity_last
    CHARACTER(len=7)    :: juldate   ! YYYYMMDD
    CHARACTER (LEN=14)  :: cdate, cdate_last
    CHARACTER (len=5)   :: statid_char
    LOGICAL             :: station_found

    REAL(KIND=dp), ALLOCATABLE :: manel(:), v_nyq(:)

    err = 0

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' for '//TRIM(infilename)//' on proc ', my_radar_id

    ! Re-construct the filename series of the PPI-files from the given input filename,
    !  which contains the basic filename, the number of PPI sweeps and the list of elevation identifiers:
    READ (infilename(16:18), '(i3.3)') ninfiles
    kk = INDEX(infilename,'|') - 1   ! length of the actual filename (anything before the first '|'
    CALL split_string (TRIM(infilename), '|', LEN(words), words, nwords)
    IF (nwords /= ninfiles+1) THEN
      WRITE(*,'(a)') 'ERROR '//TRIM(yzroutine)//': Inconsistent infilename '//TRIM(infilename)
      IF (ALLOCATED(words)) DEALLOCATE(words)
      err = 50
      RETURN
    END IF
    ALLOCATE(infilenames(ninfiles))
    DO i=1, ninfiles
      infilenames(i) = infilename(1:kk)
      WRITE(infilenames(i)(16:18), '(a3)') words(i+1)(1:3)
    END DO


    ! .. Now all filenames of the respective timestep / volume scan are known, so read their metadata:

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      err = 1
      RETURN
    ENDIF

    ! .. Initializations for the below loop over all files of the volume scan:
    ALLOCATE (manel(ninfiles)); manel(:)      = miss_value
    ALLOCATE (v_nyq(ninfiles)); v_nyq(:)      = miss_value
    cdate(:) = ' '
    cdate_last = 'YYYYMMDDhhmmss'
    quantity_last(:) = 'X'
    station_id = 999999
    station_id_last = 999999
    mii = missval_int
    naz = missval_int
    nra = missval_int
    ra_inc = miss_value
    alt_msl_true = miss_value
    lat = miss_value
    lon = miss_value
    az_inc = miss_value
    az_start = miss_value
    ! Initialize the type holding the radar meta data:
    rsmo = rsm_multitime2onetime ( get_meta_proto   ( icountry=i_meteoswiss ) )

    DO ii=1, ninfiles

      ! Open the file for the particular timestep
      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(ii)), H5F_ACC_RDONLY_F, file_id, ierr)
      IF (ierr /= 0) THEN
        WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '// &
                       TRIM(ydirradarin)//TRIM(infilenames(ii))
        err = ierr
      ENDIF


!!$ File content as seen by "h5dump -n 1  PLA1520500007U.006.V.h5":
!!$HDF5 "PLA1520500007U.006.V.h5" {
!!$FILE_CONTENTS {
!!$ group      /
!!$ attribute  /Conventions
!!$ group      /dataset1
!!$ group      /dataset1/data1
!!$ dataset    /dataset1/data1/data
!!$ attribute  /dataset1/data1/data/CLASS
!!$ attribute  /dataset1/data1/data/IMAGE_VERSION
!!$ group      /dataset1/data1/what
!!$ attribute  /dataset1/data1/what/gain
!!$ attribute  /dataset1/data1/what/nodata
!!$ attribute  /dataset1/data1/what/offset
!!$ attribute  /dataset1/data1/what/quantity
!!$ attribute  /dataset1/data1/what/undetect
!!$ group      /dataset1/what
!!$ attribute  /dataset1/what/enddate
!!$ attribute  /dataset1/what/endtime
!!$ attribute  /dataset1/what/product
!!$ attribute  /dataset1/what/startdate
!!$ attribute  /dataset1/what/starttime
!!$ group      /dataset1/where
!!$ attribute  /dataset1/where/a1gate
!!$ attribute  /dataset1/where/elangle
!!$ attribute  /dataset1/where/nbins
!!$ attribute  /dataset1/where/nrays
!!$ attribute  /dataset1/where/rscale
!!$ attribute  /dataset1/where/rstart
!!$ group      /what
!!$ attribute  /what/date
!!$ attribute  /what/object
!!$ attribute  /what/source
!!$ attribute  /what/time
!!$ attribute  /what/version
!!$ group      /where
!!$ attribute  /where/height
!!$ attribute  /where/lat
!!$ attribute  /where/lon
!!$ }


      ! .. read some attributes from different groups/datasets. The choice of the
      !    optional arguments has to match the data type of the attribute in the
      !    file:  dval = double/real/float   ,  ival = integer  ,  cval = character string
      !    A type mismatch might lead to a runtime crash!

      CALL read_attr_odim (file_id, '/where', 'height', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-alt_msl_true) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': alt_msl_true differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', alt_msl_true
        err = 6
      END IF
      alt_msl_true = dval

      CALL read_attr_odim (file_id, '/where', 'lat', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-lat) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': lat differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', lat
        err = 7
      END IF
      lat         = dval

      CALL read_attr_odim (file_id, '/where', 'lon', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-lon) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': lon differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', lon
        err = 8
      END IF
      lon         = dval

      CALL read_attr_odim (file_id, '/dataset1/where', 'nbins', ierr, ival=ival)
      nra = MAX(INT(ival), nra)

      CALL read_attr_odim (file_id, '/dataset1/where', 'nrays', ierr, ival=ival)
      IF (ii == 1) THEN
        naz_last = ival
      ELSE IF (ii > 1 .AND. ival /= naz_last) THEN
        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
             ': naz differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', ival, ' /= ', naz_last
        err = 9
        naz_last = ival
      END IF
      naz              = ival

      CALL read_attr_odim (file_id, '/dataset1/where', 'elangle', ierr, dval=manel(ii))

      IF (fid == i_vrad) THEN
        CALL read_attr_odim (file_id, '/dataset1/data1/what', 'offset', ierr, dval=dval)
        v_nyq(ii) = -dval    ! For VRAD data, Nyquist veloc. is the negative offset
      ELSE
        v_nyq(ii) = miss_value
      END IF

      quantity(:) = ' '
      CALL read_attr_odim (file_id, '/dataset1/data1/what', 'quantity', ierr, cval=quantity)
      IF (ii > 1 .AND. TRIM(quantity) /= TRIM(quantity_last)) THEN
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
             ': quantity differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(quantity_last)
        err = 10
      END IF
      IF (fid == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_vrad)) THEN
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
             ': quantity in file differs from that of the filename in ' // &
             TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_vrad)
        err = 11
      END IF
      IF (fid == i_dbzh .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
        WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
             ': quantity in file differs from that of the filename in ' // &
             TRIM(infilenames(ii))//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_dbzh)
        err = 12
      END IF
      quantity_last = quantity

      CALL read_attr_odim (file_id, '/dataset1/where', 'rscale', ierr, dval=dval)
      IF (ii > 1 .AND. ABS(dval-ra_inc) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': ra_inc differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', dval, ' /= ', ra_inc
        err = 13
      END IF
      ra_inc         = dval

!      CALL read_attr_odim (file_id, '/dataset1/where', 'rstart', ierr, dval=dval)  ! rstart = ra_inc!

      IF (naz > 0) THEN
        az_inc = 360.0_dp / naz
      ELSE
        az_inc = 1.0_dp
      END IF

!!$ "a1gate" is always 0 here, but originally it should be the azimut index of the true start of the PPI scan
!!$ and it is not useful to determine the azimut start of the data in the file itself.
!!$      CALL read_attr_odim (file_id, '/dataset1/where', 'a1gate', ierr, ival=ival)
!!$      IF (ii > 1 .AND. ival /= iaz_start_last) THEN
!!$        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
!!$             ': az_start differs from previous PPI file in ' // &
!!$             TRIM(infilenames(ii))//': ', ival, ' /= ', iaz_start_last
!!$      END IF
!!$      iaz_start_last = ival
!!$      az_start = REAL(ival, KIND=dp) + 0.5_dp*az_inc

      ! Nominal first azimut is always half the azimut increment for Swiss radars:
      az_start = 0.5_dp*az_inc

      CALL read_attr_odim (file_id, '/what', 'source', ierr, cval=cval)
      CALL get_stationID_from_HDF5_source (cval, station_id, mii, statid_char, ierr)
      IF (ierr /= 0) THEN
        err = 14
      END IF
      IF (ii > 1 .AND. station_id /= station_id_last) THEN
        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
             ': station_id differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', station_id, ' /= ', station_id_last
        err = 15
      END IF
      station_id_last = station_id

      ! station_name: overtake later from background list

      ! The nominal date (end of volume scan) is only available from the filename, in form of the
      ! Julian day of year (YYJUL), where YY are the last 2 digits of the year and JUL is the day number from
      !  1 to 366.
      ! The given dates in the attributes are the true start- and endtimes of the PPI scans!
      READ (infilenames(ii)(4:5), '(i2.2)') year
      IF (year < 90) THEN
        juldate = '20' // infilenames(ii)(4:8)
      ELSE
        juldate = '19' // infilenames(ii)(4:8)
      END IF
      cdate = jul2yyyymmdd (juldate) // infilenames(ii)(9:12) // '00'
      ! Rewind the nominal end time to the nominal start time of the volume scan:
      cdate = new_datetime(cdate, -rsmo%dt_obs(3))
      IF (ii > 1 .AND. cdate /= cdate_last) THEN
        WRITE (*,'(a,a,a,a)') 'WARNING '//TRIM(yzroutine)// &
             ': date/time differs from previous PPI file in ' // &
             TRIM(infilenames(ii))//': ', cdate, ' /= ', cdate_last
        err = 16
      END IF
      cdate_last = cdate

      CALL h5fclose_f(file_id, ierr)

    END DO

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)


    ! .. Search the actual station in the background list and use it's meta data
    !     as a background list:
    station_found = .FALSE.
    DO i = 1, SIZE(rsmm_bg)
      IF ( station_id == rsmm_bg(i)%station_id ) THEN
        ! .. Copy the base values into rs_meta_ncdf(ista).
        rsmo = rsm_multitime2onetime ( rsmm_bg(i) )
        station_found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT.station_found) THEN
      err = 3
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // 'WARNING '//TRIM(yzroutine)//': radar station ', station_id, &
           ' not found in background list rsmm_bg!'
      DEALLOCATE(manel, v_nyq)
      RETURN
    END IF

    ! .. Now partially override with the meta data from the file:
    !     (rsmo%station_name is taken from rsmm_bg)
    rsmo%station_id       = station_id
    rsmo%nrep_ncdf(fid)   = ninfiles
    rsmo%nel              = ninfiles
    rsmo%alt_msl_true     = alt_msl_true
    rsmo%ra_inc           = ra_inc
    rsmo%ra_inc_obs       = ra_inc
    rsmo%az_start         = az_start
    rsmo%az_inc           = az_inc
    rsmo%naz              = naz
    rsmo%naz_ncdf(:,fid)  = naz
    rsmo%nra              = nra
    rsmo%nra_obs          = nra
    rsmo%lon              = lon
    rsmo%lat              = lat

    ! Extract obstime from vector cdate; the only possible error
    !  status of this subroutine (cdate not valid)  equals err_ista(ista) = 5
    !  from above, which has already been checked:
    CALL put_obstimes_to_rsm_onetime ( rsmo, TRIM(infilename), (/ fid /), ydate_ini_mod, cdate, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq)
      RETURN
    END IF

    ! .. recognizes also DWD precip scans:
    CALL get_scanstrategy_ID ( rsmo, TRIM(infilename), manel, 0.1_dp, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq)
      RETURN
    ELSE
      rsmo%ext_nyq(:,:) = unused_value
      rsmo%high_nyq(:)  = unused_value
      rsmo%prf(:)       = unused_value
      IF (fid == i_vrad) THEN
        DO ii = 1, ninfiles
          DO i = 1, rsmo%nel
            IF (ABS(manel(ii)-rsmo%el_arr(i)) <= 0.1) THEN
              rsmo%ext_nyq(i,:) = v_nyq(ii)
              rsmo%high_nyq(i)  = v_nyq(ii)
              rsmo%prf(i)       = NINT(4.0d0*v_nyq(ii)/rsmo%lambda) ! Compute from v_nyq and wavelength
            END IF
          END DO
        END DO
      END IF
    END IF

    ! .. Store also the base name of this file in the rsm structure:
    rsmo%obsfile(:,fid) = TRIM(infilename)

    DEALLOCATE(manel, v_nyq)

    ! .. Clean up:
    IF (ALLOCATED(infilenames)) DEALLOCATE(infilenames)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_metadata_from_h5_mch

  !==============================================================================

  SUBROUTINE get_metadata_from_h5_italy ( infilename, rsmo, rsmm_bg, err )

    !==============================================================================
    !
    ! .. HDF5-files from Italy: multi-moment multi-sweep files
    !
    !==============================================================================

    TYPE(radar_meta_type_onetime),   INTENT(inout) :: rsmo
    TYPE(radar_meta_type),   INTENT(in)    :: rsmm_bg(:)  ! Background-list for cross-checking
    CHARACTER(LEN=*),        INTENT(in)    :: infilename
    INTEGER,                 INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_metadata_from_h5_italy'

    INTEGER            :: i, j, k, ierr
    INTEGER(kind=hid_t):: file_id, root_id, grp_id, data_id
    INTEGER            :: storage_type, nrootmemb, max_corder, ndataset, obj_typ

    INTEGER            :: fids(ndatakind)  ! file IDs
                                           ! Index 1 = i_vrad
                                           ! Index 2 = i_qualvrad
                                           ! Index 3 = i_dbzh
                                           ! Index 4 = i_qualdbzh
    INTEGER            :: ival, year, fidindlist(ndatakind)
    INTEGER            :: mii, station_id, naz, naz_last, nra, nele
    REAL(KIND=dp)      :: dval
    REAL(KIND=dp)      :: alt_msl_true, lat, lon, az_inc, az_start, ra_inc
    CHARACTER(len=500) :: cval
    CHARACTER(len=20)  :: quantity, quantity_last
    CHARACTER(len=7)   :: juldate   ! YYYYMMDD
    CHARACTER(len=3)   :: cnumber
    CHARACTER (LEN=14) :: cdate
    CHARACTER (len=4)  :: nominaltime
    CHARACTER (len=5)  :: statid_char
    LOGICAL            :: station_found
    CHARACTER(len=20), ALLOCATABLE :: dsetname(:), vardataset(:)
    CHARACTER(len=45)  :: tmpgroup
    CHARACTER(len=LEN(rsmm_bg(1)%obsfile(1,1))) :: obsfile(ndatakind)

    REAL(KIND=dp), ALLOCATABLE :: manel(:), v_nyq(:)

    err = 0

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' for '//TRIM(infilename)//' on proc ', my_radar_id

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      err = 1
      RETURN
    ENDIF

    ! Open the file for the particular timestep
    CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilename), H5F_ACC_RDONLY_F, file_id, ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '// &
                     TRIM(ydirradarin)//TRIM(infilename)
      err = ierr
    ENDIF


!!$'/where/lat'  dval
!!$'/where/lon'  dval
!!$'/where/height' dval
!!$
!!$'/what/date' cval
!!$'/what/time' cval
!!$'/what/source' cval

!!$    ! Count the datasets
!!$    call h5gn_members_f(file_id, '/', nrootmemb, ierr)
!!$    ndataset=0
!!$    do i = 1, nrootmemb
!!$      call h5gget_obj_info_idx_f(file_id,'/',i-1,gr1_name,gr1_type,errorflag)
!!$      if ((gr1_name /= 'how') .and. (gr1_name /= 'what') .and. (gr1_name /= 'where')) then
!!$        ndataset=ndataset+1
!!$      endif
!!$    enddo

!!$'/dataset1/where/elangle' --> dval, elev(1)
!!$'/dataset2/where/elangle' --> dval, elev(2)
!!$...
!!$'/dataset1/where/nbins' --> ival, nra(1)
!!$'/dataset2/where/nbins' --> ival, nra(2)
!!$...
!!$'/dataset1/where/nrays' --> ival, naz_ncdf(1), az_inc = 360.0/ival
!!$'/dataset2/where/nrays' --> ival, naz_ncdf(2) (same as 01)
!!$...
!!$'/dataset1/where/rscale' --> dval, ra_inc
!!$'/dataset2/where/rscale' --> dval, ra_inc (same as 01)
!!$...
!!$
!!$! The moments are here:
!!$'/dataset1/data1/what/quantity' --> cval --> rs_meta%obsfile(:,1-4)
!!$'/dataset1/data1/data' --> ival --> iarr3D(:,:,1)
!!$'/dataset1/data2/what/quantity' --> cval --> rs_meta%obsfile(:,1-4)
!!$'/dataset1/data2/data' --> ival --> iarr3D(:,:,2)
!!$'/dataset2/data1/what/quantity' --> cval --> rs_meta%obsfile(:,1-4)
!!$'/dataset2/data1/data' --> ival --> iarr3D(:,:,1)
!!$'/dataset2/data2/what/quantity' --> cval --> rs_meta%obsfile(:,1-4)
!!$'/dataset2/data2/data' --> ival --> iarr3D(:,:,2)
!!$...

    ! .. identifier of the different data sets:
    fids(:) = -1

    CALL h5eset_auto_f(0, ierr)  ! Switch off HDF5 error messages temporarily,
                                 !  because the below loops will deliberately create
                                 !  errors when querying non-existent datasets

    ! .. Open the root group, below which the datasets reside:
    CALL h5gopen_f(file_id, '/', root_id, ierr)
    ! .. Get the number of members in the root group, which are candidates for datasets:
    CALL h5gget_info_f(root_id, storage_type, nrootmemb, max_corder, ierr)

    IF (ierr /= 0 .OR. nrootmemb <= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': No dataset found in HDF5-file '//TRIM(infilename)
      err = 17
    END IF

    ALLOCATE (dsetname(MAX(nrootmemb, 1)))
    dsetname(:)(:)= ' '

    ! .. search for groups named dataset1 ... dataset<nrootmemb> and store
    !    the names of the existing dataset-groups, which are fewer than
    !    nrootmemb:
    nele = 0
    DO i = 1, nrootmemb
      WRITE(cnumber, '(i3)') i
      dsetname(i) = 'dataset'//TRIM(ADJUSTL(cnumber))
      CALL h5oopen_f(root_id, TRIM(dsetname(i)), grp_id, ierr)
      IF (ierr == 0) THEN
        CALL h5iget_type_f(grp_id, obj_typ, ierr)
        IF (obj_typ == H5I_GROUP_F) THEN
          nele = nele + 1
          dsetname(nele) = dsetname(i)
        END IF
      END IF
      CALL h5oclose_f(grp_id, ierr)
    ENDDO
    ! .. Close the root group:
    CALL h5gclose_f(root_id, ierr)

    ! .. find out how many variables are stored in the datasets by
    !    inspecting only the first dataset1:
    CALL h5gopen_f(file_id, '/'//TRIM(dsetname(1)), grp_id, ierr)
    CALL h5gget_info_f(grp_id, storage_type, nrootmemb, max_corder, ierr)

    IF (ierr /= 0 .OR. nrootmemb <= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': No data found in HDF5-file '//TRIM(infilename)
      err = 18
    END IF

    ALLOCATE (vardataset(MAX(nrootmemb, 1)))
    vardataset(:)(:) = ' '

    ! .. Within dataset1, search for different "data" groups data1 ... data<nrootmemb> which
    !    contain different datasets (reflectivity, radial wind, ...) and their metadata:
    ndataset = 0
    DO i = 1, nrootmemb
      WRITE(cnumber, '(i3)') i
      vardataset(i) = 'data'//TRIM(ADJUSTL(cnumber))
      CALL h5oopen_f(grp_id, TRIM(vardataset(i)), data_id, ierr)
      IF (ierr == 0) THEN
        CALL h5iget_type_f(data_id, obj_typ, ierr)
        IF (obj_typ == H5I_GROUP_F) THEN
          ndataset = ndataset + 1
          vardataset(ndataset) = vardataset(i)
        END IF
      END IF
      CALL h5oclose_f(data_id, ierr)
    END DO
    ! .. Close the data group:
    CALL h5oclose_f(grp_id, ierr)

    CALL h5eset_auto_f(1, ierr)  ! Switch on again HDF5 error messages

    ! THOMAS modification, added the following 2 lines (otherwise rsmo%obs_hdf5_varname_dbzh
    ! is not defined):
    ! Initialize the type holding the radar meta data:
    rsmo = rsm_multitime2onetime ( get_meta_proto   ( icountry=i_arpasim ) )

    ! loop over all datasets (quantity) of the first dataset group to find out which ones are present:
    DO i=1, ndataset
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/'//TRIM(vardataset(i))//'/what', &
                           'quantity', ierr, cval=quantity)
      IF (ierr == 0) THEN
        IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_vrad)) THEN
          obsfile(i_vrad)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_vrad) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
          obsfile(i_dbzh)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_dbzh) = i
        ELSE
          IF (ldebug_radsim) WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // ': quantity=', TRIM(quantity), &
               ' not implemented! Record discarded! '//TRIM(infilename)
        END IF
      ELSE
        err = ierr
      END IF
    END DO

    ! .. Initializations for the parameters to be read from the hdf5 attributes:
    ALLOCATE (manel(nele)); manel(:)      = unused_value
    ALLOCATE (v_nyq(nele)); v_nyq(:)      = unused_value
    cdate(:) = ' '
    quantity_last(:) = 'X'
    station_id = 999999
    mii = missval_int
    naz = missval_int
    nra = missval_int
    ra_inc = miss_value
    alt_msl_true = miss_value
    lat = miss_value
    lon = miss_value
    az_inc = miss_value
    az_start = miss_value
    ! Initialize the type holding the radar meta data:
    rsmo = rsm_multitime2onetime ( get_meta_proto ( icountry=i_arpasim ) )

    ! .. First the "overall" attributes:
    !===================================

    ! .. Station ID and coordinates:
    CALL read_attr_odim (file_id, '/what', 'source', ierr, cval=cval)
    ! .. Tweak missing WMO number in the "source" attribute:
    IF ( ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': source could not be read from file ' // &
           TRIM(infilename)
      err = 1
    END IF
    CALL get_stationID_from_HDF5_source (cval, station_id, mii, statid_char, ierr)
    IF (ierr /= 0) THEN
      err = 1
    END IF

    ! .. date and time, in this case the nominal start time of the scan:
    CALL read_attr_odim (file_id, '/what', 'date', ierr, cval=cval)
    IF ( ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': date could not be read from file ' // &
           TRIM(infilename)
      err = 2
    END IF
    cdate = TRIM(cval)

    CALL read_attr_odim (file_id, '/what', 'time', ierr, cval=cval)
    IF ( ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': time could not be read from file ' // &
           TRIM(infilename)
      err = 3
    END IF
    cdate = TRIM(cdate)//TRIM(cval)
    ! Round cdate down to the next full 5 min interval:
    IF (rsmo%dt_obs(3) > 0.0_dp) THEN
      ! assuming new notation, where the time increment is in dt_obs(3):
      cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(3)/60))
    ELSE IF (rsmo%dt_obs(1) > 0.0_dp) THEN
      ! assuming old notation, where the time increment is in dt_obs(1) and (2), (3) are miss_value:
      cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(1)/60))
    ELSE
      ! This should not happen, except that the user give explicitly a wrong dt_obs-triplet in the namelist:
      WRITE (*,'(a,3(1x,f0.1))') 'WARNING '//TRIM(yzroutine)// &
           ': specification of dt_obs(1:3) is wrong for station in file ' // &
           TRIM(infilename)//': dt_obs =', rsmo%dt_obs
      err = 18
    END IF

    ! THOMAS modification (VIRGINIA)
    ! Check if filename corresponds to odimtime read in the file
    ! If not, the nominal time of the file is taken as reference
    READ(infilename(14:17),'(a4)') nominaltime
    IF (cdate(9:12) /= TRIM(nominaltime)) THEN
      WRITE(*,*) "TENGO IL VALORE NOMINALE DEL FILE: ",nominaltime," (valore precedente= ",cdate(9:12),")"
      cdate(9:14) = TRIM(nominaltime)//'00'
    ENDIF


    ! station_name: overtake later from background list

    CALL read_attr_odim (file_id, '/where', 'height', ierr, dval=dval)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': alt_msl_true could not be read from file ' // &
           TRIM(infilename)
      err = 4
    END IF
    alt_msl_true = dval

    CALL read_attr_odim (file_id, '/where', 'lat', ierr, dval=dval)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': lat could not be read from file ' // &
           TRIM(infilename)
      err = 5
    END IF
    lat         = dval

    CALL read_attr_odim (file_id, '/where', 'lon', ierr, dval=dval)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': lat could not be read from file ' // &
           TRIM(infilename)
      err = 6
    END IF
    lon         = dval

    ! .. Loop over all datasets to read and cross-check all attributes:
    DO j=1, nele

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'nbins', ierr, ival=ival)
      nra = MAX(INT(ival), nra)

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'nrays', ierr, ival=ival)
      IF (j == 1) THEN
        naz_last = ival
      ELSE IF (j > 1 .AND. ival /= naz_last) THEN
        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
             ': naz differs from previous dataset in file ' // &
             TRIM(infilename)//': ', ival, ' /= ', naz_last
        err = 9
        naz_last = ival
      END IF
      naz              = ival

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'rscale', ierr, dval=dval)
      IF (j > 1 .AND. ABS(dval-ra_inc) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': ra_inc differs from previous dataset in file ' // &
             TRIM(infilename)//': ', dval, ' /= ', ra_inc
        err = 10
      END IF
      ra_inc         = dval

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'elangle', ierr, dval=manel(j))

!      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'rstart', ierr, dval=dval)  ! rstart = ra_inc!

      IF (naz > 0) THEN
        az_inc = 360.0_dp / naz
        ! just in case there are few more naz in the file than nominal, round az_inc to the next 0.02 degrees:
        az_inc = REAL(NINT(az_inc/0.02_dp)*0.02_dp, kind=dp)
      ELSE
        az_inc = 1.0_dp
      END IF

      ! Nominal first azimut is chosen as half the azimut increment.
      ! This is always the best rounded value of the true azimut of the first ray in the data,
      !  because for Italy, the first true antenna azimut with a value > 0.0 is stored at the first
      !  index.
      az_start = 0.5_dp*az_inc

      DO i=1, ndataset

        tmpgroup(:) = ' '
        tmpgroup    = '/'//TRIM(dsetname(j))//'/'//TRIM(vardataset(i))

        ! .. Check that "quantity" are consistent among "/datasetXX/dataYY":
        quantity(:) = ' '
        CALL read_attr_odim (file_id, TRIM(tmpgroup)//'/what', 'quantity', ierr, cval=quantity)
        IF (fids(i_vrad) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_vrad)) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
               TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_vrad)
          err = 12
        END IF
        IF (fids(i_dbzh) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
               TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_dbzh)
          err = 13
        END IF

        ! .. For radial wind, store the Nyquist velocity:
        IF (TRIM(quantity) == 'VRAD') THEN
          CALL read_attr_odim (file_id, TRIM(tmpgroup)//'/what', 'offset', ierr, dval=dval)
          v_nyq(j) = -dval    ! For VRAD data, Nyquist veloc. is the negative offset
        END IF


      END DO
    END DO

    ! .. Close file:
    CALL h5fclose_f(file_id, ierr)

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. Clean up memory:
    DEALLOCATE(dsetname, vardataset)


    ! .. Search the actual station in the background list and use it's meta data
    !     as a background list:
    station_found = .FALSE.
    DO i = 1, SIZE(rsmm_bg)
      IF ( station_id == rsmm_bg(i)%station_id ) THEN
        ! .. Copy the base values into rs_meta_ncdf(ista):
        rsmo =  rsm_multitime2onetime ( rsmm_bg(i) )
        station_found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT.station_found) THEN
      err = 3
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // 'WARNING '//TRIM(yzroutine)//': radar station ', station_id, &
           ' not found in background list rsmm_bg!'
      DEALLOCATE(manel, v_nyq)
      RETURN
    END IF

    ! .. Now partially override with the meta data from the file:
    !     (rsm%station_name is taken from rsmm_bg)

    rsmo%station_id       = station_id
    rsmo%nel              = nele
    rsmo%alt_msl_true     = alt_msl_true
    rsmo%ra_inc           = ra_inc
    rsmo%ra_inc_obs       = ra_inc
    rsmo%az_start         = az_start
    rsmo%az_inc           = az_inc
    rsmo%naz              = naz
    rsmo%nra              = nra
    rsmo%nra_obs          = nra
    rsmo%lon              = lon
    rsmo%lat              = lat
    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        rsmo%obsfile(:,i)   = obsfile(i)
        rsmo%nrep_ncdf(i)   = nele
        rsmo%naz_ncdf(:,i)  = naz
      END IF
    END DO

    ! Extract blocks of obstimes from cdate; the only possible error
    !  status of this subroutine (cdate not valid)  equals err_ista(ista) = 5
    !  from above, which has already been checked:
    fidindlist = -1
    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        fidindlist(i) = i
      END IF
    END DO
    CALL put_obstimes_to_rsm_onetime ( rsmo, TRIM(infilename), fidindlist, ydate_ini_mod, cdate, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq)
      RETURN
    END IF

    ! .. recognizes also DWD precip scans:
    CALL get_scanstrategy_ID ( rsmo, TRIM(infilename), manel, 0.2_dp, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq)
      RETURN
    ELSE
      rsmo%ext_nyq(:,:) = unused_value
      rsmo%high_nyq(:)  = unused_value
      rsmo%prf(:)       = unused_value
      IF (fids(1) > 0) THEN
        ! .. Radial wind data are present, so store Nyquist velocity information:
        DO j = 1, nele
          DO i = 1, rsmo%nel
            IF (ABS(manel(j)-rsmo%el_arr(i)) <= 0.1) THEN
              rsmo%ext_nyq(i,:) = v_nyq(j)
              rsmo%high_nyq(i)  = v_nyq(j)
              rsmo%prf(i)       = NINT(4.0d0*v_nyq(j)/rsmo%lambda) ! Compute from v_nyq and wavelength
            END IF
          END DO
        END DO
      END IF
    END IF

    DEALLOCATE(manel, v_nyq)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_metadata_from_h5_italy

  !==============================================================================

  SUBROUTINE get_metadata_from_h5_kitcband ( infilename, rsmo, rsmm_bg, err )

    !==============================================================================
    !
    ! .. HDF5-files from KIT Karlsruhe C-band: multi-moment multi-sweep files
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
    !==============================================================================

    TYPE(radar_meta_type_onetime),   INTENT(inout) :: rsmo
    TYPE(radar_meta_type),   INTENT(in)    :: rsmm_bg(:)  ! Background-list for cross-checking
    CHARACTER(LEN=*),        INTENT(in)    :: infilename
    INTEGER,                 INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_metadata_from_h5_kitcband'

    INTEGER            :: i, j, k, ierr
    INTEGER(kind=hid_t):: file_id, root_id, grp_id, data_id
    INTEGER            :: storage_type, nrootmemb, max_corder, ndataset, obj_typ

    INTEGER            :: fids(ndatakind)  ! file IDs
                                           ! Index 1 = i_vrad
                                           ! Index 2 = i_qualvrad
                                           ! Index 3 = i_dbzh
                                           ! Index 4 = i_qualdbzh
    INTEGER            :: ival, year, fidindlist(ndatakind)
    INTEGER            :: mii, station_id, naz, naz_last, nra, nele
    REAL(KIND=dp)      :: dval
    REAL(KIND=dp)      :: alt_msl_true, lat, lon, az_inc, az_start, ra_inc
    CHARACTER(len=500) :: cval
    CHARACTER(len=20)  :: quantity, quantity_last
    CHARACTER(len=7)   :: juldate   ! YYYYMMDD
    CHARACTER(len=3)   :: cnumber
    CHARACTER (LEN=14) :: cdate
    CHARACTER (len=4)  :: nominaltime
    CHARACTER (len=5)  :: statid_char
    LOGICAL            :: station_found
    CHARACTER(len=20), ALLOCATABLE :: dsetname(:), vardataset(:)
    CHARACTER(len=45)  :: tmpgroup
    CHARACTER(len=LEN(rsmm_bg(1)%obsfile(1,1))) :: obsfile(ndatakind)

    REAL(KIND=dp), ALLOCATABLE :: manel(:), v_nyq(:)

    err = 0

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' for '//TRIM(infilename)//' on proc ', my_radar_id

    ! Initialize HDF5 FORTRAN interface:
    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      err = 1
      RETURN
    ENDIF

    ! Open the file for the particular timestep
    CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilename), H5F_ACC_RDONLY_F, file_id, ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '// &
                     TRIM(ydirradarin)//TRIM(infilename)
      err = ierr
    ENDIF

    ! .. Read the station_id from the file:
    station_id = 999999
    CALL read_attr_odim (file_id, '/how', 'site_name', ierr, cval_varlen=cval)
    cval = TRIM(ADJUSTL(cval))
    READ (cval, *) station_id
    
    ! .. Search the actual station in the background list and use it's meta data
    !     as a background list:
    station_found = .FALSE.
    DO i = 1, SIZE(rsmm_bg)
      IF ( station_id == rsmm_bg(i)%station_id ) THEN
        ! .. Copy the base values into rs_meta_ncdf(ista):
        rsmo =  rsm_multitime2onetime ( rsmm_bg(i) )
        station_found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT.station_found) THEN
      err = 3
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // 'WARNING '//TRIM(yzroutine)//': radar station ', station_id, &
           ' not found in background list rsmm_bg!'
      DEALLOCATE(manel, v_nyq)
      RETURN
    END IF

!!$'/where/lat'  dval
!!$'/where/lon'  dval
!!$'/where/height' dval
!!$
!!$'/what/date' cval_varlen  -> 2023-07-12T12:35:03.000Z (startdate of volscan)
!!$'/how/site_name' cval_varlen  -> substitute for source: "20001"

!!$ datasets: /scan0 ... /scan13
!!$           /scanX/how
    
!!$'/scan0/how/elevation' --> dval, elev(1)
!!$'/scan1/how/elevation' --> dval, elev(2)
!!$...
!!$'/scan0/how/bin_count' --> ival, nra(1)   == nbins
!!$'/scan1/how/bin_count' --> ival, nra(2)
!!$...
!!$attribute  /scan0/how/angle_step   == az_inc
!!$attribute  /scan0/how/azi_start
!!$attribute  /scan0/how/azi_stop
!!$attribute  /scan0/how/PRF
!!$attribute  /scan0/how/ray_count         == nrays
!!$attribute  /scan0/how/radar_wave_length
!!$attribute  /scan0/how/range             == max. range
!!$attribute  /scan0/how/range_samples
!!$attribute  /scan0/how/range_step        ra_inc = range_samples*range_step
!!$attribute  /scan0/how/range_start       = 0
!!$attribute  /scan0/how/nyquist_velocity      = v_nyq (this is the original Nyquist-velocity corresponding to the lower of the two PRFs) ; also set high_nqy with that value
!!$attribute  /scan0/how/unambiguous_velocity  = ext_nyq  (this is the relevant Dual PRF Nyquist-velocity)
!!$
!!$attribute  /scan0/moment_0/moment  = name of moment
!!$attribute  /scan0/moment_0/dyn_range_min    = offset
!!$attribute  /scan0/moment_0/dyn_range_max    --> gain = (dyn_range_max-dyn_range_min) / 255

!!$dataset    /scan0/moment_0   = Zh
!!$dataset    /scan0/moment_0   = ZDR
!!$dataset    /scan0/moment_5   = Vh
!!$dataset    /scan0/moment_7   = PHIDP
!!$dataset    /scan0/moment_8   = KDP
!!$dataset    /scan0/moment_9   = RHOHV

!!$! The moments are here: BEWARE: moment_X is not a group as in ODIM, but a dataset!!!
!!$'/scan0/moment_0/moment' --> cval_varlen --> rs_meta%obsfile(:,1-10)
!!$'/scan0/moment_0' --> ival --> iarr3D(:,:,1)
!!$'/scan0/moment_1/moment' --> cval_varlen --> rs_meta%obsfile(:,1-10)
!!$'/scan0/moment_1' --> ival --> iarr3D(:,:,2)
!!$'/scan1/moment_0/moment' --> cval_varlen --> rs_meta%obsfile(:,1-10)
!!$'/scan1/moment_0' --> ival --> iarr3D(:,:,1)
!!$'/scan1/moment_1/moment' --> cval_varlen --> rs_meta%obsfile(:,1-10)
!!$'/scan1/moment_1' --> ival --> iarr3D(:,:,2)
!!$...

!!$ varlist from hdf5-file: "'Zh:1 Zv:1 UZh:1 UZv:1 ZDR:1 Vh:1 Wh:1 PHIDP:1 KDP:1 RHOHV:1 SQIh:1 CCORh:1 SNRh:1'"

!!$ Zh data scaling: val = dyn_range_min + (dyn_range_max-dyn_range_min) / 255 * ival
!!$                  if val = -32.0: val = no_echo

!!$ Vh data scaling: if ival == 0: val = no_data / miss_value
!!$                  else        : val = dyn_range_min + (dyn_range_max-dyn_range_min) / 254 * (ival-1)

    
    ! .. identifier of the different data sets:
    fids(:) = -1

    CALL h5eset_auto_f(0, ierr)  ! Switch off HDF5 error messages temporarily,
                                 !  because the below loops will deliberately create
                                 !  errors when querying non-existent datasets

    ! .. Open the root group, below which the datasets reside:
    CALL h5gopen_f(file_id, '/', root_id, ierr)
    ! .. Get the number of members in the root group, which are candidates for datasets:
    CALL h5gget_info_f(root_id, storage_type, nrootmemb, max_corder, ierr)

    IF (ierr /= 0 .OR. nrootmemb <= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': No dataset found in HDF5-file '//TRIM(infilename)
      err = 17
    END IF

    ALLOCATE (dsetname(MAX(nrootmemb, 1)))
    dsetname(:)(:)= ' '

    ! .. search for groups named scan1 ... scan<nrootmemb> and store
    !    the names of the existing dataset-groups, which are fewer than
    !    nrootmemb:
    nele = 0
    DO i = 1, nrootmemb
      WRITE(cnumber, '(i3)') i-1
      dsetname(i) = 'scan'//TRIM(ADJUSTL(cnumber))
      CALL h5oopen_f(root_id, TRIM(dsetname(i)), grp_id, ierr)
      IF (ierr == 0) THEN
        CALL h5iget_type_f(grp_id, obj_typ, ierr)
        IF (obj_typ == H5I_GROUP_F) THEN
          nele = nele + 1
          dsetname(nele) = dsetname(i)
        END IF
      END IF
      CALL h5oclose_f(grp_id, ierr)
    ENDDO
    ! .. Close the root group:
    CALL h5gclose_f(root_id, ierr)


    ! .. find out how many variables are stored in the datasets by
    !    inspecting only the first dataset1:
    CALL h5gopen_f(file_id, '/'//TRIM(dsetname(1)), grp_id, ierr)
    CALL h5gget_info_f(grp_id, storage_type, nrootmemb, max_corder, ierr)

    IF (ierr /= 0 .OR. nrootmemb <= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': No data found in HDF5-file '//TRIM(infilename)
      err = 18
    END IF

    ALLOCATE (vardataset(MAX(nrootmemb, 1)))
    vardataset(:)(:) = ' '

    ! .. Within scan1, search for different "moment_" groups moment_1 ... moment_<nrootmemb> which
    !    contain different datasets (reflectivity, radial wind, ...) and their metadata:
    ndataset = 0
    DO i = 1, nrootmemb
      WRITE(cnumber, '(i3)') i-1
      vardataset(i) = 'moment_'//TRIM(ADJUSTL(cnumber))
      CALL h5oopen_f(grp_id, TRIM(vardataset(i)), data_id, ierr)
      IF (ierr == 0) THEN
        CALL h5iget_type_f(data_id, obj_typ, ierr)
        IF (obj_typ == H5I_DATASET_F) THEN
          ndataset = ndataset + 1
          vardataset(ndataset) = vardataset(i)
        END IF
      END IF
      CALL h5oclose_f(data_id, ierr)
    END DO
    ! .. Close the data group:
    CALL h5oclose_f(grp_id, ierr)

    CALL h5eset_auto_f(1, ierr)  ! Switch on again HDF5 error messages


    ! loop over all datasets (quantity) of the first dataset group to find out which ones are present:
    DO i=1, ndataset
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/'//TRIM(vardataset(i)), &
           'moment', ierr, cval_varlen=quantity)
      IF (ierr == 0) THEN
        IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_vrad)) THEN
          obsfile(i_vrad)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_vrad) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
          obsfile(i_dbzh)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_dbzh) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_zdr)) THEN
          obsfile(i_zdr)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_zdr) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_kdp)) THEN
          obsfile(i_kdp)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_kdp) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_phidp)) THEN
          obsfile(i_phidp)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_phidp) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_rhv)) THEN
          obsfile(i_rhv)   = TRIM(infilename)//'/'//TRIM(vardataset(i))
          fids(i_rhv) = i
        ELSE
          IF (ldebug_radsim) WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // ': moment=', TRIM(quantity), &
               ' not implemented! Record discarded! '//TRIM(infilename)
        END IF
      ELSE
        err = ierr
      END IF
    END DO

    ! .. Initializations for the parameters to be read from the hdf5 attributes:
    ALLOCATE (manel(nele)); manel(:)      = unused_value
    ALLOCATE (v_nyq(nele)); v_nyq(:)      = unused_value
    cdate(:) = ' '
    quantity_last(:) = 'X'
    mii = missval_int
    naz = missval_int
    nra = missval_int
    ra_inc = miss_value
    alt_msl_true = miss_value
    lat = miss_value
    lon = miss_value
    az_inc = miss_value
    az_start = miss_value

    ! .. First the "overall" attributes:
    !===================================

    ! .. date and time, in this case the nominal start time of the scan:
    CALL read_attr_odim (file_id, '/what', 'date', ierr, cval_varlen=cval)
    IF ( ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': date could not be read from file ' // &
           TRIM(infilename)
      err = 2
    END IF
    ! 2023-07-12T12:35:03.000Z
    cval = TRIM(ADJUSTL(cval))
    cdate = cval(1:4)//cval(6:7)//cval(9:10)//cval(12:13)//cval(15:16)//cval(18:19)

    ! Round cdate down to the next full 5 min interval:
    IF (rsmo%dt_obs(3) > 0.0_dp) THEN
      ! assuming new notation, where the time increment is in dt_obs(3):
      cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(3)/60))
    ELSE IF (rsmo%dt_obs(1) > 0.0_dp) THEN
      ! assuming old notation, where the time increment is in dt_obs(1) and (2), (3) are miss_value:
      cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(1)/60))
    ELSE
      ! This should not happen, except that the user give explicitly a wrong dt_obs-triplet in the namelist:
      WRITE (*,'(a,3(1x,f0.1))') 'WARNING '//TRIM(yzroutine)// &
           ': specification of dt_obs(1:3) is wrong for station in file ' // &
           TRIM(infilename)//': dt_obs =', rsmo%dt_obs
      err = 18
    END IF

    ! station_name: overtake later from background list

    CALL read_attr_odim (file_id, '/where', 'height', ierr, dval=dval)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': alt_msl_true could not be read from file ' // &
           TRIM(infilename)
      err = 4
    END IF
    alt_msl_true = dval

    CALL read_attr_odim (file_id, '/where', 'lat', ierr, dval=dval)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': lat could not be read from file ' // &
           TRIM(infilename)
      err = 5
    END IF
    lat         = dval

    CALL read_attr_odim (file_id, '/where', 'lon', ierr, dval=dval)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
           ': lat could not be read from file ' // &
           TRIM(infilename)
      err = 6
    END IF
    lon         = dval

    ! .. Loop over all datasets to read and cross-check all attributes:
    DO j=1, nele

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'bin_count', ierr, ival=ival)
      nra = MAX(INT(ival), nra)

!!$      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'ray_count', ierr, ival=ival)
!!$      IF (j == 1) THEN
!!$        naz_last = ival
!!$        naz      = MIN(INT(ival), rsmo%naz)
!!$      ELSE IF (j > 1 .AND. ival /= naz_last) THEN
!!$        WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
!!$             ': naz differs from previous dataset in file ' // &
!!$             TRIM(infilename)//': ', ival, ' /= ', naz_last
!!$        naz_last = ival
!!$        IF (ival >= rsmo%naz) THEN
!!$          naz = MIN(INT(ival), naz)
!!$        ELSE
!!$          err = 9
!!$        END IF
!!$      END IF

      ! simplified:
      naz = rsmo%naz
              
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'range_samples', ierr, ival=ival)
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'range_step', ierr, dval=dval)
      dval = ival * dval  ! ra_inc = range_samples (number of averaged pulses) * range_step (1 pulselength)
      IF (j > 1 .AND. ABS(dval-ra_inc) > 1E-6_dp) THEN
        WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
             ': ra_inc differs from previous dataset in file ' // &
             TRIM(infilename)//': ', dval, ' /= ', ra_inc
        err = 10
      END IF
      ra_inc         = dval

      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'elevation', ierr, dval=manel(j))

!      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'range_start', ierr, dval=dval)  ! rstart = ra_inc!

      IF (naz > 0) THEN
        az_inc = 360.0_dp / naz
        ! just in case there are few more naz in the file than nominal, round az_inc to the next 0.02 degrees:
        az_inc = REAL(NINT(az_inc/0.02_dp)*0.02_dp, kind=dp)
      ELSE
        az_inc = 1.0_dp
      END IF

      ! Nominal first azimut is chosen as half the azimut increment.
      ! This is a good esimated value of the true azimut of the first ray in the data
      az_start = 0.5_dp*az_inc

      DO i=1, ndataset

        tmpgroup(:) = ' '
        tmpgroup    = '/'//TRIM(dsetname(j))//'/'//TRIM(vardataset(i))

        ! .. Check that "moment" are consistent among "/scanXX/moment_YY":
        quantity(:) = ' '
        CALL read_attr_odim (file_id, TRIM(tmpgroup), 'moment', ierr, cval_varlen=quantity, silent=.TRUE.)
        IF (ierr == 0) THEN
          IF (fids(i_vrad) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_vrad)) THEN
            WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                 ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                 TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_vrad)
            err = 12
          END IF
          IF (fids(i_dbzh) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
            WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                 ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                 TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_dbzh)
            err = 13
          END IF
          IF (fids(i_zdr) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_zdr)) THEN
            WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                 ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                 TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_zdr)
            err = 14
          END IF
          IF (fids(i_kdp) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_kdp)) THEN
            WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                 ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                 TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_kdp)
            err = 15
          END IF
          IF (fids(i_phidp) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_phidp)) THEN
            WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                 ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                 TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_phidp)
            err = 16
          END IF
          IF (fids(i_rhv) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_rhv)) THEN
            WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                 ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                 TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_rhv)
            err = 17
          END IF
        ELSE
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': missing dataset '//TRIM(tmpgroup)//' in file '//TRIM(infilename)          
        END IF
      END DO

      ! .. For radial wind, store the Nyquist velocity:
      CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'unambiguous_velocity', ierr, dval=dval)
      v_nyq(j) = dval

    END DO

    ! .. Close file:
    CALL h5fclose_f(file_id, ierr)

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. Clean up memory:
    DEALLOCATE(dsetname, vardataset)



    ! .. Now partially override with the meta data from the file:
    !     (rsm%station_name is taken from rsmm_bg)

    rsmo%station_id       = station_id
    rsmo%nel              = nele
    rsmo%alt_msl_true     = alt_msl_true
    rsmo%ra_inc           = ra_inc
    rsmo%ra_inc_obs       = ra_inc
    rsmo%az_start         = az_start
    rsmo%az_inc           = az_inc
    rsmo%naz              = naz
    rsmo%nra              = nra
    rsmo%nra_obs          = nra
    rsmo%lon              = lon
    rsmo%lat              = lat
    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        rsmo%obsfile(:,i)   = obsfile(i)
        rsmo%nrep_ncdf(i)   = nele
        rsmo%naz_ncdf(:,i)  = naz
      END IF
    END DO

    ! Extract blocks of obstimes from cdate; the only possible error
    !  status of this subroutine (cdate not valid)  equals err_ista(ista) = 5
    !  from above, which has already been checked:
    fidindlist = -1
    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        fidindlist(i) = i
      END IF
    END DO
    CALL put_obstimes_to_rsm_onetime ( rsmo, TRIM(infilename), fidindlist, ydate_ini_mod, cdate, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq)
      RETURN
    END IF

    ! .. recognizes also DWD precip scans:
    CALL get_scanstrategy_ID ( rsmo, TRIM(infilename), manel, 0.1_dp, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      DEALLOCATE(manel, v_nyq)
      RETURN
    ELSE
      rsmo%ext_nyq(:,:) = unused_value
      rsmo%high_nyq(:)  = unused_value
      rsmo%prf(:)       = unused_value
      IF (fids(1) > 0) THEN
        ! .. Radial wind data are present, so store Nyquist velocity information:
        DO j = 1, nele
          DO i = 1, rsmo%nel
            IF (ABS(manel(j)-rsmo%el_arr(i)) <= 0.1) THEN
              rsmo%ext_nyq(i,:) = v_nyq(j)
              rsmo%high_nyq(i)  = v_nyq(j)
              rsmo%prf(i)       = NINT(4.0d0*v_nyq(j)/rsmo%lambda) ! Compute from v_nyq and wavelength
            END IF
          END DO
        END DO
      END IF
    END IF

    DEALLOCATE(manel, v_nyq)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_metadata_from_h5_kitcband

  !==============================================================================

  SUBROUTINE get_metadata_from_h5_opera ( infilename, rsmo, rsmm_bg, fileformat, err )

    !==============================================================================
    !
    ! .. HDF5-files of the type that is sent to OPERA by each country. The files
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
    !    T_PA<char1-scan-ident><char1-ele-ident><int2-station-discr>_<char4-center-ident>_<YYYYMMDDhhmmss>.<ext>
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
    ! .. This reader is also used for native hdf5-files from MeteoSwiss, which
    !     are single-moment / or multi-moment / single-sweep.
    !
    !==============================================================================

    CHARACTER(LEN=*),              INTENT(in)    :: infilename
    TYPE(radar_meta_type_onetime), INTENT(inout) :: rsmo
    TYPE(radar_meta_type),         INTENT(in)    :: rsmm_bg(:)  ! Background-list for cross-checking
    CHARACTER(LEN=LEN(rsmo%obsfile_format)), INTENT(out) :: fileformat
    INTEGER,                       INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'get_metadata_from_h5_opera'

    CHARACTER(len=cobsflen), ALLOCATABLE :: infilenames(:)
    INTEGER            :: i, ii, iii, j, k, ierr, ninfiles, tmpdim
    INTEGER(kind=hid_t):: file_id
    INTEGER            :: ndataset

    INTEGER            :: fids(ndatakind)  ! file IDs
                                           ! Index 1 = i_vrad
                                           ! Index 2 = i_qualvrad
                                           ! Index 3 = i_dbzh
                                           ! Index 4 = i_qualdbzh
    INTEGER            :: ival, year, fidindlist(ndatakind)
    INTEGER            :: mii, station_id,  station_id_last, naz, naz_last, nra, nele, nele_obs
    REAL(KIND=dp)      :: dval
    REAL(KIND=dp)      :: alt_msl_true, lat, lon, az_inc, az_start, ra_inc, az_start_last
    CHARACTER(len=500) :: cval, csource
    CHARACTER(len=20)  :: quantity, quantity_last
    CHARACTER(len=7)   :: juldate   ! YYYYMMDD
    CHARACTER (LEN=14) :: cdate, cdate_last
    CHARACTER(len=60)  :: cpath
    CHARACTER (len=5)  :: statid_char
    LOGICAL            :: station_found, ele_found
    LOGICAL            :: is_multimoment, is_multisweep, is_native_format

    CHARACTER(len=20), ALLOCATABLE :: dsetname(:), vardataset(:)
    CHARACTER(len=45)  :: tmpgroup
    CHARACTER(len=4)   :: provider

    REAL(KIND=dp), ALLOCATABLE :: manel(:), v_nyq(:), high_prf(:), dualprf_ratio(:)

    err = 0

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' for '//TRIM(infilename)//' on proc ', my_radar_id

    !=====================================================================================================================
    ! Re-construct the filename series of the PPI-files from the given input filename,
    !  which contains the basic filename, the number of PPI sweeps and the list of elevation identifiers:
    !   (This should work also for ninfiles = 1)
    !=====================================================================================================================

    IF (infilename(1:4) == 'T_PA') THEN

      ! .. This is a file from the OPERA data hub:
    
      ! .. Reconstruct the series of filenames and allocate infilenames(ninfiles):
      CALL reconstruct_file_series_opera (yzroutine, infilename, infilenames, ninfiles, ierr)
      IF (ierr /= 0) THEN
        err = 50
        RETURN
      END IF

      is_native_format = .FALSE.
      provider = infilename(12:15)
      
    ELSE IF (SCAN(TRIM(infilename), 'MP') == 1 .AND. SCAN(TRIM(infilename), 'L') == 2) THEN
      
      ! .. The file starts with an M or P and the second char is an L, so it is a Swiss native h5 file,
      !     for which this reader can also be applied:

      ! .. Reconstruct the series of filenames and allocate infilenames(ninfiles):
      CALL reconstruct_file_series_mch (yzroutine, infilename, infilenames, ninfiles, ierr)
      IF (ierr /= 0) THEN
        err = 50
        RETURN
      END IF

      is_native_format = .TRUE.
      provider = 'LSSW'

    ELSE

      WRITE (*,'(a)') 'ERROR '//TRIM(yzroutine)//': Unknown file format for this reader! ' //&
           'Only OPERA h5 files starting with "T_PA" or Swiss multi-moment single-sweep implemented!'
      err = 10
      RETURN
      
    END IF

    !=====================================================================================================================
    ! Initialize HDF5 FORTRAN interface:
    !=====================================================================================================================

    CALL h5open_f(ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed to initialize the HDF5 fortran interface'
      err = 1
      RETURN
    ENDIF

    !=====================================================================================================================
    ! .. Determine the file type and content by looking at the first file:
    !=====================================================================================================================

    CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(1)), H5F_ACC_RDONLY_F, file_id, ierr)
    IF (ierr /= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '//TRIM(ydirradarin)//TRIM(infilenames(1))
      err = 51
    ENDIF

    ! .. For file "file_id", detect the number of "dataset"s in the root group (nele_obs)
    !    and allocate/build the list of dataset names (dsetname).
    !    This will allocate dsetname. yzroutine and infilename are only for annotation of error messages:
    CALL get_names_of_datasets_odim_h5(file_id, TRIM(yzroutine), TRIM(infilename), '/', 'dataset', &
         LEN(dsetname), nele_obs, dsetname, ierr)
    IF (ierr /= 0) THEN
      err = ierr
    END IF

    ! .. For file "file_id", detect the number of "data"s in the /dataset1 group (ndataset)
    !    and allocate/build the list of vardataset names (vardataset).
    !    This will allocate vardataset. yzroutine and infilename are only for annotation of error messages:
    CALL get_names_of_datasets_odim_h5(file_id, TRIM(yzroutine), TRIM(infilename), '/'//TRIM(dsetname(1)), 'data', &
         LEN(vardataset), ndataset, vardataset, ierr)
    IF (ierr /= 0) THEN
      err = ierr
    END IF

    CALL h5fclose_f(file_id, ierr)

    !=====================================================================================================================
    ! .. Determine the country of origin and the provided datasets by looking at the first file:
    !=====================================================================================================================

    ! .. identifier of the different data sets:
    fids(:) = -1

    CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(1)), H5F_ACC_RDONLY_F, file_id, ierr)
    IF (ierr /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '// &
                     TRIM(ydirradarin)//TRIM(infilenames(1))
      err = 51
    ENDIF

    CALL read_attr_odim (file_id, '/what', 'source', ierr, cval=csource)

    ! .. Get WMO-id and/or NOD (statid_char) from csource:
    !     (If no station_id (=WMO) is in the csource, station_id = 999999 and mii = 99)
    !     (If no statid_char (=NOD) is in the csource, statid_char = 'XXXXX')
    CALL get_stationID_from_HDF5_source(csource, station_id, mii, statid_char, ierr)
    IF (ierr /= 0) THEN
      err = ierr
    END IF

    ! .. Search the actual station in the background list and use it's meta data
    !     as a background list:
    station_found = .FALSE.
    DO i = 1, SIZE(rsmm_bg)
      IF ( station_id == rsmm_bg(i)%station_id ) THEN
        ! .. Copy the base values into rs_meta_ncdf(ista):
        rsmo =  rsm_multitime2onetime ( rsmm_bg(i) )
        station_found = .TRUE.
        EXIT
      END IF
    END DO

    ! loop over all datasets (quantity) of the first dataset group to find out which ones are present:
    DO i=1, ndataset
      cpath(:) = ' '
      cpath    = '/'//TRIM(dsetname(1))//'/'//TRIM(vardataset(i))//'/what'
      CALL read_attr_odim (file_id, TRIM(cpath), 'quantity', ierr, cval=quantity)
      IF (ierr == 0) THEN
        IF (quantity(1:4) == rsmo%obs_hdf5_varname_vrad(1:4)) THEN
          fids(i_vrad) = i
        ELSE IF (TRIM(quantity) == TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
          fids(i_dbzh) = i
        ELSE
          IF (ldebug_radsim) WRITE (*,'(a,a,a)') 'WARNING '//TRIM(yzroutine) // ': quantity=', TRIM(quantity), &
               ' not implemented! Record discarded! '//TRIM(infilename)
        END IF
      ELSE
        err = ierr
        WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed reading '// &
                       TRIM(cpath)//'/quantity'// &
                       ' from HDF5-file '//TRIM(ydirradarin)//TRIM(infilenames(1))
      END IF
    END DO

    CALL h5fclose_f(file_id, ierr)

    IF (err /= 0) THEN
      WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Problem reading the HDF5-file '// &
                     TRIM(ydirradarin)//TRIM(infilenames(1))      
      CALL h5close_f (ierr)
      CALL cleanup_memory
      RETURN
    END IF

    IF (.NOT.station_found) THEN
      err = 3
      WRITE (*,'(i5,1x,a,i7,a)') my_radar_id, &
           TRIM(yzroutine) // 'WARNING '//TRIM(yzroutine)//': radar station ', station_id, &
           ' not found in background list rsmm_bg! file = '//TRIM(ydirradarin)//TRIM(infilename)
      CALL h5close_f (ierr)
      CALL cleanup_memory
      RETURN
    END IF

    !=====================================================================================================================
    ! .. Determine the file format:
    !=====================================================================================================================

    fileformat(:)  = ' '
    is_multimoment = .FALSE.
    is_multisweep  = .FALSE.
    SELECT CASE (nele_obs)
    CASE ( 2 : )
      SELECT CASE (ndataset)
      CASE ( 2 : )
        IF (is_native_format) THEN
          fileformat     = c_format_h5_native_mmms
        ELSE
          fileformat     = c_format_h5_2_opera_mmms
        END IF
        is_multimoment = .TRUE.
        is_multisweep  = .TRUE.
      CASE ( 1 )
        IF (is_native_format) THEN
          fileformat     = c_format_h5_native_smms
        ELSE
          fileformat     = c_format_h5_2_opera_smms
        END IF
        is_multimoment = .FALSE.
        is_multisweep  = .TRUE.
      CASE ( : 0 )
        WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine) // ': ndataset=', ndataset, &
             ', could not determine type of ODIM-HDF5 file ' // TRIM(infilename)
        err = 47
      END SELECT
    CASE ( 1 )
      SELECT CASE (ndataset)
      CASE ( 2 : )
        IF (is_native_format) THEN
          fileformat     = c_format_h5_native_mmss
        ELSE
          fileformat     = c_format_h5_2_opera_mmss
        END IF
        is_multimoment = .TRUE.
        is_multisweep  = .FALSE.
      CASE ( 1 )
        IF (is_native_format) THEN
          fileformat     = c_format_h5_native_smss
        ELSE
          fileformat     = c_format_h5_2_opera_smss
        END IF
        is_multimoment = .FALSE.
        is_multisweep  = .FALSE.
      CASE ( : 0 )
        WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine) // ': ndataset=', ndataset, &
             ', could not determine type of ODIM-HDF5 file ' // TRIM(infilename)
        err = 48
      END SELECT
    CASE ( : 0)
      WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine) // ': nele_obs=', nele_obs, &
           ', could not determine type of ODIM-HDF5 file ' // TRIM(infilename)
      err = 49
    END SELECT

    IF (err /= 0) THEN
      CALL h5close_f (ierr)
      CALL cleanup_memory
      RETURN
    END IF

    !=====================================================================================================================
    ! .. For multisweep files, filter elevations using the default scan strategies as a whitelist. If an elevation is not
    !    in a default strategy, ignore it. In this way it is possible to ignore certain elevations if
    !    for example, the ra_inc differs from the other elevations (Netherlands, Czech Republic).
    !=====================================================================================================================

    IF (is_multisweep) THEN

      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(1)), H5F_ACC_RDONLY_F, file_id, ierr)

      nele = 0
      SELECT CASE (rsmo%icountry)
      CASE (i_netherlands)
        
        DO i=1, nele_obs
          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(i))//'/where', 'rscale', ierr, dval=ra_inc, silent=.TRUE.)
          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(i))//'/where', 'elangle', ierr, dval=dval, silent=.TRUE.)
          ! .. Accept elevation in a multi-sweep file only if it is in the list of the known nominal default strategies for the actual country:
          IF ( ele_is_present_in_default(dval, rsmo) .AND. ABS(ra_inc-rsmo%ra_inc) <= 1E-1_dp ) THEN
            nele = nele + 1
            dsetname(nele) = dsetname(i)
          END IF
        END DO
        DO i= nele+1, SIZE(dsetname)
          dsetname(i)(:) = ' '
        END DO

      CASE default

        DO i=1, nele_obs
          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(i))//'/where', 'elangle', ierr, dval=dval, silent=.TRUE.)
          ! .. Accept elevation in a multi-sweep file only if it is in the list of the known nominal default strategies for the actual country:
          IF ( ele_is_present_in_default(dval, rsmo) ) THEN
            nele = nele + 1
            dsetname(nele) = dsetname(i)
          END IF
        END DO
        DO i= nele+1, SIZE(dsetname)
          dsetname(i)(:) = ' '
        END DO

      END SELECT
      
      CALL h5fclose_f(file_id, ierr)

      IF (nele <= 0) THEN
        WRITE (*,'(a,i0,a)') 'WARNING '//TRIM(yzroutine) // ': nele=', nele, &
             ', no matching elevations found in ' // TRIM(infilename)
        err = 78
        CALL h5close_f (ierr)
        CALL cleanup_memory
        RETURN
      END IF

    ELSE  
      nele = nele_obs
    END IF

    !=====================================================================================================================
    ! .. Initializations for the below loop over all files of the volume scan:
    !=====================================================================================================================

    tmpdim = MAX(ninfiles, nele)
    ALLOCATE (manel(tmpdim));         manel(:)         = miss_value
    ALLOCATE (v_nyq(tmpdim));         v_nyq(:)         = miss_value
    ALLOCATE (high_prf(tmpdim));      high_prf(:)      = miss_value
    ALLOCATE (dualprf_ratio(tmpdim)); dualprf_ratio(:) = miss_value
    cdate(:) = ' '
    cdate_last = 'YYYYMMDDhhmmss'
    quantity_last(:) = 'X'
    station_id_last = 999999
    naz = missval_int
    nra = missval_int
    ra_inc = miss_value
    alt_msl_true = miss_value
    lat = miss_value
    lon = miss_value
    az_inc = miss_value
    az_start = miss_value
    az_start_last = miss_value


    !=====================================================================================================================
    ! .. Loop over all files of the volume scan:
    !=====================================================================================================================

    ii = 0 ! Index over datasets
    file_loop: DO iii=1, ninfiles

      ! Open the file for the particular timestep
      CALL h5fopen_f(TRIM(ydirradarin)//TRIM(infilenames(iii)), H5F_ACC_RDONLY_F, file_id, ierr)
      IF (ierr /= 0) THEN
        WRITE(*,'(a)') 'WARNING '//TRIM(yzroutine)//': Failed opening the HDF5-file '// &
                       TRIM(ydirradarin)//TRIM(infilenames(iii))
        err = ierr
      ENDIF

      ! .. Accept elevation in a single-sweep file only if it is in the list of the known nominal default strategies for the actual country:
      IF (.NOT. is_multisweep) THEN
        CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/where', 'elangle', ierr, dval=dval, silent=.TRUE.)
        IF (ierr /= 0) THEN
          ele_found = .FALSE.
        ELSE
          ele_found = ele_is_present_in_default(dval, rsmo)
        END IF
      ELSE
        ele_found = .TRUE.  ! dummy
      END IF

      IF (is_multisweep .OR. (.NOT. is_multisweep .AND. ele_found) ) THEN

        ii = ii + 1

!!$ File content as seen by "h5dump -n 1 T_PAGB52_C_LSSW_20220104101500.hdf":
!!$HDF5 "T_PAGB52_C_LSSW_20220104101500.hdf" {
!!$FILE_CONTENTS {
!!$ group      /
!!$ attribute  /Conventions
!!$ group      /dataset1
!!$ group      /dataset1/data1
!!$ dataset    /dataset1/data1/data
!!$ attribute  /dataset1/data1/data/CLASS
!!$ attribute  /dataset1/data1/data/IMAGE_VERSION
!!$ group      /dataset1/data1/how
!!$ group      /dataset1/data1/what
!!$ attribute  /dataset1/data1/what/gain
!!$ attribute  /dataset1/data1/what/nodata
!!$ attribute  /dataset1/data1/what/offset
!!$ attribute  /dataset1/data1/what/quantity
!!$ attribute  /dataset1/data1/what/undetect
!!$ group      /dataset1/data2
!!$ dataset    /dataset1/data2/data
!!$ attribute  /dataset1/data2/data/CLASS
!!$ attribute  /dataset1/data2/data/IMAGE_VERSION
!!$ group      /dataset1/data2/how
!!$ group      /dataset1/data2/what
!!$ attribute  /dataset1/data2/what/gain
!!$ attribute  /dataset1/data2/what/nodata
!!$ attribute  /dataset1/data2/what/offset
!!$ attribute  /dataset1/data2/what/quantity
!!$ attribute  /dataset1/data2/what/undetect
!!$ group      /dataset1/data3
!!$ dataset    /dataset1/data3/data
!!$ attribute  /dataset1/data3/data/CLASS
!!$ attribute  /dataset1/data3/data/IMAGE_VERSION
!!$ group      /dataset1/data3/how
!!$ group      /dataset1/data3/what
!!$ attribute  /dataset1/data3/what/gain
!!$ attribute  /dataset1/data3/what/nodata
!!$ attribute  /dataset1/data3/what/offset
!!$ attribute  /dataset1/data3/what/quantity
!!$ attribute  /dataset1/data3/what/undetect
!!$ group      /dataset1/how
!!$ attribute  /dataset1/how/NI
!!$ attribute  /dataset1/how/startazA
!!$ attribute  /dataset1/how/stopazA
!!$ group      /dataset1/what
!!$ attribute  /dataset1/what/enddate
!!$ attribute  /dataset1/what/endtime
!!$ attribute  /dataset1/what/product
!!$ attribute  /dataset1/what/startdate
!!$ attribute  /dataset1/what/starttime
!!$ group      /dataset1/where
!!$ attribute  /dataset1/where/a1gate
!!$ attribute  /dataset1/where/elangle
!!$ attribute  /dataset1/where/nbins
!!$ attribute  /dataset1/where/nrays
!!$ attribute  /dataset1/where/rscale
!!$ attribute  /dataset1/where/rstart
!!$ group      /how
!!$ attribute  /how/RXlossH
!!$ attribute  /how/RXlossV
!!$ attribute  /how/TXlossH
!!$ attribute  /how/TXlossV
!!$ attribute  /how/antgainH
!!$ attribute  /how/antgainV
!!$ attribute  /how/antspeed
!!$ attribute  /how/beamwH
!!$ attribute  /how/beamwV
!!$ attribute  /how/beamwidth
!!$ attribute  /how/highprf
!!$ attribute  /how/lowprf
!!$ attribute  /how/pulsewidth
!!$ attribute  /how/radconstH
!!$ attribute  /how/radconstV
!!$ attribute  /how/radomelossH
!!$ attribute  /how/radomelossV
!!$ attribute  /how/rpm
!!$ attribute  /how/scan_index
!!$ attribute  /how/system
!!$ attribute  /how/wavelength
!!$ group      /what
!!$ attribute  /what/date
!!$ attribute  /what/object
!!$ attribute  /what/source
!!$ attribute  /what/time
!!$ attribute  /what/version
!!$ group      /where
!!$ attribute  /where/height
!!$ attribute  /where/lat
!!$ attribute  /where/lon
!!$ }
!!$}


        ! .. read some attributes from different groups/datasets. The choice of the
        !    optional arguments has to match the data type of the attribute in the
        !    file:  dval = double/real/float   ,  ival = integer  ,  cval = character string
        !    A type mismatch might lead to a runtime crash!

        CALL read_attr_odim (file_id, '/where', 'height', ierr, dval=dval)
        IF (ii > 1 .AND. ABS(dval-alt_msl_true) > 1E-6_dp) THEN
          WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
               ': alt_msl_true differs from previous PPI file in ' // &
               TRIM(infilenames(iii))//': ', dval, ' /= ', alt_msl_true
          err = 6
        END IF
        alt_msl_true = dval

        CALL read_attr_odim (file_id, '/where', 'lat', ierr, dval=dval)
        IF (ii > 1 .AND. ABS(dval-lat) > 1E-6_dp) THEN
          WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
               ': lat differs from previous PPI file in ' // &
               TRIM(infilenames(iii))//': ', dval, ' /= ', lat
          err = 7
        END IF
        lat         = dval

        CALL read_attr_odim (file_id, '/where', 'lon', ierr, dval=dval)
        IF (ii > 1 .AND. ABS(dval-lon) > 1E-6_dp) THEN
          WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
               ': lon differs from previous PPI file in ' // &
               TRIM(infilenames(iii))//': ', dval, ' /= ', lon
          err = 8
        END IF
        lon         = dval

        IF (is_multisweep) THEN

          ! .. Loop over all datasets to read and cross-check all attributes:
          DO j=1, nele

            CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'nbins', ierr, ival=ival)
            nra = MAX(INT(ival), nra)

            CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'nrays', ierr, ival=ival)
            IF (provider == 'KITC' .AND. rsmo%icountry == i_dwd) THEN

              ! Workaround for problematic files from KIT-Cube which may contain different numbers of azimuts
              ! for the same elevations of different variables, and/or different numbers of azimuts for different
              ! elevations. Mostly this is one more azimut (361) than nominal (360). As this happens
              ! very often and these data are not used operationally, those files are not discarded right away,
              ! but in the data reader we will later use a1gate to separate out the superfluous datum.

!!$              IF (j == 1) THEN
!!$                naz_last = ival
!!$                naz      = MIN(INT(ival), rsmo%naz)
!!$              ELSE IF (j > 1 .AND. ival /= naz_last) THEN
!!$                WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
!!$                     ': naz differs from previous dataset in file ' // &
!!$                     TRIM(infilename)//': ', ival, ' /= ', naz_last
!!$                naz_last = ival
!!$                IF (ival >= rsmo%naz) THEN
!!$                  naz = MIN(INT(ival), naz)
!!$                ELSE
!!$                  err = 9
!!$                END IF
!!$              END IF
              
              ! simplified:
              naz = rsmo%naz
              
            ELSE
              
              IF (j == 1) THEN
                naz_last = ival
              ELSE IF (j > 1 .AND. ival /= naz_last) THEN
                WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
                     ': naz differs from previous dataset in file ' // &
                     TRIM(infilename)//': ', ival, ' /= ', naz_last
                err = 9
                naz_last = ival
              END IF
              naz = MAX(INT(ival), naz)

            END IF

            CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'rscale', ierr, dval=dval)
            IF (j > 1 .AND. ABS(dval-ra_inc) > 1E-1_dp) THEN
              WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
                   ': ra_inc differs from previous dataset in file ' // &
                   TRIM(infilename)//': ', dval, ' /= ', ra_inc
              err = 10
            END IF
            ra_inc = REAL(NINT(dval*10.0_dp)*0.1_dp, kind=dp)

            CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'elangle', ierr, dval=manel(j))

            !          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/where', 'rstart', ierr, dval=dval)  ! rstart = ra_inc!

            IF (naz > 0) THEN
              az_inc = 360.0_dp / naz
              ! just in case there are few more naz in the file than nominal, round az_inc to the next 0.02 degrees:
              az_inc = REAL(NINT(az_inc/0.02_dp)*0.02_dp, kind=dp)
            ELSE
              az_inc = 1.0_dp
            END IF

            ! Nominal first azimut is chosen as half the azimut increment.
            ! This is always the best rounded value of the true azimut of the first ray in the data,
            !  because for ODIM-HDF5, the first true antenna azimut with a value > 0.0 is stored at the first
            !  index.
            az_start = 0.5_dp*az_inc

            DO i=1, ndataset

              tmpgroup(:) = ' '
              tmpgroup    = '/'//TRIM(dsetname(j))//'/'//TRIM(vardataset(i))

              ! .. Check that "quantity" are consistent among "/datasetXX/dataYY":
              quantity(:) = ' '
              CALL read_attr_odim (file_id, TRIM(tmpgroup)//'/what', 'quantity', ierr, cval=quantity)
              IF (fids(i_vrad) == i .AND. quantity(1:4) /= rsmo%obs_hdf5_varname_vrad(1:4)) THEN
                WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                     ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                     TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_vrad)
                err = 12
              END IF
              IF (fids(i_dbzh) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
                WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                     ': quantity in file differs from that of the /'//TRIM(dsetname(1))//' in ' // &
                     TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_dbzh)
                err = 13
              END IF

              ! Find out the Nyquist velocity for radial wind data sets:
              !  In newer files there is an attribute '/dataset1/how/NI' (Nyquist Interval)
              !  If that is not found, we look at '/how/NI', and otherwise we resort
              ! to the negative offset of the VRAD dataset.

              IF (fids(i_vrad) == i .AND. quantity(1:4) == rsmo%obs_hdf5_varname_vrad(1:4)) THEN

                ! First try: Nyquist-velocity from /dataset(j)/data(i)/how:
                CALL read_attr_odim (file_id, TRIM(tmpgroup)//'/how', 'NI', ierr, dval=dval, silent=.TRUE.)
                IF (ierr == 0) THEN
                  v_nyq(j) = dval
                ELSE
                  ! Second try: from /dataset(j)//how:
                  CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'NI', ierr, dval=dval, silent=.TRUE.)
                  IF (ierr == 0) THEN
                    v_nyq(j) = dval
                  ELSE
                    ! Third try: from /how:
                    CALL read_attr_odim (file_id, '/how', 'NI', ierr, dval=dval, silent=.TRUE.)
                    IF (ierr == 0) THEN
                      v_nyq(j) = dval
                    ELSE
                      cpath(:) = ' '
                      cpath = '/'//TRIM(dsetname(j))//'/'//TRIM(vardataset(i))//'/what'

                      CALL read_attr_odim (file_id, TRIM(cpath), 'offset', ierr, dval=dval)
                      v_nyq(j) = -dval
                    END IF
                  END IF
                END IF

                ! First try:
                CALL read_attr_odim (file_id, TRIM(tmpgroup)//'/how', 'highprf', ierr, dval=dval, silent=.TRUE.)
                IF (ierr == 0) THEN
                  high_prf(j) = dval
                  CALL read_attr_odim (file_id, TRIM(tmpgroup)//'/how', 'lowprf', ierr, dval=dval, silent=.TRUE.)
                  IF (ierr == 0) THEN
                    dualprf_ratio(j) = high_prf(j) / dval
                  ELSE
                    dualprf_ratio(j) = 1.0_dp
                  END IF
                ELSE
                  ! Second try:
                  CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'highprf', ierr, dval=dval, silent=.TRUE.)
                  IF (ierr == 0) THEN
                    high_prf(j) = dval
                    CALL read_attr_odim (file_id, '/'//TRIM(dsetname(j))//'/how', 'lowprf', ierr, dval=dval, silent=.TRUE.)
                    IF (ierr == 0) THEN
                      dualprf_ratio(j) = high_prf(j) / dval
                    ELSE
                      dualprf_ratio(j) = 1.0_dp
                    END IF
                  ELSE
                    ! Third try:
                    CALL read_attr_odim (file_id, '/how', 'highprf', ierr, dval=dval, silent=.TRUE.)
                    IF (ierr == 0) THEN
                      high_prf(j) = dval
                      CALL read_attr_odim (file_id, '/how', 'lowprf', ierr, dval=dval, silent=.TRUE.)
                      IF (ierr == 0) THEN
                        dualprf_ratio(j) = high_prf(j) / dval
                      ELSE
                        dualprf_ratio(j) = 1.0_dp
                      END IF
                    ELSE
                      dualprf_ratio(j) = 1.0_dp
                    END IF
                  END IF
                END IF

              END IF

            END DO
          END DO

        ELSE

          ! This is a single-sweep file, either single- or multi-moment:

          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/where', 'nbins', ierr, ival=ival)
          nra = MAX(INT(ival), nra)

          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/where', 'nrays', ierr, ival=ival)
          IF (ii == 1) THEN
            naz_last = ival
          ELSE IF (ii > 1 .AND. ival /= naz_last) THEN
            WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
                 ': naz differs from previous PPI file in ' // &
                 TRIM(infilenames(iii))//': ', ival, ' /= ', naz_last
            err = 9
            naz_last = ival
          END IF
          naz = MAX(INT(ival), naz)

          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/where', 'rscale', ierr, dval=dval)
          IF (ii > 1 .AND. ABS(dval-ra_inc) > 1E-1_dp) THEN
            WRITE (*,'(a,f7.1,a,f7.1)') 'WARNING '//TRIM(yzroutine)// &
                 ': ra_inc differs from previous PPI file in ' // &
                 TRIM(infilenames(iii))//': ', dval, ' /= ', ra_inc
            err = 13
          END IF
          ra_inc = REAL(NINT(dval*10.0_dp)*0.1_dp, kind=dp)

          CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/where', 'elangle', ierr, dval=manel(ii))

          IF (naz > 0) THEN
            az_inc = 360.0_dp / naz
            ! just in case there are few more naz in the file than nominal, round az_inc to the next 0.02 degrees:
            az_inc = REAL(NINT(az_inc/0.02_dp)*0.02_dp, kind=dp)
          ELSE
            az_inc = 1.0_dp
          END IF

          ! Nominal first azimut is chosen as half the azimut increment.
          ! This is always the best rounded value of the true azimut of the first ray in the data,
          !  because for ODIM-HDF5, the first true antenna azimut with a value > 0.0 is stored at the first
          !  index.
          az_start = 0.5_dp*az_inc

          ! Check if the variables in the files are consistent to the first file in the series:
          DO i=1, ndataset
            cpath(:) = ' '
            cpath = '/'//TRIM(dsetname(1))//'/'//TRIM(vardataset(i))
            quantity(:) = ' '
            CALL read_attr_odim (file_id, TRIM(cpath)//'/what', 'quantity', ierr, cval=quantity)
            IF (fids(i_vrad) == i .AND. quantity(1:4) /= rsmo%obs_hdf5_varname_vrad(1:4)) THEN
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                   ': quantity in file differs from that of the '//TRIM(cpath)//'/what in ' // &
                   TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_vrad)
              err = 12
            END IF
            IF (fids(i_dbzh) == i .AND. TRIM(quantity) /= TRIM(rsmo%obs_hdf5_varname_dbzh)) THEN
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
                   ': quantity in file differs from that of the '//TRIM(cpath)//'/what in ' // &
                   TRIM(infilename)//': ' // TRIM(quantity) // ' /= ' // TRIM(rsmo%obs_hdf5_varname_dbzh)
              err = 13
            END IF


            ! Find out the Nyquist velocity for radial wind data sets:
            !  In newer files there is an attribute '/dataset1/how/NI' (Nyquist Interval)
            !  If that is not found, we look at '/how/NI', and otherwise we resort
            ! to the negative offset of the VRAD dataset.

            IF (fids(i_vrad) == i .AND. quantity(1:4) == rsmo%obs_hdf5_varname_vrad(1:4)) THEN

              ! First try: Nyquist-velocity from /dataset1/data(i)/how:
              CALL read_attr_odim (file_id, TRIM(cpath)//'/how', 'NI', ierr, dval=dval, silent=.TRUE.)
              IF (ierr == 0) THEN
                v_nyq(ii) = dval
              ELSE
                ! Second try: Nyquist-velocity from /dataset1/how:
                CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/how', 'NI', ierr, dval=dval, silent=.TRUE.)
                IF (ierr == 0) THEN
                  v_nyq(ii) = dval
                ELSE
                  ! Third try: from /how:
                  CALL read_attr_odim (file_id, '/how', 'NI', ierr, dval=dval, silent=.TRUE.)
                  IF (ierr == 0) THEN
                    v_nyq(ii) = dval
                  ELSE
                    CALL read_attr_odim (file_id, TRIM(cpath)//'/what', 'offset', ierr, dval=dval)
                    v_nyq(ii) = -dval
                  END IF
                END IF
              END IF

              ! First try:
              CALL read_attr_odim (file_id, TRIM(cpath)//'/how', 'highprf', ierr, dval=dval, silent=.TRUE.)
              IF (ierr == 0) THEN
                high_prf(ii) = dval
                CALL read_attr_odim (file_id, TRIM(cpath)//'/how', 'lowprf', ierr, dval=dval, silent=.TRUE.)
                IF (ierr == 0) THEN
                  dualprf_ratio(ii) = high_prf(ii) / dval
                ELSE
                  dualprf_ratio(ii) = 1.0_dp
                END IF
              ELSE
                ! Second try:
                CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/how', 'highprf', ierr, dval=dval, silent=.TRUE.)
                IF (ierr == 0) THEN
                  high_prf(ii) = dval
                  CALL read_attr_odim (file_id, '/'//TRIM(dsetname(1))//'/how', 'lowprf', ierr, dval=dval, silent=.TRUE.)
                  IF (ierr == 0) THEN
                    dualprf_ratio(ii) = high_prf(ii) / dval
                  ELSE
                    dualprf_ratio(ii) = 1.0_dp
                  END IF
                ELSE
                  ! Third try: from /how:
                  CALL read_attr_odim (file_id, '/how', 'highprf', ierr, dval=dval, silent=.TRUE.)
                  IF (ierr == 0) THEN
                    high_prf(ii) = dval
                    CALL read_attr_odim (file_id, '/how', 'lowprf', ierr, dval=dval, silent=.TRUE.)
                    IF (ierr == 0) THEN
                      dualprf_ratio(ii) = high_prf(ii) / dval
                    ELSE
                      dualprf_ratio(ii) = 1.0_dp
                    END IF
                  ELSE
                    dualprf_ratio(ii) = 1.0_dp
                  END IF
                END IF
              END IF

            END IF

          END DO

        END IF

        CALL read_attr_odim (file_id, '/what', 'source', ierr, cval=cval)
        CALL get_stationID_from_HDF5_source (cval, station_id, mii, statid_char, ierr)
        IF (ierr /= 0) THEN
          err = 14
        END IF
        IF (ii > 1 .AND. station_id /= station_id_last) THEN
          WRITE (*,'(a,i7,a,i7)') 'WARNING '//TRIM(yzroutine)// &
               ': station_id differs from previous PPI file in ' // &
               TRIM(infilenames(iii))//': ', station_id, ' /= ', station_id_last
          err = 15
        END IF
        station_id_last = station_id


        ! station_name: overtake later from background list

        ! date and time:
        CALL read_attr_odim (file_id, '/what', 'date', ierr, cval=cval)
        IF ( ierr /= 0) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': date could not be read from file ' // &
               TRIM(infilename)
          err = 16
        END IF
        cdate(:) = ' '
        cdate = TRIM(cval)

        CALL read_attr_odim (file_id, '/what', 'time', ierr, cval=cval)
        IF ( ierr /= 0) THEN
          WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)// &
               ': time could not be read from file ' // &
               TRIM(infilenames(iii))
          err = 17
        END IF
        cdate = TRIM(cdate)//TRIM(cval)

        SELECT CASE (rsmo%icountry)
        CASE (i_meteoswiss, i_dwd, i_france, i_arpasim, i_belgium, i_denmark, i_poland, i_czech, i_netherlands, i_slovakia)
          ! .. cdate is start time or end time incl. seconds, round down to the next full dt_obs minutes
          !     to get the start time without seconds:
          IF (rsmo%dt_obs(3) > 0.0_dp) THEN
            ! assuming new notation, where the time increment is in dt_obs(3):
            cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(3)/60))
          ELSE IF (rsmo%dt_obs(1) > 0.0_dp) THEN
            ! assuming old notation, where the time increment is in dt_obs(1) and (2), (3) are miss_value:
            cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(1)/60))
          ELSE
            ! This should not happen, except that the user give explicitly a wrong dt_obs-triplet in the namelist:
            WRITE (*,'(a,3(1x,f0.1))') 'WARNING '//TRIM(yzroutine)// &
                 ': specification of dt_obs(1:3) is wrong for station in file ' // &
                 TRIM(infilenames(iii))//': dt_obs =', rsmo%dt_obs
            err = 18
          END IF
        CASE default
          ! .. assume the same for other countries, too:
          IF (rsmo%dt_obs(3) > 0.0_dp) THEN
            ! assuming new notation, where the time increment is in dt_obs(3):
            cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(3)/60))
          ELSE IF (rsmo%dt_obs(1) > 0.0_dp) THEN
            ! assuming old notation, where the time increment is in dt_obs(1) and (2), (3) are miss_value:
            cdate = round_datestr_min(cdate, -NINT(rsmo%dt_obs(1)/60))
          ELSE
            ! This should not happen, except that the user give explicitly a wrong dt_obs-triplet in the namelist:
            WRITE (*,'(a,3(1x,f0.1))') 'WARNING '//TRIM(yzroutine)// &
                 ': specification of dt_obs(1:3) is wrong for station in file ' // &
                 TRIM(infilenames(iii))//': dt_obs =', rsmo%dt_obs
            err = 18
          END IF
        END SELECT

        IF (ii > 1 .AND. cdate /= cdate_last) THEN
          WRITE (*,'(a,a,a,a)') 'WARNING '//TRIM(yzroutine)// &
               ': date/time differs from previous PPI file in ' // &
               TRIM(infilenames(iii))//': ', cdate, ' /= ', cdate_last
          err = 19
        END IF
        cdate_last = cdate

      END IF  ! is_multisweep .OR. (.NOT. is_multisweep .AND. ele_is_found )

      ! .. Close file:
      CALL h5fclose_f(file_id, ierr)

    END DO file_loop     ! Loop over ninfiles

    ! .. Close FORTRAN interface:
    CALL h5close_f (ierr)

    ! .. Re-set ninfiles and tmpdim:
    ninfiles = ii
    tmpdim = MAX(ninfiles, nele)

    ! Fix specialities of some countries:
    SELECT CASE (rsmo%icountry)
    CASE (i_denmark)
      ! Denmark has alternating volume scans with ~240 and ~480 range bins.
      ! We take a constant 480 in order to be able to digest these data.
      nra = 480
      ! Also it has alternating 8.4 and 8.5 elevations, which is too similar
      !  compared to the 0.1 deg tolerance. So we set 8.4 to 8.5:
      DO i=1, tmpdim
        IF (ABS(manel(i) - 8.4_dp) <= 0.5_dp) THEN
          manel(i) = 8.5_dp
        END IF
      END DO
    END SELECT

    ! .. Now partially override the station meta data with the meta data from the file:
    !     (rsm%station_name is taken from rsmm_bg)
    rsmo%station_id         = station_id
    rsmo%obsfile_format(1)  = fileformat
    IF (is_multisweep) THEN
      rsmo%nel              = nele
    ELSE
      rsmo%nel              = ninfiles
    END IF
    rsmo%alt_msl_true     = alt_msl_true
    rsmo%ra_inc           = ra_inc
    rsmo%ra_inc_obs       = ra_inc
    rsmo%az_start         = az_start
    rsmo%az_inc           = az_inc
    rsmo%naz              = naz
    rsmo%nra              = nra
    rsmo%nra_obs          = nra
    rsmo%lon              = lon
    rsmo%lat              = lat
    IF (is_multimoment) THEN
      ! all variables are in one file:
      DO i=1, ndatakind
        IF (fids(i) > 0) THEN
          rsmo%obsfile(:,i) = TRIM(infilename)//'/'//TRIM(vardataset(fids(i)))
        END IF
      END DO
    ELSE
      ! each variable is in a different file:
      DO i=1, ndatakind
        IF (fids(i) > 0) THEN
          rsmo%obsfile(:,i) = TRIM(infilename)
        END IF
      END DO
    END IF
    DO i=1, ndatakind
      IF (fids(i) > 0) THEN
        IF (is_multisweep) THEN
          rsmo%nrep_ncdf(i)     = nele
          rsmo%naz_ncdf(:,i)    = naz
        ELSE
          rsmo%nrep_ncdf(i)     = ninfiles
          rsmo%naz_ncdf(:,i)    = naz
        END IF
      END IF
    END DO


    ! Extract blocks of obstimes from cdate; the only possible error
    !  status of this subroutine (cdate not valid)  equals err_ista(ista) = 5
    !  from above, which has already been checked:
    fidindlist = -1
    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        fidindlist(i) = i
      END IF
    END DO

    CALL put_obstimes_to_rsm_onetime ( rsmo, TRIM(infilename), fidindlist, ydate_ini_mod, cdate, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      CALL cleanup_memory
      RETURN
    END IF

    ! .. recognizes also DWD precip scans:
    CALL get_scanstrategy_ID ( rsmo, TRIM(infilename), manel(1:tmpdim), 0.1_dp, ierr )
    IF (ierr /= 0) THEN
      err = ierr
      CALL cleanup_memory
      RETURN
    ELSE
      rsmo%ext_nyq(:,:) = unused_value
      rsmo%high_nyq(:)  = unused_value
      rsmo%prf(:)       = unused_value
      IF (fids(i_vrad) > 0) THEN
        ! .. Radial wind data are present, so store Nyquist velocity information:
        DO j = 1, tmpdim
          DO i = 1, rsmo%nel
            IF (ABS(manel(j)-rsmo%el_arr(i)) <= 0.1_dp) THEN
              rsmo%ext_nyq(i,:) = v_nyq(j)
              rsmo%high_nyq(i)  = v_nyq(j)
              rsmo%prf(i)       = NINT(4.0d0*v_nyq(j)/rsmo%lambda) ! Compute from v_nyq and wavelength
            END IF
          END DO
        END DO
      END IF
    END IF

    ! .. Clean up memory:
    CALL cleanup_memory

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    SUBROUTINE cleanup_memory
      IF (ALLOCATED(manel))         DEALLOCATE(manel)
      IF (ALLOCATED(v_nyq))         DEALLOCATE(v_nyq)
      IF (ALLOCATED(high_prf))      DEALLOCATE(high_prf)
      IF (ALLOCATED(dualprf_ratio)) DEALLOCATE(dualprf_ratio)
      IF (ALLOCATED(dsetname))      DEALLOCATE(dsetname)
      IF (ALLOCATED(vardataset))    DEALLOCATE(vardataset)
    END SUBROUTINE cleanup_memory

    ! .. Check whether the elevation angle elev is contained in any of the default nominal scan strategies:
    FUNCTION ele_is_present_in_default(elev, rsmo) RESULT (l)

      REAL(kind=dp), INTENT(in) :: elev
      TYPE(radar_meta_type_onetime) :: rsmo
      LOGICAL :: l

      INTEGER :: k

      l = .FALSE.
      DO k=1, nscanstrategies_max
        IF (rsmo%nel_default(k) > 0) THEN
          IF ( ANY( ABS(rsmo%el_arr_default(1:rsmo%nel_default(k),k) - elev) <= 0.1_dp ) ) THEN
            l = .TRUE.
          END IF
        END IF
      END DO

    END FUNCTION ele_is_present_in_default

  END SUBROUTINE get_metadata_from_h5_opera

  !==============================================================================

  SUBROUTINE get_names_of_datasets_odim_h5(file_id, yzroutine, infilename, rootname, dsetbase, strlen, &
                                           ndatasets, dsetnames, err)

    INTEGER(hid_t),   INTENT(in)  :: file_id
    CHARACTER(len=*), INTENT(in)  :: yzroutine, infilename ! to annotate error/warning messages
    CHARACTER(len=*), INTENT(in)  :: rootname ! root of the tree to search for datasets
    CHARACTER(len=*), INTENT(in)  :: dsetbase ! basename of the datasets under root
    INTEGER,          INTENT(in)  :: strlen
    INTEGER,          INTENT(out) :: ndatasets
    CHARACTER(len=strlen), ALLOCATABLE, INTENT(out) :: dsetnames(:)
    INTEGER,          INTENT(out) :: err

    INTEGER          :: i, storage_type, max_corder, nrootmemb, obj_typ, ierr
    INTEGER(hid_t)   :: root_id, grp_id
    CHARACTER(len=3) :: cnumber

    err = 0
    
    CALL h5eset_auto_f(0, ierr)  ! Switch off HDF5 error messages temporarily,
                                 !  because the below loops will deliberately create
                                 !  errors when querying non-existent datasets

    ! .. Open the root group, below which the datasets reside:
    CALL h5gopen_f(file_id, TRIM(rootname), root_id, ierr)
    ! .. Get the number of members in the root group, which are candidates for datasets:
    CALL h5gget_info_f(root_id, storage_type, nrootmemb, max_corder, ierr)

    IF (ierr /= 0 .OR. nrootmemb <= 0) THEN
      WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//': No '//TRIM(dsetbase)//' found in HDF5-file(s) '//TRIM(infilename)
      err = 17
    END IF

    ALLOCATE (dsetnames(MAX(nrootmemb, 1)))
    dsetnames(:)(:)= ' '

    ! .. search for groups named dataset1 ... dataset<nrootmemb> and store
    !    the names of the existing dataset-groups, which are fewer than
    !    nrootmemb:
    ndatasets = 0
    DO i = 1, nrootmemb
      WRITE(cnumber, '(i3)') i
      dsetnames(i) = TRIM(dsetbase)//TRIM(ADJUSTL(cnumber))
      CALL h5oopen_f(root_id, TRIM(dsetnames(i)), grp_id, ierr)
      IF (ierr == 0) THEN
        CALL h5iget_type_f(grp_id, obj_typ, ierr)
        IF (obj_typ == H5I_GROUP_F) THEN
          ndatasets = ndatasets + 1
          dsetnames(ndatasets) = dsetnames(i)
        END IF
      END IF
      CALL h5oclose_f(grp_id, ierr)
    ENDDO
    ! .. Close the root group:
    CALL h5gclose_f(root_id, ierr)

    CALL h5eset_auto_f(1, ierr)  ! Switch on again HDF5 error messages

  END SUBROUTINE get_names_of_datasets_odim_h5
  
  SUBROUTINE read_attr_odim (file_id, grpname, attrname, ierr, ival, ivec, rval, dval, rvec, dvec, cval, cval_varlen, silent)

    ! .. Input/Output variables:
    INTEGER(hid_t),   INTENT(in)  :: file_id
    CHARACTER(len=*), INTENT(in)  :: grpname     ! full path to group/dataset
    CHARACTER(len=*), INTENT(in)  :: attrname    ! attribute name within group
    INTEGER,          INTENT(out) :: ierr
    LOGICAL,          INTENT(in),    OPTIONAL :: silent
    INTEGER,          INTENT(out),   OPTIONAL :: ival
    REAL(KIND=dp),    INTENT(out),   OPTIONAL :: rval  ! for 32 bit reals, but returned as 64 bit
    REAL(KIND=dp),    INTENT(out),   OPTIONAL :: dval
    CHARACTER(len=*), INTENT(inout), OPTIONAL :: cval         ! for fixed length strings
    CHARACTER(len=*), INTENT(inout), OPTIONAL :: cval_varlen  ! for variable length strings
    INTEGER,          INTENT(out),  ALLOCATABLE, OPTIONAL :: ivec(:)
    REAL(KIND=dp),    INTENT(out),  ALLOCATABLE, OPTIONAL :: rvec(:) ! for 32 bit reals, but returned as 64 bit
    REAL(KIND=dp),    INTENT(out),  ALLOCATABLE, OPTIONAL :: dvec(:)

    ! .. Local variables:
    INTEGER            :: ierror, attr_data_type, obj_typ, ntrue, ndims
    INTEGER(hid_t)     :: grp_id, attr_id, type_id, dsp_id, h5_type
    CHARACTER(len=500) :: catt_value
    INTEGER(hsize_t)   :: attr_dims(1)
    INTEGER(hsize_t)   :: size_cattr
    REAL(KIND=sp)      :: ratt_value
    REAL(KIND=dp)      :: datt_value
    INTEGER            :: iatt_value
    INTEGER(hsize_t), ALLOCATABLE :: datadims(:), maxdatadims(:)
    REAL(KIND=sp),    ALLOCATABLE :: ratt_vec(:)
    LOGICAL            :: issilent
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata ! Read buffer for variable length string
    TYPE(C_PTR) :: f_ptr

    ierr = 0

    attr_data_type = 0
    ntrue = 0
    IF (PRESENT(ival)) THEN
      ntrue = ntrue + 1
      attr_data_type = 1
    END IF
    IF (PRESENT(dval)) THEN
      ntrue = ntrue + 1
      attr_data_type = 2
    END IF
    IF (PRESENT(cval)) THEN
      ntrue = ntrue + 1
      attr_data_type = 3
    END IF
    IF (PRESENT(ivec)) THEN
      ntrue = ntrue + 1
      attr_data_type = 4
    END IF
    IF (PRESENT(dvec)) THEN
      ntrue = ntrue + 1
      attr_data_type = 5
    END IF
    IF (PRESENT(cval_varlen)) THEN
      ntrue = ntrue + 1
      attr_data_type = 6
    END IF
    IF (PRESENT(rval)) THEN
      ntrue = ntrue + 1
      attr_data_type = 7
    END IF
    IF (PRESENT(rvec)) THEN
      ntrue = ntrue + 1
      attr_data_type = 8
    END IF
    IF (ntrue > 1) THEN
      ! More than one output variable present, but this is not allowed.
      !  read_attr_odim() should be called for one specific output datum and data type only!
      WRITE(*,'(a)') 'ERROR read_attr_odim:'// &
           ' more than one output value ival/rval/dval/cval/cval_varlen/ivec/rvec/dvec'// &
           ' specified, which is not allowed!'
      ierr = 1
      RETURN
    END IF

    IF (PRESENT(silent)) THEN
      issilent = silent
    ELSE
      issilent = .FALSE.
    END IF

    IF (issilent) THEN
      CALL h5eset_auto_f(0, ierror)  ! Switch off HDF5 error messages temporarily
    END IF

    ! .. Access the group/dataset/datatype of the attribute in grpname directly and return its grp_id:
    CALL h5oopen_f(file_id, TRIM(grpname), grp_id, ierror)
    IF (ierror /= 0) THEN
      IF (.NOT. issilent) THEN
        WRITE (*,'(a)') 'WARNING read_attr_odim(): group/dataset/datatype '//TRIM(grpname)// &
             ' not found in file when looking for attribute '//TRIM(attrname)
      END IF
      ierr = ierror
      CALL h5oclose_f(grp_id, ierror)
      CALL prepare_return()
      RETURN
    END IF


    CALL h5iget_type_f(grp_id, obj_typ, ierror)

!!$
!!$                                        H5I_FILE_F
!!$                                        H5I_GROUP_F
!!$                                        H5I_DATATYPE_F
!!$                                        H5I_DATASPACE_F
!!$                                        H5I_DATASET_F
!!$                                        H5I_ATTR_F
!!$                                        H5I_BADID_F

    ! .. Read the value of the attribute for the object types supported by h5oopen_f():

    IF ( ANY( obj_typ == (/H5I_GROUP_F, H5I_DATASET_F, H5I_DATATYPE_F/) ) ) THEN

      CALL h5aopen_by_name_f(grp_id, '.', TRIM(attrname), attr_id, ierror, H5P_DEFAULT_F)
      IF (ierror /= 0) THEN
        IF (.NOT. issilent) THEN
          WRITE(*,'(a)') 'WARNING read_attr_odim(): attr '//TRIM(attrname)//' not found in '// &
               TRIM(grpname)//' in file'
        END IF
        ierr = ierror
        CALL h5aclose_f(attr_id, ierror)
        CALL h5oclose_f(grp_id, ierror)
        CALL prepare_return()
        RETURN
      END IF

      CALL h5aget_type_f(attr_id, type_id, ierror)
      CALL h5tget_native_type_f(type_id, 2, h5_type, ierror)


      IF ( ANY(attr_data_type == (/4, 5, 8/)) ) THEN
        ! .. These are the array data types:
        CALL h5aget_space_f               (attr_id, dsp_id, ierror)
        CALL h5sget_simple_extent_ndims_f (dsp_id,  ndims, ierror )
        ALLOCATE (datadims   (ndims))
        ALLOCATE (maxdatadims(ndims))
        CALL h5sget_simple_extent_dims_f  (dsp_id, datadims, maxdatadims, ierror)
      ELSE
        ! .. This is a scaler:
        ndims       = 0
      END IF

!!$        H5T_NATIVE_CHAR
!!$        H5T_NATIVE_SHORT
!!$        H5T_NATIVE_INT
!!$        H5T_NATIVE_LONG
!!$        H5T_NATIVE_LLONG
!!$
!!$        H5T_NATIVE_UCHAR
!!$        H5T_NATIVE_USHORT
!!$        H5T_NATIVE_UINT
!!$        H5T_NATIVE_ULONG
!!$        H5T_NATIVE_ULLONG
!!$
!!$        H5T_NATIVE_FLOAT  /  H5T_NATIVE_REAL (fortran interface)
!!$        H5T_NATIVE_DOUBLE
!!$        H5T_NATIVE_LDOUBLE
!!$
!!$        H5T_NATIVE_B8
!!$        H5T_NATIVE_B16
!!$        H5T_NATIVE_B32
!!$        H5T_NATIVE_B64

!!$ string data type, but unclear how to query it:
!!$   H5T_FORTRAN_S1

      IF (attr_data_type == 1 .AND. ndims == 0) THEN
        CALL h5aread_f(attr_id, h5_type, iatt_value, attr_dims, ierror)
        ival = INT(iatt_value)
        ierr = ierror
      ELSE IF (attr_data_type == 2 .AND. ndims == 0) THEN
        CALL h5aread_f(attr_id, h5_type, datt_value, attr_dims, ierror)
        dval = DBLE(datt_value)
        ierr = ierror
      ELSE IF (attr_data_type == 3 .AND. ndims == 0) THEN
        catt_value(:) = REPEAT(' ',LEN(catt_value))
        CALL h5ltget_attribute_string_f(file_id, grpname, attrname, &
             catt_value, ierror)
        ierr = ierror
        ! Work-around for the case that in the hdf5-file a string attribute has not the correct length N+1 (for the add. \0)
        ! but the length N, in which case strange things happen when using h5ltget_attribute_string_f() for reading,
        ! i.e., if h5dump gives for example
!!$   DATATYPE  H5T_STRING {
!!$     STRSIZE 12;
!!$     STRPAD H5T_STR_NULLTERM;
!!$     CSET H5T_CSET_ASCII;
!!$     CTYPE H5T_C_S1;
!!$   }
!!$   DATASPACE  SCALAR
!!$   DATA {
!!$   (0): "ODIM_H5/V2_2"
!!$   }
        ! but it should be
!!$     STRSIZE 13;
        ! :
        cval(:) = REPEAT(' ',LEN(cval))
        IF (ierr == 0) THEN
          CALL h5aget_storage_size_f(attr_id, size_cattr, ierror)
          IF (size_cattr > 0) THEN
            cval = TRIM(catt_value(1:size_cattr))
          END IF
          ierr = ierror
        END IF
      ELSE IF (attr_data_type == 4 .AND. ndims == 1) THEN
        IF (ALLOCATED(ivec)) DEALLOCATE(ivec)
        ALLOCATE (ivec(datadims(1)))
        CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, ivec, datadims, ierror)
        ierr = ierror
      ELSE IF (attr_data_type == 5 .AND. ndims == 1) THEN
        IF (ALLOCATED(dvec)) DEALLOCATE(dvec)
        ALLOCATE (dvec(datadims(1)))
        CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, dvec, datadims, ierror)
        ierr = ierror
      ELSE IF (attr_data_type == 6 .AND. ndims == 0) THEN
        ! variable length string:
!!$         DATATYPE  H5T_STRING {
!!$            STRSIZE H5T_VARIABLE;
!!$            STRPAD H5T_STR_NULLTERM;
!!$            CSET H5T_CSET_ASCII;
!!$            CTYPE H5T_C_S1;
!!$         }
        ALLOCATE (rdata(LEN(cval_varlen)))
        f_ptr = C_LOC(rdata(1))
        CALL h5aread_f(attr_id, type_id, f_ptr, ierror)
        ierr = ierror
        cval_varlen(:) = ' '
        cval_varlen = C_to_F_string(rdata(1), LEN(cval_varlen))
        DEALLOCATE (rdata)
      ELSE IF (attr_data_type == 7 .AND. ndims == 0) THEN
        CALL h5aread_f(attr_id, h5_type, ratt_value, attr_dims, ierror)
        rval = DBLE(ratt_value)
        ierr = ierror
      ELSE IF (attr_data_type == 8 .AND. ndims == 1) THEN
        ALLOCATE (ratt_vec(datadims(1)))
        CALL h5aread_f(attr_id, h5_type, ratt_vec, datadims, ierror)
        ierr = ierror
        IF (ALLOCATED(rvec)) DEALLOCATE(rvec)
        ALLOCATE (rvec(datadims(1)))
        rvec = DBLE(ratt_vec)
      ELSE
        IF (.NOT. issilent) THEN
          WRITE (*,'(a,i0,a)') 'WARNING read_attr_odim(): ',  attr_data_type, &
               ' is unknown data type for attribute '//TRIM(attrname)//' in '//TRIM(grpname)
        END IF
        ierr = 4
      END IF

      IF ( ANY(attr_data_type == (/4, 5, 8/)) ) THEN
        CALL H5sclose_f(dsp_id, ierror)
      END IF
      CALL H5tclose_f(h5_type, ierror)
      CALL H5tclose_f(type_id, ierror)
      CALL h5aclose_f(attr_id, ierror)

      IF (ALLOCATED(datadims))    DEALLOCATE(datadims)
      IF (ALLOCATED(maxdatadims)) DEALLOCATE(maxdatadims)
      IF (ALLOCATED(ratt_vec))    DEALLOCATE(ratt_vec)

      IF (ierr /= 0) THEN
        IF (.NOT. issilent) THEN
          WRITE (*,'(a,i0,a,i0)') 'WARNING read_attr_odim(): Failed to read attribute '//TRIM(attrname)//' in '//TRIM(grpname)// &
               ' with ndims=', ndims, ' and attr_data_type=', attr_data_type
        END IF
      END IF

    ELSE

      IF (.NOT. issilent) THEN
        WRITE (*,'(a,i0,a)') 'WARNING read_attr_odim(): obj_typ = ', obj_typ, ' for which to read attributes is not yet implemented'
      END IF
      ierr = 5

    END IF

    ! .. Close the group/dataset/datatype of the attribute:
    CALL h5oclose_f(grp_id, ierror)

    CALL prepare_return()
    
  CONTAINS

    SUBROUTINE prepare_return
      IF (issilent) THEN
        CALL h5eset_auto_f(1, ierror)  ! Switch on HDF5 error messages again
      END IF
    END SUBROUTINE prepare_return

    FUNCTION C_to_F_string(c_string_pointer, maxlen) RESULT(f_string)

      INTEGER,                INTENT(in) :: maxlen
      TYPE(c_ptr),            INTENT(in) :: c_string_pointer
      CHARACTER(len=maxlen)              :: f_string
      CHARACTER(kind=c_char), DIMENSION(:), POINTER :: char_array_pointer => NULL()
      INTEGER :: i, length

      f_string(:) = ' '
      length = LEN(f_string)

      CALL C_F_POINTER( c_string_pointer, char_array_pointer, (/ INT(maxlen,kind=c_int) /) )
      IF (.NOT.ASSOCIATED(char_array_pointer)) THEN
        f_string = 'NULL'
        RETURN
      END IF

      DO i = 1, maxlen
        IF (char_array_pointer(i) == c_null_char .OR. i > length) THEN
          EXIT
        END IF
        f_string(i:i) = char_array_pointer(i)
      END DO

    END FUNCTION C_to_F_string
    
  END SUBROUTINE read_attr_odim

  SUBROUTINE read_dataset_odim (file_id, dsetname, ierr, &
                                ivec, dvec, iarr2d, darr2d)

    ! .. Input/Output variables:
    INTEGER(hid_t), INTENT(in)    :: file_id
    CHARACTER(len=*), INTENT(in)  :: dsetname    ! full path and dataset name within group
    INTEGER,          INTENT(out) :: ierr
    INTEGER,          INTENT(out), ALLOCATABLE, OPTIONAL :: ivec(:)
    REAL(KIND=dp), INTENT(out), ALLOCATABLE, OPTIONAL :: dvec(:)
    INTEGER,          INTENT(out), ALLOCATABLE, OPTIONAL :: iarr2d(:,:)
    REAL(KIND=dp), INTENT(out), ALLOCATABLE, OPTIONAL :: darr2d(:,:)

    ! .. Local variables:
    INTEGER            :: ierror, dset_data_type, ntrue, ndims
    INTEGER(hid_t)     :: dsp_id, dset_id, type_id, h5_type
    INTEGER            :: obj_typ
    INTEGER(hsize_t), ALLOCATABLE :: datadims(:), maxdatadims(:)

    ierr = 0

    ntrue = 0
    IF (PRESENT(ivec)) THEN
      ntrue = ntrue + 1
      dset_data_type = 1
    END IF
    IF (PRESENT(dvec)) THEN
      ntrue = ntrue + 1
      dset_data_type = 2
    END IF
    IF (PRESENT(iarr2d)) THEN
      ntrue = ntrue + 1
      dset_data_type = 3
    END IF
    IF (PRESENT(darr2d)) THEN
      ntrue = ntrue + 1
      dset_data_type = 4
    END IF
    IF (ntrue > 1) THEN
      ! More than one output variable present, but this is not allowed.
      !  read_dataset_odim() should be called for one specific output datum and data type only!
      WRITE(*,'(a)') 'WARNING read_dataset_odim():'// &
                     ' more than one output value ivec/dvec/iarr2d/darr2d specified, which is not allowed!'
      ierr = 1
      RETURN
    END IF


    ! .. Access the dataset of the attribute in grpname directly and return its grp_id:
    CALL h5oopen_f(file_id, TRIM(dsetname), dset_id, ierror)
    IF (ierror /= 0) THEN
      WRITE(*,'(a)') 'WARNING read_dataset_odim(): group/dataset/datatype '//TRIM(dsetname)// &
                     ' not found in file'
      ierr = ierror
      CALL h5oclose_f (dset_id, ierror)
      RETURN
    END IF


    CALL h5iget_type_f(dset_id, obj_typ, ierror)  ! typ must be = H5I_DATASET_F

!!$
!!$                                        H5I_FILE_F
!!$                                        H5I_GROUP_F
!!$                                        H5I_DATATYPE_F
!!$                                        H5I_DATASPACE_F
!!$                                        H5I_DATASET_F
!!$                                        H5I_ATTR_F
!!$                                        H5I_BADID_F


    ! .. Read the value of the attribute for the object types supported by h5oopen_f():

    IF ( obj_typ == H5I_DATASET_F ) THEN

      CALL h5dget_type_f        (dset_id, type_id, ierror)
      CALL h5tget_native_type_f (type_id, 2, h5_type, ierror)
      CALL h5dget_space_f       (dset_id, dsp_id, ierror)
      CALL h5sget_simple_extent_ndims_f (dsp_id, ndims, ierror)
      ALLOCATE(datadims(ndims))
      ALLOCATE(maxdatadims(ndims))
      CALL h5sget_simple_extent_dims_f (dsp_id, datadims, maxdatadims, ierror)

      IF (dset_data_type == 1 .AND. ndims == 1) THEN

        IF (ALLOCATED(ivec)) DEALLOCATE(ivec)
        ALLOCATE( ivec(datadims(1)) )
        CALL h5ltread_dataset_f (dset_id, '.', H5T_NATIVE_INTEGER, ivec, datadims, ierror)
        ierr = ierror

      ELSE IF (dset_data_type == 2 .AND. ndims == 1) THEN

        IF (ALLOCATED(dvec)) DEALLOCATE(dvec)
        ALLOCATE( dvec(datadims(1)) )
        CALL h5ltread_dataset_f (dset_id, '.', H5T_NATIVE_DOUBLE, dvec, datadims, ierror)
        ierr = ierror

      ELSE IF (dset_data_type == 3 .AND. ndims == 2) THEN

        IF (ALLOCATED(iarr2d)) DEALLOCATE(iarr2d)
        ALLOCATE( iarr2d(datadims(1),datadims(2)) )
        CALL h5ltread_dataset_f (dset_id, '.', H5T_NATIVE_INTEGER, iarr2d, datadims, ierror)
        ierr = ierror

      ELSE IF (dset_data_type == 4 .AND. ndims == 2) THEN

        IF (ALLOCATED(darr2d)) DEALLOCATE(darr2d)
        ALLOCATE( darr2d(datadims(1),datadims(2)) )
        CALL h5ltread_dataset_f (dset_id, '.', H5T_NATIVE_DOUBLE, darr2d, datadims, ierror)
        ierr = ierror
        
      ELSE
        WRITE (*,'(a,i0,a)') 'WARNING read_dataset_odim(): ', dset_data_type, &
             ' is unknown data type/dimension for dataset '//TRIM(dsetname)
        ierr = 4
      END IF

      IF (ierr /= 0) THEN
        WRITE (*,'(a,i0,a,i0)') 'WARNING read_dataset_odim(): Failed to read dataset '//TRIM(dsetname)// &
             ' with ndims=', ndims, ' and dset_data_type=', dset_data_type
      END IF

      CALL h5sclose_f (dsp_id,  ierror)
      CALL h5tclose_f (h5_type, ierror)
      CALL h5tclose_f (type_id, ierror)

    ELSE

      WRITE (*,'()') 'WARNING read_dataset_odim(): '//TRIM(dsetname)//' is not of type dataset!'
      ierr = 5

    END IF

    ! .. Close the dataset of the attribute:
    CALL h5oclose_f (dset_id, ierror)

  END SUBROUTINE read_dataset_odim

  SUBROUTINE get_stationID_from_HDF5_source(srcstring, station_id, country_id, statid_char, ierr)

    CHARACTER(len=*), INTENT(in)   :: srcstring
    INTEGER,          INTENT(out)  :: station_id, country_id, ierr
    CHARACTER(len=5), INTENT(out)  :: statid_char

    INTEGER                        :: nwords, nsubwords, i
    CHARACTER(len=20), ALLOCATABLE :: words(:), subwords(:)    ! for results of subroutine split_string()
    CHARACTER(len=20)              :: wmoword, plcword, nodword
    LOGICAL                        :: have_wmo, have_plc, have_nod

    ierr = 0
    station_id = 999999
    country_id = 99
    statid_char = 'XXXXX'
    have_wmo = .FALSE.
    have_nod = .FALSE.
    have_plc = .FALSE.

    CALL split_string(TRIM(srcstring), ',', LEN(words), words, nwords)
    IF (nwords <= 0) THEN
      WRITE (*,'(a)') 'WARNING get_stationID_from_HDF5_source(): Failed to split source string '// &
                      TRIM(srcstring)
      ierr = 1
    END IF

    DO i=1, nwords
      CALL split_string(TRIM(words(i)), ':', LEN(subwords(i)), subwords, nsubwords)
      IF (nsubwords /= 2) THEN
        ierr = 2
      ELSE
        SELECT CASE (TRIM(subwords(1)))
        CASE ('WMO')
          wmoword = subwords(2)
          have_wmo = .TRUE.
        CASE ('NOD')
          nodword = subwords(2)
          have_nod = .TRUE.
        CASE ('PLC')
          plcword = subwords(2)
          have_plc = .TRUE.
        END SELECT
      END IF
      IF (ALLOCATED(subwords)) DEALLOCATE(subwords)
    END DO

    IF (ALLOCATED(words)) DEALLOCATE (words)

    IF (ierr == 2) THEN
      IF (have_wmo .OR. have_nod .OR. have_plc) THEN
        ! The source string contains additional komma-separated words with wrong content, but at least one of the
        ! station identifiers has been found. We can work with this information and reset the error status to 0.
        ierr = 0
      ELSE
        ! None of the station identifiers could be found in srcstring. This is problematic and we issue a WARNING
        ! and do not reset the error status.
        WRITE(*,'(a)') 'WARNING get_stationID_from_HDF5_source():'//&
                       ' Failed to get station identifier from source string '//TRIM(srcstring)
      END IF
    END IF
        
    IF (have_wmo) THEN
      READ (wmoword(1:5), *) station_id    ! full WMO ID
      READ (wmoword(1:2), *) country_id    ! country ID
      IF (have_nod) THEN
        statid_char = TRIM(nodword)
      END IF
    ELSE
      ! Some countries do not or not always deliver WMO IDs but only NODs. This can be repaired here:
      IF (have_nod) THEN
        statid_char = TRIM(nodword)
        SELECT CASE (nodword(1:2))
        CASE ('it')
          ! Try to guess station_id and country_id from the station short name:
          SELECT CASE (TRIM(nodword(3:5)))
          CASE ('spc')
            ! Italy, San Pietro Capofiume:
            station_id = 16144
            country_id = 16
          CASE ('gat')
            ! Italy, Gattatico (no official WMO station yet):
            ! THOMAS modification: Modifica probabilmenteirrilevante, serve per determinare l'ID 
            ! del radar e la nazione di  appartenenza a partire da "plcword" (station short name)
            station_id = 16199
            country_id = 16
          CASE default
            WRITE(*,'(a)') 'WARNING get_stationID_from_HDF5_source():'// &
                           ' could not determine IT station_id from name '//TRIM(plcword)
            ierr = 41
          END SELECT
        CASE ('fr')
          ! France
        CASE ('fi')
          ! Finnland
        CASE ('pl')
          ! Poland
        CASE ('cz')
          ! Czech Republic
        CASE ('sk')
          ! Slovakia
        CASE ('dk')
          ! Denmark
          SELECT CASE (TRIM(nodword(3:5)))
          CASE ('sin')
            station_id = 6034
            country_id = 6
          CASE ('rom')
            station_id = 6096
            country_id = 6
          CASE ('vir')
            station_id = 6103
            country_id = 6
          CASE ('ste')
            station_id = 6173
            country_id = 6
          CASE ('bor')
            station_id = 6194
            country_id = 6
          CASE default
            WRITE(*,'(a)') 'WARNING get_stationID_from_HDF5_source():'// &
                           ' could not determine DK station_id from name '//TRIM(plcword)
            ierr = 41
          END SELECT
        CASE ('nl')
          ! Netherlands
          SELECT CASE (nodword(3:5))
          CASE ('dhl')
            station_id = 6234
            country_id = 6
          CASE ('hrw')
            station_id = 6356
            country_id = 6
          END SELECT
        CASE ('gr')
          ! Greece
        CASE ('rs')
          ! Serbia
        END SELECT
      ELSE IF (have_plc) THEN
        ! Try to guess station_id and country_id from the station short name:
        SELECT CASE (TRIM(plcword))
        CASE ('itspc')
          ! Italy, San Pietro Capofiume:
          station_id = 16144
          country_id = 16
        CASE ('itgat')
          ! Italy, Gattatico (no official WMO station yet):
          ! THOMAS modification: Modifica probabilmenteirrilevante, serve per determinare l'ID 
          ! del radar e la nazione di  appartenenza a partire da "plcword" (station short name)
          station_id = 16199
          country_id = 16
        CASE default
          WRITE(*,'(a)') 'WARNING get_stationID_from_HDF5_source():'// &
                         ' could not determine plc station_id from name '//TRIM(plcword)
          ierr = 4
        END SELECT
      ELSE
        WRITE(*,'(a)') 'WARNING get_stationID_from_HDF5_source():'// &
                       ' could not determine non-wmo/non-ndo station_id from source '//TRIM(srcstring)
        ierr = 3
      END IF
    END IF

  END SUBROUTINE get_stationID_from_HDF5_source

#endif    /* HDF5_RADAR_INPUT  */

  !==============================================================================

  SUBROUTINE append_obstimes_to_rsm ( rsm, cfileident, fids, zydate_ini_mod, cdate, err )

    !------------------------------------------------------------------------------
    !
    ! Description: Determines the different observation times given in the vector cdate(:)
    !              and compute internal model time in seconds since model start
    !              and store them in rsm.
    !
    !              It is assumed that cdate(:) is sorted according to blocks of
    !               same obstimes in strict ascending order.
    !
    !              If rsm does not contain any obstimes yet, the new obs times
    !               from cdate(:) are put to rsm%obs_times_obs(:), otherwise they are
    !               appended to the existing list.
    !              Likewise, the start- and end-indices of the blocks of same obstime
    !               are put to rsm%obs_startrec(:) and rsm%obs_endrec(:).
    !
    !              cfileident is a string to identify the file from which
    !               manel was read, and is used for error messages.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! In, InOut, Outgoing:

    TYPE(radar_meta_type),   INTENT(inout) :: rsm
    CHARACTER(LEN=*),        INTENT(in)    :: cfileident
    INTEGER,                 INTENT(in)    :: fids(:)   ! file ident nos. for rsm%obs_startrec(:,fid)
    CHARACTER(LEN=14),       INTENT(in)    :: zydate_ini_mod
    CHARACTER(LEN=14),       INTENT(in)    :: cdate(:)
    INTEGER,                 INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    INTEGER                                :: i, n, isecdif, ireport, iu, io
    REAL(KIND=dp),           ALLOCATABLE   :: sec_ob(:)
    INTEGER,                 ALLOCATABLE   :: records(:), endrecs(:)
    CHARACTER(LEN=14)                      :: otime, timecheck
    CHARACTER(LEN=14),       ALLOCATABLE   :: cdate_loc(:)
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'append_obstimes_to_rsm'
    CHARACTER(LEN=80)                      :: yzerrmsg

    INTEGER                                ::  &
         ierrf(2,6)          ! error status for time difference routine

    err = 0

    n = SIZE(cdate)
    IF (n <= 0) THEN
      err = 9
      WRITE (*,'(a,i6.6)') &
           TRIM(yzroutine)//' WARNING: no obs times found in NetCDF file ' // &
           TRIM(cfileident) // ' for station ID ', rsm%station_id
      RETURN
    END IF

    ALLOCATE( records(n), endrecs(n), sec_ob(n), cdate_loc(n) )

    ireport    = 1
    CALL diff_seconds ( zydate_ini_mod, cdate(1), isecdif, ierrf )
    records(1) = 1
    endrecs(1) = 1
    otime      = cdate(1)
    cdate_loc(1) = otime
    sec_ob(1)  = isecdif ! Difference to model start time

    ! Everytime that cdate changes, a new obs time is recorded:
    DO i = 2, n

      timecheck  = cdate(i)
      IF ( timecheck /= otime ) THEN
        ireport               = ireport + 1
        cdate_loc(ireport)    = timecheck
        CALL diff_seconds ( zydate_ini_mod, timecheck, isecdif, ierrf )
        sec_ob(ireport)       = isecdif ! Difference to model start time
        otime                 = timecheck
        records(ireport)      = i
        endrecs(ireport)      = i
      ELSE
        endrecs(ireport)      = i
      END IF

    END DO

    DO i = 2, ireport
      ! Check that different obs_times_obs are in monotonic ascending order.
      !  Otherwise, the data are most likely corrupt:
      IF (sec_ob(i) < sec_ob(i-1)) THEN
        err = 10
        WRITE (*,'(a,i6.6,a)') &
             TRIM(yzroutine)//' ERROR: obs times not in ascending order in NetCDF file ' // &
             TRIM(cfileident) // ' for station ID ', rsm%station_id, '! Perhaps corrupt data!'
        RETURN
      END IF
    END DO
        
    ! If the meta data do not contain any obstimes yet, start from index 0,
    !  otherwise append to the existing list:
    IF (rsm%nobs_times_obs < 0) rsm%nobs_times_obs = 0
    iu                        = rsm%nobs_times_obs + 1
    io                        = rsm%nobs_times_obs + ireport
    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        rsm%obs_startrec(iu:io,fids(i)) = records(1:ireport)
        rsm%obs_endrec(iu:io,fids(i))   = endrecs(1:ireport)
      END IF
    END DO
    rsm%obs_times_obs(iu:io)    = REAL(sec_ob(1:ireport), KIND=dp)
    rsm%obs_cdate(iu:io)        = cdate_loc(1:ireport)
    rsm%nobs_times_obs          = rsm%nobs_times_obs + ireport

    DEALLOCATE( records, endrecs, sec_ob, cdate_loc )

  END SUBROUTINE append_obstimes_to_rsm

  !==============================================================================
  !==============================================================================  

  SUBROUTINE reconstruct_file_series_mch (caller, infilename_series, infilenames, ninfiles, err)

    !------------------------------------------------------------------------------
    !
    ! Re-construct the filename series of the PPI-files of MeteoSwiss
    !  from the given input filename,
    !  which contains the basic filename, the number of PPI sweeps and the list
    !  of elevation identifiers.
    !
    !------------------------------------------------------------------------------
    
    CHARACTER(len=*),        INTENT(in)                 :: caller
    CHARACTER(len=*),        INTENT(in)                 :: infilename_series
    CHARACTER(len=cobsflen), INTENT(inout), ALLOCATABLE :: infilenames(:)
    INTEGER,                 INTENT(out)                :: ninfiles, err
    
    INTEGER                        :: i, k, nwords
    CHARACTER(len=20), ALLOCATABLE :: words(:)    ! for results of subroutine split_string()

    err = 0

    READ (infilename_series(16:18), '(i3.3)') ninfiles
    k = INDEX(infilename_series,'|') - 1   ! length of the actual filename (anything before the first '|'
    CALL split_string (TRIM(infilename_series), '|', LEN(words), words, nwords)
    IF (nwords /= ninfiles+1) THEN
      WRITE(*,*) 'WARNING '//TRIM(caller)//': Inconsistent infilename '//TRIM(infilename_series)
      IF (ALLOCATED(words)) DEALLOCATE(words)
      err = 1
      RETURN
    END IF
    ALLOCATE(infilenames(ninfiles))
    DO i = 1, ninfiles
      infilenames(i) = infilename_series(1:k)
      WRITE(infilenames(i)(16:18), '(a3)') words(i+1)(1:3)
    END DO

    IF (ALLOCATED(words)) DEALLOCATE(words)

  END SUBROUTINE reconstruct_file_series_mch

  !==============================================================================
  !==============================================================================  
  
  SUBROUTINE reconstruct_file_series_opera (caller, infilename_series, infilenames, ninfiles, err)

    !------------------------------------------------------------------------------
    !
    ! Re-construct the filename series of the PPI-files of the OPERA data hub (T_PA*)
    !  from the given input filename,
    !  which contains the basic filename, the number of PPI sweeps and the list
    !  of elevation identifiers.
    !
    !------------------------------------------------------------------------------
    
    CHARACTER(len=*),        INTENT(in)                 :: caller
    CHARACTER(len=*),        INTENT(in)                 :: infilename_series
    CHARACTER(len=cobsflen), INTENT(inout), ALLOCATABLE :: infilenames(:)
    INTEGER,                 INTENT(out)                :: ninfiles, err
    
    INTEGER                        :: i, k, nwords
    CHARACTER(len=20), ALLOCATABLE :: words(:)    ! for results of subroutine split_string()

    err = 0

    CALL split_string (TRIM(infilename_series), '|', LEN(words), words, nwords)
    IF (nwords > 1) THEN
      ! There is more than 1 file in the series:
      !  e.g. T_PAZF40_C_LFPW_20220104104922.h5|A-4522|B-4602|C-4652|D-4742|E-4832|F-4922
      ninfiles = IACHAR(infilename_series(6:6)) - IACHAR('A') + 1
      k = INDEX(infilename_series,'|') - 1   ! length of the actual filename (anything before the first '|'
      IF (nwords /= ninfiles+1) THEN
        WRITE(*,'(a)') 'WARNING '//TRIM(caller)//': Inconsistent infilename '//TRIM(infilename_series)
        IF (ALLOCATED(words)) DEALLOCATE(words)
        err = 1
        RETURN
      END IF
      ALLOCATE(infilenames(ninfiles))
      infilenames(:)(:) = ' '
      DO i = 1, ninfiles
        infilenames(i) = infilename_series(1:k)
        WRITE(infilenames(i)(6:6), '(a1)') words(i+1)(1:1)  ! The elevation letter
        WRITE(infilenames(i)(27:30), '(a)') words(i+1)(3:6) ! The minutes and seconds
      END DO
    ELSE
      ! There is only one filename in the series:
      ninfiles = 1
      ALLOCATE(infilenames(ninfiles))
      infilenames(:)(:) = ' '
      infilenames(1) = TRIM(infilename_series)
    END IF

    IF (ALLOCATED(words)) DEALLOCATE(words)

  END SUBROUTINE reconstruct_file_series_opera

  !==============================================================================
  !==============================================================================  
  
  SUBROUTINE put_obstimes_to_rsm_onetime ( rsm, cfileident, fids, zydate_ini_mod, cdate, err )

    !------------------------------------------------------------------------------
    !
    ! Description: Determines the observation time given in the string cdate
    !              and compute internal model time in seconds since model start
    !              and store them in rsm.
    !
    !              It is assumed that cdate(:) is sorted according to blocks of
    !               same obstimes in strict ascending order.
    !
    !              If rsm does not contain any obstimes yet, the new obs times
    !               from cdate(:) are put to rsm%obs_times_obs(:), otherwise they are
    !               appended to the existing list.
    !              Likewise, the start- and end-indices of the blocks of same obstime
    !               are put to rsm%obs_startrec(:) and rsm%obs_endrec(:).
    !
    !              cfileident is a string to identify the file from which
    !               manel was read, and is used for error messages.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! In, InOut, Outgoing:

    TYPE(radar_meta_type_onetime),   INTENT(inout) :: rsm
    CHARACTER(LEN=*),        INTENT(in)    :: cfileident
    INTEGER,                 INTENT(in)    :: fids(:)   ! file ident nos. for rsm%obs_startrec(:,fid)
    CHARACTER(LEN=14),       INTENT(in)    :: zydate_ini_mod
    CHARACTER(LEN=14),       INTENT(in)    :: cdate
    INTEGER,                 INTENT(out)   :: err

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    INTEGER                     :: i, n, isecdif
    REAL(KIND=dp)               :: sec_ob
    CHARACTER(LEN=*), PARAMETER :: yzroutine = 'put_obstimes_to_rsm_onetime'
    CHARACTER(LEN=80)           :: yzerrmsg

    INTEGER                     ::  &
         ierrf(2,6)          ! error status for time difference routine

    err = 0

    IF ( LEN_TRIM(cdate) == 0 .OR. TRIM(cdate) == 'YYYYMMDDhhmmss' ) THEN
      err = 9
      WRITE (*,'(a,i6.6)') &
           TRIM(yzroutine)//' WARNING: no obs time found in file ' // &
           TRIM(cfileident) // ' for station ID ', rsm%station_id
      rsm%nobs_times_obs = 0
      RETURN
    END IF

    CALL diff_seconds ( zydate_ini_mod, cdate, isecdif, ierrf )
    sec_ob     = isecdif ! Difference to model start time


    DO i=1, SIZE(fids)
      IF (fids(i) > 0) THEN
        rsm%obs_startrec(1:1,fids(i)) = 1
        rsm%obs_endrec(1:1,fids(i))   = 1
      END IF
    END DO
    rsm%obs_times_obs(1:1)    = REAL(sec_ob, KIND=dp)
    rsm%obs_cdate(1:1)        = cdate
    rsm%nobs_times_obs        = 1    

  END SUBROUTINE put_obstimes_to_rsm_onetime

  !==============================================================================

  ! merge rsm_in to the existing rsm_out and, if lcheck_inputrecords = .true.,
  !  compare the records for same obs_times_obs. rsm_in comes from reading of a
  !  data file which may contain more than one datakind (vr, z, qv, qz).
  ! The routine merges the names of the data files as well as the obs_times_obs
  !  into one meta data structure.
  ! Some cross-checks of the consistency of the obs_times_obs and startrec/endrec's between the different
  !  datakinds are preformed.

  ! This is the version if both rsm_in and rsm_out can have multiple observation times:
  SUBROUTINE merge_metadata_m( rsm_in, rsm_out, err )

    IMPLICIT NONE

    TYPE(radar_meta_type), INTENT(in)    :: rsm_in
    TYPE(radar_meta_type), INTENT(inout) :: rsm_out
    INTEGER,               INTENT(out)   :: err

    TYPE(radar_meta_type)    :: rsm_tmp
    INTEGER                  :: i, ii, m, n, itime
    INTEGER, DIMENSION(SIZE(rsm_in%obsfile(1,:))) :: fid_in, fid_last

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'merge_metadata_m'
    CHARACTER(LEN=80)        :: yzerrmsg

    err = 0
    m = rsm_out%nobs_times_obs
    fid_last(:) = missval_int
    n = rsm_in%nobs_times_obs
    merge_loop: DO ii = 1, n

      fid_in(:) = missval_int
      DO i=1, SIZE(rsm_in%obsfile(ii,:))
        IF (TRIM(rsm_in%obsfile(ii,i)) /= obsfile_missingname) THEN
          fid_in(i) = i
          fid_last(i) = i
!!$ UB          EXIT
        END IF
      END DO

      itime = missval_int
      DO i=1, m
        IF ( ABS(rsm_out%obs_times_obs(i)-rsm_in%obs_times_obs(ii)) <= 1e-12_dp ) THEN
          itime = i
          EXIT
        END IF
      END DO

      IF (itime < 0 .AND. ANY(fid_in(:) > 0)) THEN
        ! This is a new obs_time, so set itime to m to append it to the list below:
        m = m + 1
        IF (m <= nobstimes_max) THEN
          itime = m
          rsm_out%nobs_times_obs = m
        ELSE
          err = 3
          WRITE (*,'(a,i5,a,i6.6,a)') &
            'WARNING '//TRIM(yzroutine)//': too many obs_times_obs (', m,') for radar '// &
            TRIM(rsm_out%station_name)//' ',rsm_out%station_id,' '//TRIM(rsm_out%scanname)// &
            '! Increase nobstimes_max!'
          EXIT merge_loop
        ENDIF
      END IF

      IF (itime > 0 .AND. ANY(fid_in(:) > 0)) THEN
        ! Store records and filenames at the appropriate time index itime
        !  (either a new obs_time or an existing one):
        DO i=1, SIZE(rsm_in%obsfile(ii,:))
          IF (TRIM(rsm_in%obsfile(ii,i)) /= obsfile_missingname) THEN
            ! Check, whether we already have identified matching obs data for time itime and radar
            !   parameter. If so (else branch), we have encountered ambiguous data and will mark it
            !   as such (using lobs_avail for tracking). After merging is completed (outside of
            !   merge_metadata), these data will be set to be skipped.
            ! JM: Argh, doesn't work properly for the multi-time data here. Hence switched off here
            !   again. So far, we don't have data that required this check (ie polarimetric
            !   moments) in multi-time format, so it's not a problem yet. Once, we do, we have to
            !   revisit this...
            !IF (TRIM(rsm_out%obsfile(itime,i)) == obsfile_missingname) THEN
              rsm_out%obsfile(itime,i)      = rsm_in%obsfile(ii,i)
              rsm_out%naz_ncdf(itime,i)     = rsm_in%naz_ncdf(ii,i)
              rsm_out%obs_startrec(itime,i) = rsm_in%obs_startrec(ii,i)
              rsm_out%obs_endrec  (itime,i) = rsm_in%obs_endrec  (ii,i)
            !ELSE
            !  rsm_out%lobs_avail(itime,i)    = .TRUE.
            !  err = 5
            !  WRITE(*,'(a,i2,a,i6.6,a)') &
            !    'WARNING '//TRIM(yzroutine)//': multiple matching obs files for a parameter (fid #',i,&
            !    ') for radar '// &
            !    TRIM(rsm_out%station_name)//' ',rsm_out%station_id,' '//TRIM(rsm_out%scanname)//':'
            !  WRITE(*,'(a)') '  '//TRIM(rsm_out%obsfile(itime,i))
            !  WRITE(*,'(a)') '  '//TRIM(rsm_in%obsfile(ii,i))
            !END IF
          END IF
        END DO
        rsm_out%obsfile_format(itime)     = rsm_in%obsfile_format(ii)
        rsm_out%obs_times_obs(itime)      = rsm_in%obs_times_obs(ii)
        rsm_out%obs_cdate(itime)          = rsm_in%obs_cdate(ii)
        DO i=1, rsm_out%nel
          IF (rsm_out%ext_nyq(i,itime) < miss_threshold .AND. rsm_in%ext_nyq(i,ii) >= miss_threshold) THEN
            rsm_out%ext_nyq(i,itime) = rsm_in%ext_nyq(i,ii)
          END IF
        END DO
      END IF

    END DO merge_loop

    IF (ANY(fid_last(:) > 0)) THEN
      DO i=1, SIZE(rsm_in%obsfile(ii,:))
        IF (fid_last(i) > 0) THEN
          IF (rsm_out%nrep_ncdf(fid_last(i)) < 0) THEN
            ! .. nrep_ncdf for the datakind fid_last has not yet been set for the first time:
            rsm_out%nrep_ncdf(fid_last(i)) = rsm_in%nrep_ncdf(fid_last(i))
          ELSE
            ! .. otherwise, increment nrep_ncdf the datakind fid_last:
            rsm_out%nrep_ncdf(fid_last(i)) = rsm_out%nrep_ncdf(fid_last(i)) + MAX(rsm_in%nrep_ncdf(fid_last(i)), 0)
          END IF
        END IF
      END DO
    END IF

    !============================================================
    ! Merge other scalar parameters that are not in every file:

! ra_start has the wrong value of 0.0 in the data files, so do not overtake it!
!   value is taken from the background list!
!!$    IF (rsm_out%ra_start <= 0.0_dp .AND. rsm_in%ra_start > 0.0_dp) THEN
!!$      rsm_out%ra_start = rsm_in%ra_start
!!$    END IF
! az_start is not necessarily the center of the first output azimut sampling interval, so do not overtake it:
!   value is taken from the background list!
!!$    IF (rsm_out%az_start <= -1.0_dp .AND. rsm_in%az_start > -1.0_dp) THEN
!!$     rsm_out%az_start = rsm_in%az_start
!!$    END IF

    DO i=1, rsm_out%nel
      IF (rsm_out%high_nyq(i) < miss_threshold .AND. rsm_in%high_nyq(i) >= miss_threshold) THEN
        rsm_out%high_nyq(i) = rsm_in%high_nyq(i)
      END IF
      IF (rsm_out%dualprf_ratio(i) < miss_threshold .AND. rsm_in%dualprf_ratio(i) >= miss_threshold) THEN
        rsm_out%dualprf_ratio(i) = rsm_in%dualprf_ratio(i)
      END IF
      IF (rsm_out%prf(i) < miss_threshold .AND. rsm_in%prf(i) >= miss_threshold) THEN
        rsm_out%prf(i) = rsm_in%prf(i)
      END IF
    END DO
    IF (rsm_out%rngate_len <= 0.0_dp .AND. rsm_in%rngate_len > 0.0_dp) THEN
      rsm_out%rngate_len = rsm_in%rngate_len
    END IF
    IF (rsm_out%num_gates <= 0 .AND. rsm_in%num_gates > 0) THEN
      rsm_out%num_gates = rsm_in%num_gates
    END IF
    IF (rsm_out%num_pulses  <= 0 .AND. rsm_in%num_pulses  > 0) THEN
      rsm_out%num_pulses  = rsm_in%num_pulses
    END IF

  END SUBROUTINE merge_metadata_m

  ! This is the version if rsm_in has only one observation time and rsm_out has multiple observation times:
  SUBROUTINE merge_metadata_o( rsm_in, rsm_out, err )

    IMPLICIT NONE

    TYPE(radar_meta_type_onetime), INTENT(in)    :: rsm_in
    TYPE(radar_meta_type),         INTENT(inout) :: rsm_out
    INTEGER,                       INTENT(out)   :: err

    TYPE(radar_meta_type)    :: rsm_tmp
    INTEGER                  :: i, m, itime
    INTEGER                  :: fid_in(SIZE(rsm_in%obsfile(1,:)))

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'merge_metadata_o'
    CHARACTER(LEN=80)        :: yzerrmsg

    err = 0
    m = rsm_out%nobs_times_obs

    IF (rsm_in%nobs_times_obs == 1) THEN
     
      fid_in(:) = missval_int
      DO i=1, SIZE(rsm_in%obsfile(1,:))
        IF (TRIM(rsm_in%obsfile(1,i)) /= obsfile_missingname) THEN
          fid_in(i) = i
        END IF
      END DO

      itime = missval_int
      DO i=1, m
        IF ( ABS(rsm_out%obs_times_obs(i)-rsm_in%obs_times_obs(1)) <= 1e-12_dp ) THEN
          itime = i
          EXIT
        END IF
      END DO

      IF (itime < 0 .AND. ANY(fid_in(:) > 0)) THEN
        ! This is a new obs_time, so set itime to m to append it to the list below:
        m = m + 1
        IF (m <= nobstimes_max) THEN
          itime = m
          rsm_out%nobs_times_obs = m
        ELSE
          err = 3
          WRITE (*,'(a,i5,a,i6.6,a)') &
            'WARNING '//TRIM(yzroutine)//': too many obs_times_obs (', m,') for radar '// &
            TRIM(rsm_out%station_name)//' ',rsm_out%station_id,' '//TRIM(rsm_out%scanname)// &
            '! Increase nobstimes_max!'
        ENDIF
      END IF

      IF (itime > 0 .AND. ANY(fid_in(:) > 0) .AND. err == 0) THEN
        ! Store records and filenames at the appropriate time index itime
        !  (either a new obs_time or an existing one):
        DO i=1, SIZE(rsm_in%obsfile(1,:))
          IF (TRIM(rsm_in%obsfile(1,i)) /= obsfile_missingname) THEN
            ! Check, whether we already have identified matching obs data for time itime and radar
            !   parameter. If so (else branch), we have encountered ambiguous data and will mark it
            !   as such (using lobs_avail for tracking). After merging is completed (outside of
            !   merge_metadata), these data will be set to be skipped.
            IF (TRIM(rsm_out%obsfile(itime,i)) == obsfile_missingname) THEN
              rsm_out%obsfile(itime,i)      = rsm_in%obsfile(1,i)
              rsm_out%naz_ncdf(itime,i)     = rsm_in%naz_ncdf(1,i)
              rsm_out%obs_startrec(itime,i) = rsm_in%obs_startrec(1,i)
              rsm_out%obs_endrec  (itime,i) = rsm_in%obs_endrec  (1,i)
            ELSE
              rsm_out%lobs_avail(itime,i)    = .TRUE.
              err = 5
              WRITE(*,'(a,i2,a,i6.6,a)') &
                'WARNING '//TRIM(yzroutine)//': multiple matching obs files for a parameter (fid #',i,&
                ') for radar '// &
                TRIM(rsm_out%station_name)//' ',rsm_out%station_id,' '//TRIM(rsm_out%scanname)//':'
              WRITE(*,'(a)') '  '//TRIM(rsm_out%obsfile(itime,i))
              WRITE(*,'(a)') '  '//TRIM(rsm_in%obsfile(1,i))
              !RETURN
            END IF
          END IF
        END DO
        rsm_out%obsfile_format(itime)     = rsm_in%obsfile_format(1)
        rsm_out%obs_times_obs(itime)      = rsm_in%obs_times_obs(1)
        rsm_out%obs_cdate(itime)          = rsm_in%obs_cdate(1)
        DO i=1, rsm_out%nel
          IF (rsm_out%ext_nyq(i,itime) < miss_threshold .AND. rsm_in%ext_nyq(i,1) >= miss_threshold) THEN
            rsm_out%ext_nyq(i,itime) = rsm_in%ext_nyq(i,1)
          END IF
        END DO
      END IF

    END IF

    IF (ANY(fid_in(:) > 0)) THEN
      DO i=1, SIZE(rsm_in%obsfile(1,:))
        IF (fid_in(i) > 0) THEN
          IF (rsm_out%nrep_ncdf(fid_in(i)) < 0) THEN
            ! .. nrep_ncdf for the datakind fid_last has not yet been set for the first time:
            rsm_out%nrep_ncdf(fid_in(i)) = rsm_in%nrep_ncdf(fid_in(i))
          ELSE
            ! .. otherwise, increment nrep_ncdf the datakind fid_last:
            rsm_out%nrep_ncdf(fid_in(i)) = rsm_out%nrep_ncdf(fid_in(i)) + MAX(rsm_in%nrep_ncdf(fid_in(i)), 0)
          END IF
        END IF
      END DO
    END IF

    !============================================================
    ! Merge other scalar parameters that are not in every file:

! ra_start has the wrong value of 0.0 in the data files, so do not overtake it!
!   value is taken from the background list!
!!$    IF (rsm_out%ra_start <= 0.0_dp .AND. rsm_in%ra_start > 0.0_dp) THEN
!!$      rsm_out%ra_start = rsm_in%ra_start
!!$    END IF
! az_start is not necessarily the center of the first output azimut sampling interval, so do not overtake it:
!   value is taken from the background list!
!!$    IF (rsm_out%az_start <= -1.0_dp .AND. rsm_in%az_start > -1.0_dp) THEN
!!$     rsm_out%az_start = rsm_in%az_start
!!$    END IF

    DO i=1, rsm_out%nel
      IF (rsm_out%high_nyq(i) < miss_threshold .AND. rsm_in%high_nyq(i) >= miss_threshold) THEN
        rsm_out%high_nyq(i) = rsm_in%high_nyq(i)
      END IF
      IF (rsm_out%dualprf_ratio(i) < miss_threshold .AND. rsm_in%dualprf_ratio(i) >= miss_threshold) THEN
        rsm_out%dualprf_ratio(i) = rsm_in%dualprf_ratio(i)
      END IF
      IF (rsm_out%prf(i) < miss_threshold .AND. rsm_in%prf(i) >= miss_threshold) THEN
        rsm_out%prf(i) = rsm_in%prf(i)
      END IF
    END DO
    IF (rsm_out%rngate_len <= 0.0_dp .AND. rsm_in%rngate_len > 0.0_dp) THEN
      rsm_out%rngate_len = rsm_in%rngate_len
    END IF
    IF (rsm_out%num_gates <= 0 .AND. rsm_in%num_gates > 0) THEN
      rsm_out%num_gates = rsm_in%num_gates
    END IF
    IF (rsm_out%num_pulses  <= 0 .AND. rsm_in%num_pulses  > 0) THEN
      rsm_out%num_pulses  = rsm_in%num_pulses
    END IF

  END SUBROUTINE merge_metadata_o

  !==============================================================================

#endif    /* NETCDF            */


END MODULE radar_obs_meta_read
