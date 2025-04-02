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


MODULE radar_output3d_ascii

!------------------------------------------------------------------------------
!
! Description: Top-level procedure of the radar forward operator EMVORADO for processing
!              of the volume data output in different formats (ascii, netcdf, grib2).
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
  
  USE radar_interface, ONLY : &
       abort_run,          &
       get_model_time_sec, &
       get_model_time_ddhhmmss, &
       get_obs_time_tolerance, &
       get_datetime_ini, &
       it_is_time_for_radar

  USE radar_utilities, ONLY : new_datetime,               &
                              diff_seconds,               &
                              get_free_funit

  USE radar_data, ONLY :         &
       cmaxlen, noutstreams_max, &
       radar_meta_type,          &
       t_dbzcalc_params,         &
       ydate_ini_mod,            &
       idom,                     &
       my_radar_id

  USE radar_data_io, ONLY : MAX_LEN_DESC, t_grib, t_grib2_modelspec
  
  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, &
       loutdbz, loutpolstd, loutpolall, lextdbz, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, &
       lvoldata_output, &
       ydirradarout, &
       t_voldata_ostream, voldata_ostream, &
       lmds_z, lmds_vr

  USE radar_output_utils,  ONLY : get_fileprefix_ascii_output, get_next_key_from_pos, &
                                  replace_substr_with_value, eps_ascii

  USE radar_output_volscan, ONLY : prepare_write_volscan, write_cdfin_from_volscan, &
                                   write_grib_from_volscan
  
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

  PUBLIC ::  output3d_ascii_radar
  
  !==============================================================================
  ! Module variables

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !==============================================================================
  !-----------------------------------------------------------------------------------------------
  ! Output of a 3D polar radar data set as a 3D-ASCII-file or 3D Fortran binary file on the NEC SX:
  !
  ! ATTENTION: The index order for writing is now j,i,k! Formerly it was i,j,k!

  !-----------------------------------------------------------------------------------------------
  !  Wrapper-routine to enable a loop over different output streams
  !

  SUBROUTINE output3d_ascii_radar(itime, ni, nj, nk, zfield, &
       fieldname, fieldcomment, fieldunit, ylevtyp, rmeta, dbzmeta, is_aliased)

    !.. Input/Output parameters:

    INTEGER, INTENT(in) :: itime, ni, nj, nk
    REAL(KIND=dp), INTENT(in) :: zfield(ni,nj,nk)
    CHARACTER(len=*), INTENT(in) :: fieldname, fieldcomment, fieldunit, ylevtyp
    TYPE(radar_meta_type),  INTENT(in), OPTIONAL :: rmeta  ! is needed for usage within EMVORAOD, but might be left out if called from elsewhere
    TYPE(t_dbzcalc_params), INTENT(in), OPTIONAL :: dbzmeta ! is needed for usage within EMVORAOD, but might be left out if called from elsewhere
    LOGICAL, OPTIONAL, INTENT(in) :: is_aliased           ! is needed for usage within EMVORAOD, but might be left out if called from elsewhere

    !.. Local variables:
    integer :: i

    DO i=1, noutstreams_max
      ! An output stream is "active" if its output_list contains any non-empty entries (output variables):
      IF ( ANY( LEN_TRIM(voldata_ostream(i)%output_list(:)) > 0) ) THEN
        CALL output3d_ascii_radar_loop(itime, ni, nj, nk, zfield, &
             fieldname, fieldcomment, fieldunit, ylevtyp, voldata_ostream(i), rmeta, dbzmeta, is_aliased)
      END IF
    END DO

    
  END SUBROUTINE output3d_ascii_radar

  !-----------------------------------------------------------------------------------------------
  !  Routine for output of one output stream
  !
  
  SUBROUTINE output3d_ascii_radar_loop(itime, ni, nj, nk, zfield, &
       fieldname, fieldcomment, fieldunit, ylevtyp, volstream, rmeta, dbzmeta, is_aliased)

#ifdef WITH_ZLIB
    USE, INTRINSIC :: ISO_C_BINDING
#endif

    IMPLICIT NONE

    !.. Input/Output parameters:

    INTEGER, INTENT(in) :: itime, ni, nj, nk
    REAL(KIND=dp), INTENT(in) :: zfield(ni,nj,nk)
    TYPE(t_voldata_ostream) , INTENT(in) :: volstream
    CHARACTER(len=*), INTENT(in) :: fieldname, fieldcomment, fieldunit, ylevtyp
    TYPE(radar_meta_type),  INTENT(in), OPTIONAL :: rmeta  ! is needed for usage within EMVORAOD, but might be left out if called from elsewhere
    TYPE(t_dbzcalc_params), INTENT(in), OPTIONAL :: dbzmeta ! is needed for usage within EMVORAOD, but might be left out if called from elsewhere
    LOGICAL, OPTIONAL, INTENT(in) :: is_aliased           ! is needed for usage within EMVORAOD, but might be left out if called from elsewhere

    !.. Local variables:

    CHARACTER(len=*), PARAMETER :: yzroutine ='output3d_ascii_radar'

    REAL(KIND=dp)  :: time_mod  ! time in seconds since model start
    CHARACTER(len=300) :: ausdateiname
    INTEGER            :: i, j, k, kk, ios, funit, nk_out
    INTEGER            :: ind_ele_outlist(1:nk)

    CHARACTER (LEN= 17)          :: yforrange      ! Forecast range as character
    CHARACTER (len=MAX_LEN_DESC) :: header1, header2, aliasstr, configstr

    CHARACTER (len=3)  :: cnk
    CHARACTER (len=29) :: zacttime
    CHARACTER (len=80) :: cfieldname, tmpfield
    CHARACTER (len=cmaxlen) :: errmsg
    INTEGER            :: ierr

    INTEGER, PARAMETER :: SP = selected_real_kind(6) ! 4byte kind parameter
    REAL(SP), ALLOCATABLE :: zfield_glob_32bit(:,:)
    REAL(SP) :: d32bit = 0.0
    CHARACTER(len=1) :: dimen
    CHARACTER(len=5) :: cnrbuf
    LOGICAL     :: name_is_in_list, is_supobed, write_error, force_overwrite

#ifdef WITH_ZLIB
    INTEGER(kind=c_int)     :: gzw, gzc, gzierror
    type (c_ptr)            :: funit_gz = c_null_ptr
    CHARACTER(kind=C_CHAR, len=:), ALLOCATABLE :: lbuffer  ! ALLOCATABLE chars are an F2003 feature!
    TYPE(c_ptr)             :: gzErrorMessage_Cptr
    INTEGER                 :: len_lbuf, pos_lbuf
    CHARACTER(len=5)        :: mode_gz = 'w'
    CHARACTER(len=40)       :: fstr
#endif

    ! To pass informations on model run and date/time to volscan output (cdfin and gribs):
    TYPE(t_grib)                :: grib_input_for_datetime
    TYPE(t_grib2_modelspec)     :: model_locinfo

#ifdef WITH_ZLIB
    INTERFACE
      FUNCTION gzOpen(path,mode) BIND(C, NAME='gzopen')
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(kind=c_char) :: path(*),mode(*)
        TYPE(c_ptr) :: gzOpen
      END FUNCTION gzOpen
    END INTERFACE

    INTERFACE
      FUNCTION gzGets(file,buf,len) BIND(C, NAME='gzgets')
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(kind=c_char) :: buf(*)
        TYPE(c_ptr) :: gzGets
        INTEGER(kind=c_int),VALUE :: file,len
      END FUNCTION gzGets
    END INTERFACE

    INTERFACE
      FUNCTION gzPuts(file,buf) BIND(C, NAME='gzputs')
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(kind=c_char) :: buf(*)
        INTEGER(kind=c_int) :: gzPuts
        TYPE(c_ptr),VALUE :: file
      END FUNCTION gzPuts
    END INTERFACE

    INTERFACE
      FUNCTION gzWrite(file,buf,len) BIND(C, NAME='gzwrite')
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(kind=c_char) :: buf(*)
        INTEGER(kind=c_int) :: gzWrite
        INTEGER(kind=c_int),VALUE :: len
        TYPE(c_ptr),VALUE :: file
      END FUNCTION gzWrite
    END INTERFACE

    INTERFACE
      FUNCTION gzClose(file) BIND(C, NAME='gzclose')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(kind=c_int) :: gzClose
        TYPE(c_ptr),VALUE :: file
      END FUNCTION gzClose
    END INTERFACE

    INTERFACE
      FUNCTION gzError(file,ierror) BIND(C, NAME='gzerror')
        USE, INTRINSIC :: ISO_C_BINDING
!        CHARACTER(kind=c_char) :: gzError
        TYPE(c_ptr) :: gzError
        INTEGER(kind=c_int) :: ierror
        TYPE(c_ptr),VALUE :: file
      END FUNCTION gzError
    END INTERFACE
#endif

    !--------------------------------------------------------------------------------------------
    ! .. switch to deactive ASCII-output:

    IF (.NOT.lvoldata_output) RETURN

    !--------------------------------------------------------------------------------------------
    ! 1) Check if output is requested of the current field by comparing with volstream%output_list:

    !   a) Extract the "true" field name from input character fieldname,
    !      which is the string before the "_id-...", if the latter occurs:
    cfieldname(:) = ' '
    DO i=1, LEN_TRIM(fieldname)
      IF (i <= LEN_TRIM(fieldname)-3) THEN
        IF (fieldname(i:i+3) /= '_id-') THEN
          cfieldname(i:i) = fieldname(i:i)
        ELSE
          EXIT
        END IF
      ELSE
        cfieldname(i:i) = fieldname(i:i)
      END IF
    END DO

#ifndef AUXOUT_OFFLINE
    !   b) Compare with list elements:
    name_is_in_list = .FALSE.
    DO i=1, SIZE(volstream%output_list)
      IF ( TRIM(volstream%output_list(i)) == TRIM(cfieldname) ) THEN
        name_is_in_list = .TRUE.
      END IF
    END DO
    IF ( .NOT. name_is_in_list ) THEN
      IF ( ldebug_radsim ) THEN
        WRITE (*,'(a)') 'INFO radar_src: '//TRIM(cfieldname)// &
             ' is not in namelist-parameter ''volstream%output_list'' and output is therefore skipped.'
      END IF
      RETURN
    END IF
#endif

    !--------------------------------------------------------------------------------------------
    ! 2) Check if this is a requested output time:

    ! Check if the substring "sup" is part of cfieldname. If yes, it is a superobed field and
    !  the time checking should be done with the rmeta%obs_times_fdbk rather than
    !  rmeta%obs_times_voldata:

    is_supobed = .FALSE.
    DO i=1, LEN_TRIM(cfieldname)-2
      IF (cfieldname(i:i+2) == 'sup') is_supobed = .TRUE.
    END DO

    IF (PRESENT(rmeta)) THEN
      IF (is_supobed) THEN
        IF ( .NOT.it_is_time_for_radar (rmeta%obs_times_fdbk) ) THEN
          RETURN
        END IF
      ELSE
        IF ( .NOT.it_is_time_for_radar (rmeta%obs_times_voldata) ) THEN
          RETURN
        END IF
      END IF
    END IF


    !--------------------------------------------------------------------------------------------
    ! 3) Filename for the output-file:

    ausdateiname(:) = ' '

    time_mod     = get_model_time_sec()

    SELECT CASE (TRIM(volstream%FORMAT))

    CASE ('cdfin','cdfin-mulmom')

      IF (.NOT.PRESENT(rmeta)) THEN
        CALL abort_run (my_radar_id, 10071, &
             'ERROR: problem in call to '//TRIM(yzroutine)//': '''// &
             TRIM(volstream%FORMAT)//''' output for '//TRIM(cfieldname)// &
             ' not possible due to missing OPTIONAL argument ''rmeta''! Fix the source code!', &
             TRIM(yzroutine)//', writing output file')
      END IF

!!$   Default naming convention:  "cdfin_z_id-010908_201307282200_201307282255_volscan"     or
!!$                               "cdfin_z_id-010908_201307282200_201307282255_precipscan"

      CALL ostream_create_filename (volstream, rmeta, time_mod, cfieldname, ausdateiname, zacttime, yforrange)
      IF (volstream%output_subdir(1:1) == '/') THEN
        ! If volstream%output_subdir starts with a /, it is interpreted as an absolute path name:
        ausdateiname = TRIM(volstream%output_subdir)//TRIM(ausdateiname)
      ELSE
        ! Otherwise, it is a subdir to ydirradarout, which has already a '/' at the end:
        IF (LEN_TRIM(volstream%output_subdir) == 0) THEN
          ausdateiname = TRIM(ydirradarout)//TRIM(ausdateiname)
        ELSE
          ausdateiname = TRIM(ydirradarout)//TRIM(volstream%output_subdir)//'/'//TRIM(ausdateiname)
        END IF
      END IF
      
    CASE ('grib2','grib2-mulmom')

      IF (.NOT.PRESENT(rmeta)) THEN
        CALL abort_run (my_radar_id, 10071, &
             'ERROR: problem in call to '//TRIM(yzroutine)//': '''// &
             TRIM(volstream%FORMAT)//''' output for '//TRIM(cfieldname)// &
             ' not possible due to missing OPTIONAL argument ''rmeta''! Fix the source code!', &
             TRIM(yzroutine)//', writing output file')
      END IF

!!$   Default naming convention:  "grib2_id-010908_201307282200_201307282255_volscan"     or
!!$                               "grib2_id-010908_201307282200_201307282255_precipscan"

      CALL ostream_create_filename (volstream, rmeta, time_mod, cfieldname, ausdateiname, zacttime, yforrange)
      IF (volstream%output_subdir(1:1) == '/') THEN
        ! If volstream%output_subdir starts with a /, it is interpreted as an absolute path name:
        ausdateiname = TRIM(volstream%output_subdir)//TRIM(ausdateiname)
      ELSE
        ! Otherwise, it is a subdir to ydirradarout, which has already a '/' at the end::
        IF (LEN_TRIM(volstream%output_subdir) == 0) THEN
          ausdateiname = TRIM(ydirradarout)//TRIM(ausdateiname)
        ELSE
          ausdateiname = TRIM(ydirradarout)//TRIM(volstream%output_subdir)//'/'//TRIM(ausdateiname)
        END IF
      END IF
      
    CASE default

      ! Time strings for the global config string below:
      yforrange(:)   = ' '
      yforrange      = get_model_time_ddhhmmss()
      
      zacttime(:)    = ' '
      zacttime(1:14) = new_datetime( ydate_ini_mod, rmeta%obs_times(itime) )

      ! Only basename, extention will be added below for each format:

      IF (volstream%output_subdir(1:1) == '/') THEN
        ! If volstream%output_subdir starts with a /, it is interpreted as an absolute path name:
        ausdateiname = TRIM(volstream%output_subdir)// &
             TRIM(ADJUSTL(fieldname))//'_'//TRIM(ydate_ini_mod)//'_'//TRIM(yforrange)//'_'//ylevtyp
      ELSE
        ! Otherwise, it is a subdir to ydirradarout, which has already a '/' at the end::
        IF (LEN_TRIM(volstream%output_subdir) == 0) THEN
          ausdateiname = TRIM(ydirradarout)// &
               TRIM(ADJUSTL(fieldname))//'_'//TRIM(ydate_ini_mod)//'_'//TRIM(yforrange)//'_'//ylevtyp
        ELSE
          ausdateiname = TRIM(ydirradarout)//TRIM(volstream%output_subdir)//'/'// &
               TRIM(ADJUSTL(fieldname))//'_'//TRIM(ydate_ini_mod)//'_'//TRIM(yforrange)//'_'//ylevtyp
        END IF
      END IF

    END SELECT

    !--------------------------------------------------------------------------------------------
    ! 4) Prepare elevation meta data for output:

    IF (PRESENT(rmeta)) THEN
      IF (is_supobed) THEN
        nk_out = rmeta%nel_fdbk
        ! Return, if no elevations are desired for output:
        IF (nk_out <= 0) RETURN
        ind_ele_outlist(1:nk_out) = rmeta%ind_ele_fdbk(1:nk_out)
      ELSE
        nk_out = rmeta%nel_voldata
        ! Return, if no elevations are desired for output:
        IF (nk_out <= 0) RETURN
        ind_ele_outlist(1:nk_out) = rmeta%ind_ele_voldata(1:nk_out)
      END IF
    ELSE
      nk_out = nk
      ind_ele_outlist(1:nk) = (/ (k,k=1,nk) /)
    END IF

    !--------------------------------------------------------------------------------------------
    ! 5) Write metadata and data to file:

    CALL get_free_funit (funit)
    IF (funit == -1) THEN
      CALL abort_run (my_radar_id, 10071, &
           'ERROR: problem in '//TRIM(yzroutine)//': no free file unit available! Abort!', &
           TRIM(yzroutine)//', opening output file')
    END IF

    ! Prepare the header information depending on the presence of radar meta structure
    !   rmeta:
    header1(:) = ' '
    header2(:) = ' '
    configstr(:) = ' '
    IF (PRESENT(rmeta) .AND. PRESENT(dbzmeta)) THEN
      ! if metadata structure rmeta is present, write some extra
      !   radar metadata into the file header:      
      WRITE (configstr, &
           '(5(a,"=",a,1x),1(a,"=",i6.6,1x),'//&
           '1(a,"=",L1,"(",i2.2,"x",i2.2,")"),10(1x,a,"=",L1),'//&
           '2(1x,a,"=",f0.8),2(1x,a,"=",f0.1),3(1x,a,"=",f0.3),'//&
           '1(1x,a,"=",i2.2))') &
           'parameter', TRIM(ADJUSTL(cfieldname)), &
           'time', TRIM(ADJUSTL(zacttime)), &
           'inidate_model', TRIM(ADJUSTL(ydate_ini_mod)), &
           'forecasttime_model', TRIM(ADJUSTL(yforrange)), &
           'station_name', TRIM(ADJUSTL(rmeta%station_name)), &
           'station_id', rmeta%station_id, &
           'lsmooth', lsmooth, rmeta%ngpsm_h, rmeta%ngpsm_v, &
           'lonline', lonline, &
           'lsode', lsode, &
           'lextdbz', lextdbz, &
           'loutpolstd', loutpolstd, &
           'loutpolall', loutpolall, &
           'llookup_mie', dbzmeta%llookup_mie, &
           'lweightdbz', lweightdbz, &
           'lmds_z', lmds_z, &
           'lmds_vr', lmds_vr, &
           'lfall', lfall, &
           'radar_lon', rmeta%lon, &
           'radar_lat', rmeta%lat, &
           'radar_alt_msl_mod', rmeta%alt_msl, &
           'radar_alt_msl_true', rmeta%alt_msl_true, &
           'ra_inc', rmeta%ra_inc, &
           'az_start', rmeta%az_start, &
           'az_inc', rmeta%az_inc, &
           'itype_refl', dbzmeta%itype_refl

      IF (PRESENT(is_aliased)) THEN
        WRITE (configstr, '(a,1x,a,"=",f0.3,1x,a,"=",L1)') &
             TRIM(configstr), 'vnyq(ele1,time1)', rmeta%ext_nyq(1,1), 'is_aliased', is_aliased
      END IF
      WRITE (header1, '(a,1x,a)') &
           '# ASCII '//TRIM(ADJUSTL(fieldcomment))//' ['//TRIM(ADJUSTL(fieldunit))//']', &
           TRIM(configstr)
      WRITE (cnk(1:3), '(i3)') nk_out
      WRITE (header2, '("|",'//TRIM(cnk)//'(1x,f0.2))') &
           (rmeta%el_arr(ind_ele_outlist(k)), k=1,nk_out)
    ELSE
      ! no extra information in header:
      WRITE (header1, '(a)') &
           '# ASCII '//TRIM(ADJUSTL(fieldcomment))//' ['//TRIM(ADJUSTL(fieldunit))//']'
    END IF


    SELECT CASE (TRIM(volstream%format))

    CASE ('f90-binary')

!!$ could be changed to form='unformatted', access='stream', in which case no record separators are
!!$  written inbetween the single output records. This would make it totally compatible with C-style binary files.

      IF (ldebug_radsim) THEN
        WRITE (*,'(A)') TRIM(ausdateiname)//'.bin'
      END IF

      OPEN(unit=funit, file=TRIM(ausdateiname)//'.bin', status='replace', &
           form='unformatted', iostat=ios)
      IF (ios /= 0) THEN
        CALL abort_run (my_radar_id, 10072, &
             'ERROR: problem in '//TRIM(yzroutine)//': error opening '//TRIM(ausdateiname)//'.bin', &
             TRIM(yzroutine)//', opening output file')
      ENDIF

      cnrbuf(:) = ' '
      WRITE (cnrbuf,'(i5.5)') LEN_TRIM(header1)
      WRITE (funit) cnrbuf
      WRITE (funit) TRIM(header1)

      IF (LEN_TRIM(header2) > 0) THEN
        dimen = ACHAR(3+100)
        WRITE (funit) dimen
        WRITE (cnrbuf,'(i5.5)') LEN_TRIM(header2)
        WRITE (funit) cnrbuf
        WRITE (funit) nj, ni, nk_out
        WRITE (funit) TRIM(header2)
      ELSE
        dimen = ACHAR(3)
        WRITE (funit) dimen
        WRITE (funit) nj, ni, nk_out  ! nj=range, ni=azimut, nk_out=elevation
      END IF

      ALLOCATE(zfield_glob_32bit(nj,ni))
      zfield_glob_32bit = 0.0
      DO k=1,nk_out
        kk = ind_ele_outlist(k)
        zfield_glob_32bit = TRANSPOSE(REAL(zfield(:,:,kk), kind=KIND(d32bit)))
        WRITE (funit) zfield_glob_32bit
      END DO
      DEALLOCATE(zfield_glob_32bit)

      CLOSE (funit)

    CASE ('ascii')

      WRITE (*,'(A)') TRIM(ausdateiname)//'.dat'
      OPEN(unit=funit, file=TRIM(ausdateiname)//'.dat', status='replace', iostat=ios)
      IF (ios /= 0) THEN
        CALL abort_run (my_radar_id, 10072, &
             'ERROR: problem in '//TRIM(yzroutine)//': error opening '//TRIM(ausdateiname)//'.dat', &
             TRIM(yzroutine)//', opening output file')
      ENDIF

      WRITE (funit, '(a)') TRIM(header1)
      IF (LEN_TRIM(header2) > 0) THEN
        WRITE (funit, '(i4,1x,i4,1x,i4,1x,a)') nj, ni, nk_out, TRIM(header2)
      ELSE
        WRITE (funit, '(i4,1x,i4,1x,i4)') nj, ni, nk_out
      END IF

      DO k=1,nk_out  ! elevations
        kk = ind_ele_outlist(k)
        DO i=1,ni    ! azimut
          DO j=1,nj  ! range
            IF (ABS(zfield(i,j,kk)) >= eps_ascii) THEN
              WRITE (funit, '(es12.5)') zfield(i,j,kk)
            ELSE
              WRITE (funit, '(i1)') 0
            END IF
          END DO
        END DO
      END DO

      CLOSE (funit)

#ifdef WITH_ZLIB

! compile with -DWITH_ZLIB -lz
    CASE ('ascii-gzip')

      !--------------------------------------------------------------------------------------------
      ! write in one long allocatable string using implicit loop, then call gzputs

      IF (ldebug_radsim) THEN
        WRITE (*,'(A)') TRIM(ausdateiname)//'.dat.gz'
      END IF

      funit_gz = gzopen(TRIM(ausdateiname)//'.dat.gz'//c_null_char, TRIM(mode_gz)//c_null_char)
      IF (.NOT.C_ASSOCIATED(funit_gz)) THEN
        CALL abort_run (my_radar_id, 10072, &
             'ERROR: problem in '//TRIM(yzroutine)//': error opening '//TRIM(ausdateiname)//'.dat.gz', &
             TRIM(yzroutine)//', opening output file with gzopen()')
      end if

      len_lbuf = LEN_TRIM(header1) + LEN_TRIM(header2) + 3*5 + nj*ni*nk_out*(12+1) + 1 + 100
      ALLOCATE(CHARACTER(len=len_lbuf) :: lbuffer)

      !--------------------------------------------------------------------------------------------
      ! .. Put header information to string lbuffer:

      lbuffer(1:len_lbuf) = ' '  ! need explicit address range, otherwise re-allocation with lenght 1!
      lbuffer(1:LEN_TRIM(header1)+1) = TRIM(header1)//c_new_line  ! need explicit address range, otherwise re-allocation
      pos_lbuf = LEN_TRIM(lbuffer) + 1
      IF (LEN_TRIM(header2) > 0) THEN
        WRITE (lbuffer(pos_lbuf:len_lbuf), '(i4,1x,i4,1x,i4,1x,a)') nj, ni, nk_out, TRIM(header2)
      ELSE
        WRITE (lbuffer(pos_lbuf:len_lbuf), '(i4,1x,i4,1x,i4)') nj, ni, nk_out
      END IF
      pos_lbuf = LEN_TRIM(lbuffer) + 1
      lbuffer(pos_lbuf:pos_lbuf) = c_new_line
      pos_lbuf = LEN_TRIM(lbuffer) + 1

      !--------------------------------------------------------------------------------------------
      ! .. Put data to lbuffer:
      fstr(:) = ' '
      WRITE (fstr, '("(",i0,"(es12.5,A))")') ni*nj*nk_out
      WRITE( lbuffer(pos_lbuf:len_lbuf), TRIM(fstr) )  &
           (((checkzero(zfield(i,j,ind_ele_outlist(k))), c_new_line, j=1,nj), i=1,ni), k=1,nk_out)

      !--------------------------------------------------------------------------------------------
      ! .. Put lbuffer to file and clean up:
      errmsg(:) = ' '
      gzw = gzputs(funit_gz, TRIM(lbuffer)//c_null_char)
      gzErrorMessage_Cptr = gzError(funit_gz,gzierror)
      IF (gzierror /= 0) THEN
        errmsg = C_to_F_string(gzErrorMessage_Cptr, LEN(errmsg))
        CALL abort_run (my_radar_id, 10072, &
             'ERROR: problem in '//TRIM(yzroutine)//': error writing '//TRIM(ausdateiname)//'.dat.gz', &
             TRIM(yzroutine)//', call to gzPuts(): '//TRIM(errmsg))
      END IF

      gzc = gzclose(funit_gz)

      DEALLOCATE (lbuffer)


      ! A version with explicitly setting small values to 0, not 0.0000E+00 with an outer loop over writing each
      ! number in to lbuffer(pos_lbuf:pos_lbuf+12) has been tested, but is much slower and did not
      ! lead to significantly better compression and data reduction.

#endif

   CASE ('cdfin','cdfin-mulmom') ! Output of DWD's cdfin-Format (NetCDF, simular to results of burfx2netcdf from BUFR-data of CIRRUS database)

#ifdef NETCDF
     IF (PRESENT(rmeta)) THEN
       IF (ldebug_radsim) THEN
         WRITE (*,'(A)') TRIM(ausdateiname)
       END IF
       CALL prepare_write_volscan (rmeta, itime, TRIM(cfieldname), grib_input_for_datetime, model_locinfo, ierr, errmsg)
       
       IF (ABS(volstream%content_dt) < 1e-4_dp .AND. TRIM(volstream%FORMAT) /= 'cdfin-mulmom') THEN
         ! If output file changes every timestep and every moment is in a separate file (single-moment-single-volume-file),
         !  we can safely overwrite any pre-existing files with same name, because they must come from a previous
         !  model run and the user did forget to delete them. In other situations this is not so easy to determine.
         force_overwrite = .TRUE.
       ELSE
         force_overwrite = .FALSE.
       END IF

       CALL write_cdfin_from_volscan(rmeta, itime, grib_input_for_datetime, model_locinfo, &
            &                        TRIM(ausdateiname), force_overwrite, TRIM(cfieldname), &
            &                        TRIM(configstr), TRIM(fieldunit), TRIM(fieldcomment), zfield, &
            &                        ind_ele_outlist, nk_out)
     END IF
#else
      CALL abort_run (my_radar_id, 10083, &
           'ERROR in '//TRIM(yzroutine)//': output format '//TRIM(volstream%format)// &
           ' needs compilation with -DNETCDF and netcdf4 library!', &
           TRIM(yzroutine)//', ''cdfin''-format needs ''-DNETCDF''')
#endif

   CASE ('grib2','grib2-mulmom') ! Output of DWD's grib2-Format

#ifdef GRIBAPI
     IF (PRESENT(rmeta)) THEN
       IF (ldebug_radsim) THEN
         WRITE (*,'(A)') TRIM(ausdateiname)
       END IF
       CALL prepare_write_volscan (rmeta, itime, TRIM(cfieldname), grib_input_for_datetime, model_locinfo, ierr, errmsg)

       IF (ABS(volstream%content_dt) < 1e-4_dp .AND. TRIM(volstream%format) /= 'grib2-mulmom') THEN
         ! If output file changes every timestep and every moment is in a separate file (single-moment-single-volume-file),
         !  we can safely overwrite any pre-existing files with same name, because they must come from a previous
         !  model run and the user did forget to delete them. In other situations this is not so easy to determine.
         force_overwrite = .TRUE.
       ELSE
         force_overwrite = .FALSE.
       END IF
       
       SELECT CASE (TRIM(cfieldname))
       CASE('zrsim','zrobs')
         CALL write_grib_from_volscan(rmeta, itime, grib_input_for_datetime, model_locinfo, &
                                      TRIM(ausdateiname), force_overwrite, &
                                      TRIM(cfieldname), TRIM(volstream%grib2_packingtype), &
                                      TRIM(header1), zfield, &
                                      ind_ele_outlist, nk_out, ierr, errmsg)
       CASE default
         ! other variables not yet available in grib2 definitions
         ierr = 0
         errmsg(:) = ' '
       END SELECT
       IF (ierr /= 0) THEN
         CALL abort_run (my_radar_id, 10483, TRIM(errmsg), TRIM(yzroutine))
       END IF
     END IF
#else
      CALL abort_run (my_radar_id, 10094, &
           'ERROR in '//TRIM(yzroutine)//': output format '//TRIM(volstream%format)// &
           ' needs compilation with -DGRIBAPI and eccodes library!', &
           TRIM(yzroutine)//', ''cdfin''-format needs ''-DGRIBAPI''')
#endif

   CASE default

     CALL abort_run (my_radar_id, 10073, &
          'ERROR: problem in '//TRIM(yzroutine)//': output format '//TRIM(volstream%format)// &
          ' not implemented', &
          TRIM(yzroutine)//', namelist parameter ''volstream%format''')

   END SELECT

  CONTAINS

#ifdef WITH_ZLIB
    PURE REAL(KIND=dp) FUNCTION checkzero(x)
      REAL(KIND=dp), INTENT(in) :: x
      IF (ABS(x) >= eps_ascii) THEN
        checkzero = x
      ELSE
        checkzero = 0.0_dp
      END IF
      RETURN
    END FUNCTION checkzero

    FUNCTION C_to_F_string(c_string_pointer, maxlen) RESULT(f_string)

!!$      USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr,c_f_pointer,c_char,c_null_char

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

#endif

  END SUBROUTINE output3d_ascii_radar_loop

  !=========================================================================
  !
  ! Function for checking if a parameter name represents a geometry-related
  ! parameter.
  !
  !=========================================================================

  FUNCTION is_geometry_parameter(ncdf_shortname) RESULT (flag)
    CHARACTER(len=*), INTENT(in) :: ncdf_shortname
    LOGICAL                      :: flag
    INTEGER                      :: i
    ! List of all implemented geometry related vars (need equal string lengths!):
    CHARACTER(len=*), PARAMETER  :: testlist(7) = ['losim   ', 'losupsim', 'lasim   ', &
         'lasupsim', 'hrsim   ', 'ersim   ', 'adsim   ']

    flag = .FALSE.
    DO i=1, SIZE(testlist)
      flag = flag .OR. (TRIM(ncdf_shortname) == TRIM(testlist(i)))
    END DO
  END FUNCTION is_geometry_parameter
  
  !=========================================================================
  !
  ! Procedure for creating output file names depending on the file format
  ! and patterns given in the t_voldata_ostream type
  !
  !=========================================================================

  SUBROUTINE ostream_create_filename (volstream, rmeta, time_mod, cfieldname, ofilename, opt_timerange, opt_vvzrange)

    TYPE(t_voldata_ostream) , INTENT(in)    :: volstream
    TYPE(radar_meta_type), INTENT(in)       :: rmeta
    REAL(KIND=dp), INTENT(in)               :: time_mod   ! time in seconds since model start
    CHARACTER(len=*), INTENT(in)            :: cfieldname ! parameter name
    CHARACTER(len=*), INTENT(out)           :: ofilename  ! output file name
    CHARACTER(len=*), INTENT(out), OPTIONAL :: opt_timerange  ! absolute range-of-time string for usage
                                                              !  in the calling PROGRAM ("YYYYMMDDhhmmss-YYYYMMDDhhmmss")
    CHARACTER(len=*), INTENT(out), OPTIONAL :: opt_vvzrange   ! relative range-of-time string for usage
                                                              !  in the calling PROGRAM ("ddhhmmss-ddhhmmss")

    CHARACTER(len=*), PARAMETER :: yzroutine = 'emvorado:ostream_create_filename'
    CHARACTER(len=14) :: content_starttime, content_endtime, model_starttime
    CHARACTER(len=8)  :: ddhhmmss_starttime, ddhhmmss_endtime
    CHARACTER(len=30) :: formatprefix, scantypesuffix
    CHARACTER(len=LEN(volstream%file_pattern)) :: fpat
    CHARACTER(len=20) :: key, cstation_id
    REAL(KIND=dp)     :: content_tstart, content_tend, tol_obs_times
    INTEGER           :: i, j, ku, ko, pos, isecdiff, ierr_diff(2,6)
    LOGICAL           :: is_singlemom, is_singletime
    LOGICAL           :: have_time, have_id, have_scantype, have_varname


    ofilename(:) = ' '

    !--------------------------------------------------------------------------
    ! .. Create filename in a way that times are organized in batches of length
    !     "content_dt" (sedonds) relative to a reference time "content_tref"
    !     (seconds since model start):
    
    model_starttime = get_datetime_ini()

    IF (volstream%content_dt >= 1e-4_dp) THEN
      content_tstart = FLOOR( (time_mod-volstream%content_tref+1e-6_dp) / volstream%content_dt ) * &
           volstream%content_dt + volstream%content_tref
      is_singletime = .FALSE.
    ELSE
      content_tstart = time_mod
      is_singletime = .TRUE.
    END IF
    content_tend = content_tstart + volstream%content_dt

    tol_obs_times = get_obs_time_tolerance( idom )
    content_starttime = 'YYYYMMDDHHmmss'
    content_endtime   = 'YYYYMMDDHHmmss'
    j = 1
    DO i=1, rmeta%nobs_times
      IF (rmeta%obs_times(i) >= content_tstart - tol_obs_times ) THEN
        content_starttime = rmeta%obs_cdate(i)
        content_endtime   = content_starttime  ! to make sure there is the correct endtime in case of content_dt = 0.0
        j = i
        EXIT
      END IF
    END DO

    DO i=j, rmeta%nobs_times
      IF (rmeta%obs_times(i) < content_tend - tol_obs_times) THEN
        content_endtime = rmeta%obs_cdate(i)
      END IF
    END DO

    ! Forecast timestamps from absolute timestamps:
    CALL diff_seconds ( model_starttime, content_starttime, isecdiff, ierr_diff )
    ddhhmmss_starttime = get_model_time_ddhhmmss( REAL(isecdiff,dp) )
    CALL diff_seconds ( model_starttime, content_endtime, isecdiff, ierr_diff )
    ddhhmmss_endtime   = get_model_time_ddhhmmss( REAL(isecdiff,dp) )

    IF (PRESENT(opt_timerange)) THEN
      opt_timerange(:) = ' '
      IF (is_singletime) THEN
        opt_timerange = content_starttime
      ELSE
        opt_timerange = TRIM(content_starttime)//'-'//TRIM(content_endtime)
      END IF
    END IF
    IF (PRESENT(opt_vvzrange)) THEN
      opt_vvzrange(:) = ' '
      IF (is_singletime) THEN
        opt_vvzrange = ddhhmmss_starttime
      ELSE
        opt_vvzrange = TRIM(ddhhmmss_starttime)//'-'//TRIM(ddhhmmss_endtime)
      END IF
    END IF
    
    !--------------------------------------------------------------------------
    ! .. Some format-depended settings:
    
    formatprefix(:) = ' '
    scantypesuffix(:) = ' '
    SELECT CASE (TRIM(volstream%FORMAT))
    CASE ('cdfin','cdfin-mulmom')
      formatprefix = 'cdfin'
      IF (TRIM(volstream%FORMAT) == 'cdfin') THEN
        is_singlemom = .TRUE.
      ELSE
        is_singlemom = .FALSE.
      END IF
    CASE ('grib2','grib2-mulmom')
      formatprefix = 'grib2'
      IF (TRIM(volstream%FORMAT) == 'grib2') THEN
        is_singlemom = .TRUE.
      ELSE
        is_singlemom = .FALSE.
      END IF
    CASE default
      CALL abort_run (my_radar_id, 19073, &
           'ERROR '//TRIM(yzroutine)//': output format '//TRIM(volstream%FORMAT)// &
           ' not implemented', &
           TRIM(yzroutine)//', input parameter ''volstream%format''')
    END SELECT

    IF (LEN_TRIM(volstream%file_pattern) <= 0) THEN

      !--------------------------------------------------------------------------
      ! Set the default file name pattern depending on the file format

      ! Only basename, scantype suffix will be added below:
      IF (is_singlemom) THEN
        WRITE (ofilename, '(a,"_",a,"_id-",i6.6,2("_",a))') TRIM(formatprefix), TRIM(cfieldname), rmeta%station_id,  &
             content_starttime(1:12), content_endtime(1:12)
      ELSE
        ! multi-moment files, but different files for geometry (lon, lat, height, local elevation), obs and sim:
        IF (is_geometry_parameter(TRIM(cfieldname))) THEN
          WRITE (ofilename, '(a,"_",a,"_id-",i6.6,2("_",a))') TRIM(formatprefix), 'allgeom', rmeta%station_id,  &
               content_starttime(1:12), content_endtime(1:12)            
        ELSE IF (INDEX(TRIM(cfieldname),'obs') > 0) THEN
          WRITE (ofilename, '(a,"_",a,"_id-",i6.6,2("_",a))') TRIM(formatprefix), 'allobs', rmeta%station_id,  &
               content_starttime(1:12), content_endtime(1:12)
        ELSE
          WRITE (ofilename, '(a,"_",a,"_id-",i6.6,2("_",a))') TRIM(formatprefix), 'allsim', rmeta%station_id,  &
               content_starttime(1:12), content_endtime(1:12)
        END IF
      END IF

      SELECT CASE (TRIM(rmeta%scanname))
      CASE ('PRECIP')
        scantypesuffix = 'precipscan'
      CASE default
        scantypesuffix = 'volscan'
      END SELECT

      ofilename = TRIM(ofilename)//'_'//TRIM(scantypesuffix)

    ELSE

      !--------------------------------------------------------------------------
      ! Set the user defined file name pattern

      ! .. List of valid <keys> for filename patterns:
      !
      !       - <stationid>    : will result in e.g. 'id-123456'
      !       - <varname>      : will result in e.g. 'zrsim' or 'all', depending on format
      !       - <tmodelini>    : will result in YYYYMMDDhhmmss  (absolute datetime of model start)
      !       - <tstart>       : will result in YYYYMMDDhhmmss  (absolute datetime of start of time window of file content)
      !       - <tend>         : will result in YYYYMMDDhhmmss  (absolute datetime of end of time window of file content)
      !       - <tact>         : will result in YYYYMMDDhhmmss  (absolute datetime of actual model time)
      !       - <tvvzstart>    : will result in DDhhmmss  (forecast lead time of start of time window of file content)
      !       - <tvvzend>      : will result in DDhhmmss  (forecast lead time of end of time window of file content)
      !       - <tvvzact>      : will result in DDhhmmss  (forecast lead time of actual model time)
      !       - <scantype>     : will result in 'volscan' or 'precipscan' by default, but can be defined by pat_scantype_XXX below
      !
      !    Depending on the format, some of them are optional, but most are mandatory.
      
      fpat = volstream%file_pattern ! copy original pattern

      have_time     = .FALSE.
      have_id       = .FALSE.
      have_scantype = .FALSE.
      have_varname  = .FALSE.
      
      IF (is_singletime) THEN
        pos = 1
        DO
          ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
          CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
          IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
          ! .. and, if it is a time key, replace by its actual value:
          SELECT CASE (TRIM(key))
          CASE ('<tmodelini>')
            ! optional:
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(model_starttime))
          CASE ('<tact>','<tstart>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(content_starttime))
          CASE ('<tend>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(content_endtime))
          CASE ('<tvvzact>','<tvvzstart>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(ddhhmmss_starttime))
          CASE ('<tvvzend>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(ddhhmmss_endtime))
          END SELECT
        END DO
      ELSE
        pos = 1
        DO
          ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
          CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
          IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
          ! .. and, if it is a time key, replace by its actual value:
          SELECT CASE (TRIM(key))
          CASE ('<tmodelini>')
            ! optional:
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(model_starttime))
          CASE ('<tstart>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(content_starttime))
          CASE ('<tend>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(content_endtime))
          CASE ('<tvvzstart>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(ddhhmmss_starttime))
          CASE ('<tvvzend>')
            have_time = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(ddhhmmss_endtime))
          END SELECT
        END DO
      END IF

      IF (is_singlemom) THEN
        pos = 1
        DO
          ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
          CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
          IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
          ! .. and, if it is a varname key, replace by its actual value:
          SELECT CASE (TRIM(key))
          CASE ('<varname>')
            have_varname = .TRUE.
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(cfieldname))
          END SELECT
        END DO
        
      ELSE
        IF (TRIM(formatprefix) == 'grib2') THEN
          have_varname = .TRUE. ! for grib2-mulmom, varname is not mandatory in filenames
        END IF
        pos = 1
        DO
          ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
          CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
          IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
          ! .. and, if it is a varname key, replace by its actual value:
          SELECT CASE (TRIM(key))
          CASE ('<varname>')
            have_varname = .TRUE.
            IF (is_geometry_parameter(TRIM(cfieldname))) THEN       
              CALL replace_substr_with_value (fpat, ku, ko, 'allgeom')
            ELSE IF (INDEX(TRIM(cfieldname),'obs') > 0) THEN
              CALL replace_substr_with_value (fpat, ku, ko, 'allobs')
            ELSE
              CALL replace_substr_with_value (fpat, ku, ko, 'allsim')
            END IF
          END SELECT
        END DO
      END IF

      ! .. Station-ID and scantype are always mandatory:
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. and, if it is a stationid or scantype key, replace by its actual value:
        SELECT CASE (TRIM(key))
        CASE ('<stationid>')
          have_id = .TRUE.
          cstation_id(:) = ' '
          WRITE (cstation_id, '("id-",i6.6)') rmeta%station_id
          CALL replace_substr_with_value (fpat, ku, ko, TRIM(cstation_id))
        CASE ('<scantype>')
          have_scantype = .TRUE.
          SELECT CASE (TRIM(rmeta%scanname))
          CASE ('PRECIP')
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(volstream%pat_scantype_precipscan))
          CASE default
            CALL replace_substr_with_value (fpat, ku, ko, TRIM(volstream%pat_scantype_volscan))
          END SELECT
        END SELECT
      END DO

      ! At last, check if any other keys or "<" or ">" are still remaining in the fpat. If yes, it is a wrong key:
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. if there is still a key left in fpat, it is a wrong key:
        CALL abort_run (my_radar_id, 19074, &
             'ERROR '//TRIM(yzroutine)//': volstream%file_pattern "'//TRIM(volstream%file_pattern)// &
             '" contains unknown key '//TRIM(key), &
             TRIM(yzroutine)//', input parameter ''volstream%file_pattern''')
      END DO

      ! .. if still some remaining single '<' or '>' characters are in fpat, there is an error:
      IF (INDEX(fpat,'<') > 0 .OR. INDEX(fpat,'>') > 0) THEN
        CALL abort_run (my_radar_id, 19075, &
             'ERROR '//TRIM(yzroutine)//': volstream%file_pattern "'//TRIM(volstream%file_pattern)// &
             '" contains garbled "<" or ">" characters', &
             TRIM(yzroutine)//', input parameter ''volstream%file_pattern''')
      END IF

      ! .. check if all mandatory keys have been specified in volstream%file_pattern:
      IF (.NOT.have_time) THEN
        CALL abort_run (my_radar_id, 19076, &
             'ERROR '//TRIM(yzroutine)//': volstream%file_pattern "'//TRIM(volstream%file_pattern)// &
             '" does not contain a <key> for time (<tact>, <tvvzact>, <tstart>, <tend>, <tvvzstart>, or <tvvzend>)', &
             TRIM(yzroutine)//', input parameter ''volstream%file_pattern''')
      END IF

      IF (.NOT.have_id) THEN
        CALL abort_run (my_radar_id, 19077, &
             'ERROR '//TRIM(yzroutine)//': volstream%file_pattern "'//TRIM(volstream%file_pattern)// &
             '" does not contain the mandatory key <stationid>', &
             TRIM(yzroutine)//', input parameter ''volstream%file_pattern''')
      END IF

      IF (.NOT.have_scantype) THEN
        CALL abort_run (my_radar_id, 19078, &
             'ERROR '//TRIM(yzroutine)//': volstream%file_pattern "'//TRIM(volstream%file_pattern)// &
             '" does not contain the mandatory key <scantype>', &
             TRIM(yzroutine)//', input parameter ''volstream%file_pattern''')
      END IF

      IF (.NOT.have_varname) THEN
        CALL abort_run (my_radar_id, 19079, &
             'ERROR '//TRIM(yzroutine)//': volstream%file_pattern "'//TRIM(volstream%file_pattern)// &
             '" does not contain the mandatory key <varname>', &
             TRIM(yzroutine)//', input parameter ''volstream%file_pattern''')
      END IF

      ! .. at this stage, fpat contains a valid filename:
      ofilename = TRIM(fpat)
      
    END IF

  END SUBROUTINE ostream_create_filename

END MODULE radar_output3d_ascii
