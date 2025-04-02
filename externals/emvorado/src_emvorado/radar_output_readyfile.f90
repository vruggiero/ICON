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

MODULE radar_output_readyfile

  !------------------------------------------------------------------------------
  !
  ! Description: Methods of the radar forward operator EMVORADO for writing
  !              "ready" files.
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
  
  USE radar_data, ONLY : cmaxlen
  
  USE radar_interface, ONLY : &
       abort_run,          &
       get_model_time_ddhhmmss, &
       get_datetime_act, &
       get_datetime_ini

  USE radar_utilities,    ONLY : get_free_funit
  USE radar_output_utils, ONLY : get_next_key_from_pos, replace_substr_with_value
  
!================================================================================
!================================================================================

  IMPLICIT NONE

!================================================================================
!================================================================================

  !==============================================================================
  ! Public and Private Subroutines

  PRIVATE

  PUBLIC ::  write_ready_radar

  !==============================================================================
  ! Module variables

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !============================================================================

  SUBROUTINE write_ready_radar ( path, file_pattern, content, ierror )

    CHARACTER(len=*), INTENT(in)  :: path, file_pattern, content
    INTEGER        ,  INTENT(out) :: ierror

    INTEGER :: funit, iret
    CHARACTER(len=cmaxlen) :: ofilename, rdy_filename, tmp_filename
    CHARACTER(len=*), PARAMETER :: yzroutine = 'emvorado::write_ready_radar'
    CHARACTER(len=*), PARAMETER :: tmp_prefix = ".."
    CHARACTER(len=14) :: model_starttime, model_validtime
    CHARACTER(len=8)  :: ddhhmmss_validtime
    CHARACTER(len=LEN(file_pattern)) :: fpat
    CHARACTER(len=20) :: key
    INTEGER           :: ku, ko, pos
    LOGICAL           :: have_time

    INTERFACE
      FUNCTION private_rename(old_filename, new_filename) RESULT(iret) BIND(C, NAME='rename')
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char
        INTEGER(c_int) :: iret
        CHARACTER(c_char), DIMENSION(*), INTENT(in) :: old_filename
        CHARACTER(c_char), DIMENSION(*), INTENT(in) :: new_filename
      END FUNCTION private_rename
    END INTERFACE

    ofilename(:) = ' '
    ierror = 0
    
    !--------------------------------------------------------------------------
    ! .. Create filename in a way that times are organized in batches of length
    !     "content_dt" (sedonds) relative to a reference time "content_tref"
    !     (seconds since model start):

    ddhhmmss_validtime = get_model_time_ddhhmmss()
    model_starttime    = get_datetime_ini()
    model_validtime    = get_datetime_act()
    
    ofilename(:) = ' '

    IF (LEN_TRIM(file_pattern) <= 0) THEN

      ! Default name:
      !ofilename = TRIM(path)//'READY_EMVORADO_'//TRIM(model_validtime)
      ofilename = 'READY_EMVORADO_'//TRIM(model_validtime)

    ELSE
      
      !--------------------------------------------------------------------------
      ! Set the user defined file name pattern

      ! .. List of valid <keys> for filename patterns:
      !
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
        END SELECT
      END DO

      ! At last, check if any other keys or "<" or ">" are still remaining in the fpat. If yes, it is a wrong key:
      pos = 1
      DO
        ! .. get next key '<...>' from fpat starting at position pos and increment pos to prepare for the next search:
        CALL get_next_key_from_pos (fpat, pos, key, ku, ko)
        IF (ku < 0 .OR. ko < 0 .OR. pos >= LEN_TRIM(fpat)) EXIT
        ! .. if there is still a key left in fpat, it is a wrong key:
        WRITE (*, '(a)') &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" contains unknown key '//TRIM(key)
        ierror = 1
        RETURN
      END DO

      ! .. if still some remaining single '<' or '>' characters are in fpat, there is an error:
      IF (INDEX(fpat,'<') > 0 .OR. INDEX(fpat,'>') > 0) THEN
        WRITE (*, '(a)') &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" contains garbled "<" or ">" characters'
        ierror = 2
        RETURN
      END IF

      ! .. check if all mandatory keys have been specified in file_pattern:
      IF (.NOT.have_time) THEN
        WRITE (*, '(a)') &
             'ERROR '//TRIM(yzroutine)//': file_pattern "'//TRIM(file_pattern)// &
             '" does not contain a <key> for actual time (<tact>, <tvvzact>)'
        ierror = 3
        RETURN
      END IF

      ! .. at this stage, fpat contains a valid filename:
      !ofilename = TRIM(path)//TRIM(fpat)
      ofilename = TRIM(fpat)

    END IF
    
    ! Actually create ready file.
    !
    ! This procedure is carried out in two steps: First, a file with
    ! the name "tmp_prefix+rdy_filename" is created. After the file
    ! has been closed, it is then renamed into "rdy_filename" in a
    ! second step.
    ! This detour is necessary when another process polls the output
    ! directory and relies on a "complete" ready file.
    tmp_filename = TRIM(path)//tmp_prefix//TRIM(ofilename)
    rdy_filename = TRIM(path)//TRIM(ofilename)
    CALL get_free_funit( funit )

    OPEN(unit=funit, file=TRIM(tmp_filename), status='replace', form='formatted', iostat=iret)
    IF (iret == 0) THEN
      WRITE (funit,*) content
      CLOSE (funit)
    ELSE
      WRITE (*, '(a)') 'ERROR '//TRIM(yzroutine)//': could not open '//TRIM(tmp_filename)
      ierror = iret
      RETURN
    END IF

    iret = util_rename(TRIM(tmp_filename), TRIM(rdy_filename))
    IF (iret /= 0) THEN
      WRITE (*, '(a)') 'ERROR '//TRIM(yzroutine)//': could not rename '//TRIM(tmp_filename)// &
           ' to '//TRIM(rdy_filename)
      ierror = iret
      RETURN
    END IF

  CONTAINS

    FUNCTION util_rename(old_filename, new_filename) RESULT(iret)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_null_char
      INTEGER :: iret
      CHARACTER(len=*), INTENT(in) :: old_filename
      CHARACTER(len=*), INTENT(in) :: new_filename
      iret = private_rename(TRIM(old_filename)//c_null_char, TRIM(new_filename)//c_null_char)
    END FUNCTION util_rename

  END SUBROUTINE write_ready_radar

END MODULE radar_output_readyfile
