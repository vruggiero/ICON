! @copyright Copyright (C) 2017 Deutsches Klimarechenzentrum GmbH (DKRZ)
!
! @author Jörg Behrens <behrens@dkrz.de>
!         Hendryk Bockelmann <bockelmann@dkrz.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#if HAVE_FCONFIG_H
#  include <fconfig.h>
#endif

!> the interface/module of sct for all FORTRAN users
MODULE sct
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_LONG, C_FLOAT, C_DOUBLE, C_NULL_CHAR
#ifdef HAVE_MPI
  USE MPI
#endif
  IMPLICIT NONE

  !
  ! fortran wrapper module for the c implementation of sct
  !

  PRIVATE

  PUBLIC :: sct_init, sct_new_timer, sct_del_timer, sct_start, sct_stop,           &
       &    sct_reset_timer, sct_reset_all, sct_resolution, sct_val, sct_report,   &
       &    sct_new_context, sct_context_start, sct_context_stop, sct_active,      &
       &    sct_last_dt, sct_stop_all, sct_set_callstats, sct_set_eventcounters,   &
       &    sct_set_nestedtimers, sct_add_report_attribute

#ifdef HAVE_MPI
  INTEGER, PARAMETER, PUBLIC :: SCT_COMM_WORLD = MPI_COMM_WORLD
#else
  INTEGER, PARAMETER, PUBLIC :: SCT_COMM_WORLD = 0
#endif

  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_WITHOUT_CALLSTATS = 0
  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_WITH_CALLSTATS = 1
  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_DEFAULT_CALLSTATS = 2

  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_SP_MERGE_SIMPLE = 1
  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_SP_SERIAL_ONLY = 2
  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_SP_PARALLEL_ONLY = 4
  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_SP_SELECT_ALL = 8

  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_SELECT_ALL = -1
  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_REDUCE_ALL = -2

  INTEGER(C_INT), PARAMETER, PUBLIC :: SCT_GETENV = -255

  !> generic interface for functions new_global_context and new_general_context
  INTERFACE sct_new_context
    MODULE PROCEDURE new_global_context
    MODULE PROCEDURE new_general_context
  END INTERFACE sct_new_context

  !> generic interface for sct_add_report_attribute_* family
  INTERFACE sct_add_report_attribute
    MODULE PROCEDURE sct_add_report_attribute_int_f
    MODULE PROCEDURE sct_add_report_attribute_long_f
    MODULE PROCEDURE sct_add_report_attribute_float_f
    MODULE PROCEDURE sct_add_report_attribute_double_f
    MODULE PROCEDURE sct_add_report_attribute_string_f
  END INTERFACE sct_add_report_attribute

  !> interface for C-function ::sct_set_callstats
  INTERFACE
    !> en-/dis-able callstats output
    !! \param[in] val 0 = disable output, 1 = enable output
    SUBROUTINE sct_set_callstats(val)  BIND(C, NAME='sct_set_callstats')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: val
    END SUBROUTINE sct_set_callstats
  END INTERFACE

  !> interface for C-function ::sct_set_eventcounters
  INTERFACE
    !> en-/dis-able eventcounters output
    !! \param[in] val 0 = disable output, 1 = enable output
    SUBROUTINE sct_set_eventcounters(val)  BIND(C, NAME='sct_set_eventcounters')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: val
    END SUBROUTINE sct_set_eventcounters
  END INTERFACE

  !> interface for C-function ::sct_set_nestedtimers
  INTERFACE
    !> en-/dis-able nestedtimers output
    !! \param[in] val 0 = disable output, 1 = enable output
    SUBROUTINE sct_set_nestedtimers(val)  BIND(C, NAME='sct_set_nestedtimers')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: val
    END SUBROUTINE sct_set_nestedtimers
  END INTERFACE

  !> interface for C-function ::sct_del_timer
  INTERFACE
    !> delete timer if not active
    !! \param[in] itimer handle of timer to be deleted
    SUBROUTINE sct_del_timer(itimer)  BIND(C, NAME='sct_del_timer')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    END SUBROUTINE sct_del_timer
  END INTERFACE

  !> interface for C-function ::sct_start
  INTERFACE
    !> start timer if not active
    !! \param[in] itimer handle of timer to be started
    SUBROUTINE sct_start(itimer)  BIND(C, NAME='sct_start')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    END SUBROUTINE sct_start
  END INTERFACE

  !> interface for C-function ::sct_stop
  INTERFACE
    !> stop timer if active
    !! \param[in] itimer handle of timer to be stopped
    SUBROUTINE sct_stop(itimer)  BIND(C, NAME='sct_stop')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    END SUBROUTINE sct_stop
  END INTERFACE

  !> interface for C-function ::sct_stop_all
  INTERFACE
    !> stop all active timer
    SUBROUTINE sct_stop_all()  BIND(C, NAME='sct_stop_all')
      IMPLICIT NONE
    END SUBROUTINE sct_stop_all
  END INTERFACE

  !> interface for C-function ::sct_context_start
  INTERFACE
    !> initialise a context switch
    !! \param[in] icon handle of context to be used
    !! \detail from default context to context icon
    SUBROUTINE sct_context_start(icon)  BIND(C, NAME='sct_context_start')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: icon
    END SUBROUTINE sct_context_start
  END INTERFACE

  !> interface for C-function ::sct_context_stop
  INTERFACE
    !> finalise a context switch
    !! \param[in] icon handle of context to stop
    !! \detail from context icon back to default context
    SUBROUTINE sct_context_stop(icon)  BIND(C, NAME='sct_context_stop')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: icon
    END SUBROUTINE sct_context_stop
  END INTERFACE

  !> interface for C-function ::sct_reset_timer
  INTERFACE
    !> reset all statistics of timer
    !! \param[in] itimer handle of timer to be reset
    SUBROUTINE sct_reset_timer(itimer)  BIND(C, NAME='sct_reset_timer')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    END SUBROUTINE sct_reset_timer
  END INTERFACE

  !> interface for C-function ::sct_reset_all
  INTERFACE
    !> reset all statistics of all timer
    SUBROUTINE sct_reset_all()  BIND(C, NAME='sct_reset_all')
      IMPLICIT NONE
    END SUBROUTINE sct_reset_all
  END INTERFACE

  !> interface for C-function ::sct_val
  INTERFACE
    !> get actual value of timer (accumulated time of all calls)
    !! \param[in] itimer handle of timer to be evaluated
    !! \return actual timer value
    REAL(C_DOUBLE) FUNCTION sct_val(itimer)  BIND(C, NAME='sct_val')
      IMPORT:: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    END FUNCTION sct_val
  END INTERFACE

  !> interface for C-function ::sct_last_dt
  INTERFACE
    !>  get value of last timed interval
    !! \param[in] itimer handle of timer to be evaluated
    !! \return time of last timed interval
    REAL(C_DOUBLE) FUNCTION sct_last_dt(itimer)  BIND(C, NAME='sct_last_dt')
      IMPORT:: C_INT, C_DOUBLE
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    END FUNCTION sct_last_dt
  END INTERFACE

  !> interface for C-function ::sct_resolution
  INTERFACE
    !> get timer resolution of underlying system
    !! \return timer resolution
    REAL(C_DOUBLE) FUNCTION sct_resolution()  BIND(C, NAME='sct_resolution')
      IMPORT:: C_INT, C_DOUBLE
      IMPLICIT NONE
    END FUNCTION sct_resolution
  END INTERFACE

CONTAINS

  !> Fortran wrapper for sct_add_report_attribute_int
  SUBROUTINE sct_add_report_attribute_int_f(key, val)
    CHARACTER(len=*), INTENT(in) :: key
    INTEGER(C_INT), INTENT(in) :: val

    INTERFACE
      SUBROUTINE sct_add_report_attribute_int(key, val)  BIND(C, NAME='sct_add_report_attribute_int')
        IMPORT:: C_CHAR, C_INT
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: key
        INTEGER(C_INT), VALUE, INTENT(in) :: val
      END SUBROUTINE sct_add_report_attribute_int
    END INTERFACE

    CALL sct_add_report_attribute_int(TRIM(ADJUSTL(key)) // c_null_char, val)

  END SUBROUTINE sct_add_report_attribute_int_f

  !> Fortran wrapper for sct_add_report_attribute_long
  SUBROUTINE sct_add_report_attribute_long_f(key, val)
    CHARACTER(len=*), INTENT(in) :: key
    INTEGER(C_LONG), INTENT(in) :: val

    INTERFACE
      SUBROUTINE sct_add_report_attribute_long(key, val)  BIND(C, NAME='sct_add_report_attribute_long')
        IMPORT:: C_CHAR, C_LONG
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: key
        INTEGER(C_LONG), VALUE, INTENT(in) :: val
      END SUBROUTINE sct_add_report_attribute_long
    END INTERFACE

    CALL sct_add_report_attribute_long(TRIM(ADJUSTL(key)) // c_null_char, val)

  END SUBROUTINE sct_add_report_attribute_long_f

  !> Fortran wrapper for sct_add_report_attribute_float
  SUBROUTINE sct_add_report_attribute_float_f(key, val)
    CHARACTER(len=*), INTENT(in) :: key
    REAL(C_FLOAT), INTENT(in) :: val

    INTERFACE
      SUBROUTINE sct_add_report_attribute_float(key, val)  BIND(C, NAME='sct_add_report_attribute_float')
        IMPORT:: C_CHAR, C_FLOAT
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: key
        REAL(C_FLOAT), VALUE, INTENT(in) :: val
      END SUBROUTINE sct_add_report_attribute_float
    END INTERFACE

    CALL sct_add_report_attribute_float(TRIM(ADJUSTL(key)) // c_null_char, val)

  END SUBROUTINE sct_add_report_attribute_float_f

  !> Fortran wrapper for sct_add_report_attribute_double
  SUBROUTINE sct_add_report_attribute_double_f(key, val)
    CHARACTER(len=*), INTENT(in) :: key
    REAL(C_DOUBLE), INTENT(in) :: val

    INTERFACE
      SUBROUTINE sct_add_report_attribute_double(key, val)  BIND(C, NAME='sct_add_report_attribute_double')
        IMPORT:: C_CHAR, C_DOUBLE
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: key
        REAL(C_DOUBLE), VALUE, INTENT(in) :: val
      END SUBROUTINE sct_add_report_attribute_double
    END INTERFACE

    CALL sct_add_report_attribute_double(TRIM(ADJUSTL(key)) // c_null_char, val)

  END SUBROUTINE sct_add_report_attribute_double_f

  !> Fortran wrapper for sct_add_report_attribute_string
  SUBROUTINE sct_add_report_attribute_string_f(key, val)
    CHARACTER(len=*), INTENT(in) :: key
    CHARACTER(len=*), INTENT(in) :: val

    INTERFACE
      SUBROUTINE sct_add_report_attribute_string(key, val)  BIND(C, NAME='sct_add_report_attribute_string')
        IMPORT:: C_CHAR, C_DOUBLE
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: key
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: val
      END SUBROUTINE sct_add_report_attribute_string
    END INTERFACE

    CALL sct_add_report_attribute_string(TRIM(ADJUSTL(key)) // c_null_char, TRIM(ADJUSTL(val)) // c_null_char)

  END SUBROUTINE sct_add_report_attribute_string_f

  !> \brief initialisation function of sct
  !! \param[in] timer_max maximum number of timers to be used
  !! \param[in] default_context_name name of the context to be used
  !! \param[in] default_comm communicator used for this context
  !! \details This function needs to be called in the very beginning before any timer can be used. It sets the whole needed environment.
  SUBROUTINE sct_init(timer_max, default_context_name, default_comm)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: timer_max
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: default_context_name
    INTEGER, OPTIONAL, INTENT(in) :: default_comm

    INTEGER :: dcomm
    CHARACTER(len=127) :: dcn

    INTERFACE
      SUBROUTINE my_init_timer_aux(timer_max, default_context_name, default_comm)  BIND(C, NAME='sct_init_f')
        IMPORT:: C_CHAR, C_INT
        IMPLICIT NONE
        INTEGER(C_INT), VALUE :: timer_max, default_comm
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: default_context_name
      END SUBROUTINE my_init_timer_aux
    END INTERFACE

    IF (PRESENT(default_context_name)) THEN
      dcn = default_context_name
    ELSE
      dcn = 'default context'
    ENDIF

    IF (PRESENT(default_comm)) THEN
      dcomm = default_comm
    ELSE
      dcomm = SCT_COMM_WORLD
    ENDIF

    CALL my_init_timer_aux(timer_max, TRIM(ADJUSTL(dcn)) // c_null_char, dcomm)

  END SUBROUTINE sct_init

  !> \brief check whether timer is in use
  !! \param[in] itimer handle of timer to be checked
  !! \return true if timer is already started, false is timer is stopped
  LOGICAL FUNCTION sct_active(itimer)
    INTEGER, INTENT(in) :: itimer

    INTERFACE
      INTEGER(C_INT) FUNCTION my_timer_active_c(itimer)  BIND(C, NAME='sct_active')
        IMPORT:: C_INT
        IMPLICIT NONE
        INTEGER(C_INT), VALUE, INTENT(in) :: itimer
      END FUNCTION my_timer_active_c
    END INTERFACE

    IF (my_timer_active_c(itimer) /= 0) THEN
      sct_active = .TRUE.
    ELSE
      sct_active = .FALSE.
    ENDIF
  END FUNCTION sct_active

  !> \brief defines a new global context
  !! \param[in] label name of new context
  !! \return handle for a new context; handle is always > 0 since 0 is used for the default_context
  INTEGER FUNCTION new_global_context(label)
    CHARACTER(len=*), INTENT(in) :: label

    INTERFACE
      INTEGER(C_INT) FUNCTION new_global_context_aux(label)  BIND(C, NAME='sct_new_global_context')
        IMPORT:: C_CHAR, C_INT
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: label
      END FUNCTION new_global_context_aux
    END INTERFACE

    new_global_context = new_global_context_aux(TRIM(ADJUSTL(label)) // c_null_char)

  END FUNCTION new_global_context

  !> \brief defines a new general context, not neccessarily global
  !! \param[in] label name of new context
  !! \param[in] comm communicator to be used for reduction within this context
  !! \return handle for a new context; handle is always > 0 since 0 is used for the default_context
  INTEGER FUNCTION new_general_context(label, comm)
    CHARACTER(len=*), INTENT(in) :: label
    INTEGER, INTENT(in) :: comm

    INTERFACE
      INTEGER(C_INT) FUNCTION new_general_context_aux(label, comm)  BIND(C, NAME='sct_new_context_f')
        IMPORT:: C_CHAR, C_INT
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: label
        INTEGER(C_INT), INTENT(in) :: comm
      END FUNCTION new_general_context_aux
    END INTERFACE

    new_general_context = new_general_context_aux(TRIM(ADJUSTL(label)) // c_null_char, comm )

  END FUNCTION new_general_context

  !> \brief initialise a new timer
  !! \param[in] label name of timer to be used for report
  !! \return handle of new timer
  INTEGER FUNCTION sct_new_timer(label)
    CHARACTER(len=*), INTENT(in) :: label

    INTERFACE
      INTEGER(C_INT) FUNCTION my_new_timer_aux(label)  BIND(C, NAME='sct_new_timer')
        IMPORT:: C_CHAR, C_INT
        IMPLICIT NONE
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: label
      END FUNCTION my_new_timer_aux
    END INTERFACE

    sct_new_timer = my_new_timer_aux(TRIM(ADJUSTL(label)) // c_null_char )

  END FUNCTION sct_new_timer

  !> \brief print timer results
  !! \param[in] timer_choice SCT_SELECT_ALL or timer_id for output
  !! \param[in] proc_choice SCT_SELECT_ALL or SCT_REDUCE_ALL MPI-tasks for output
  !! \param[in] thread_choice SCT_SELECT_ALL or SCT_REDUCE_ALL OpenMP-threads for output
  !! \param[in] sp_merging merging operation for serial part of OpenMP timer
  SUBROUTINE sct_report(timer_choice, proc_choice, thread_choice, sp_merging)
    INTEGER, OPTIONAL, INTENT(in) :: timer_choice, proc_choice, thread_choice, sp_merging
    INTEGER :: u, ios
    LOGICAL :: opened

    INTEGER :: timer_choice2, proc_choice2, thread_choice2, sp_merging2

    INTERFACE
      SUBROUTINE my_report_aux(timer_choice, proc_choice, thread_choice, sp_merging)  BIND(C, NAME='sct_single_report')
        IMPORT:: C_INT
        IMPLICIT NONE
        INTEGER(C_INT), VALUE, INTENT(in) :: timer_choice, proc_choice, thread_choice, sp_merging
      END SUBROUTINE my_report_aux
    END INTERFACE

    ! code to check which units we might want to flush

    DO u=0, 102
      INQUIRE(unit=u,opened=opened)
      IF (opened) THEN
        FLUSH(unit=u, iostat=ios)
      ENDIF
    ENDDO

    ! arguments

    IF (PRESENT(timer_choice)) THEN
      timer_choice2 = timer_choice
    ELSE
      timer_choice2 = SCT_SELECT_ALL
    ENDIF

    IF (PRESENT(proc_choice)) THEN
      proc_choice2 = proc_choice
    ELSE
      proc_choice2 = SCT_GETENV
    ENDIF

    IF (PRESENT(thread_choice)) THEN
      thread_choice2 = thread_choice
    ELSE
      thread_choice2 = SCT_GETENV
    ENDIF

    IF (PRESENT(sp_merging)) THEN
      sp_merging2 = sp_merging
    ELSE
      sp_merging2 = SCT_GETENV
    ENDIF

    CALL my_report_aux(timer_choice2, proc_choice2, thread_choice2, sp_merging2)

  END SUBROUTINE sct_report

END MODULE sct

