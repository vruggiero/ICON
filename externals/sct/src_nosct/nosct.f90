! Maintainer: Joerg Behrens <behrens@dkrz.de>
!             Hendryk Bockelmann <bockelmann@dkrz.de>
! Copyright: Copyright (C) Deutsches Klimarechenzentrum GmbH (DKRZ)
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


!> the "fake" interface/module of sct for all FORTRAN users
!> does not provide any timers but enables the user to use the same interfaces
MODULE sct
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE, C_NULL_CHAR

  IMPLICIT NONE

  !
  ! fortran dummy module for disabled sct
  !

  PRIVATE

  PUBLIC :: sct_init, sct_new_timer, sct_del_timer, sct_start, sct_stop,           &
       &    sct_reset_timer, sct_resolution, sct_val, sct_report, sct_new_context, &
       &    sct_context_start, sct_context_stop, sct_active, sct_last_dt,          &
       &    sct_stop_all, sct_set_callstats, sct_set_eventcounters,                &
       &    sct_set_nestedtimers

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

  INTERFACE sct_new_context
    MODULE PROCEDURE new_global_context
    MODULE PROCEDURE new_general_context
  END INTERFACE sct_new_context

CONTAINS

  SUBROUTINE sct_set_callstats(val)
    IMPLICIT NONE
    INTEGER, INTENT(in):: val
  END SUBROUTINE sct_set_callstats

  SUBROUTINE sct_set_eventcounters(val)
    IMPLICIT NONE
    INTEGER, INTENT(in):: val
  END SUBROUTINE sct_set_eventcounters

  SUBROUTINE sct_set_nestedtimers(val)
    IMPLICIT NONE
    INTEGER, INTENT(in):: val
  END SUBROUTINE sct_set_nestedtimers

  SUBROUTINE sct_init(timer_max, default_context_name, default_comm)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: timer_max
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: default_context_name
    INTEGER, OPTIONAL, INTENT(in) :: default_comm
  END SUBROUTINE sct_init

  SUBROUTINE sct_del_timer(itimer)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: itimer
  END SUBROUTINE sct_del_timer

  LOGICAL FUNCTION sct_active(itimer)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: itimer
    sct_active = .FALSE.
  END FUNCTION sct_active

  SUBROUTINE sct_start(itimer)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: itimer
  END SUBROUTINE sct_start

  SUBROUTINE sct_stop(itimer)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: itimer
  END SUBROUTINE sct_stop

  SUBROUTINE sct_stop_all()
    IMPLICIT NONE
  END SUBROUTINE sct_stop_all

  SUBROUTINE sct_context_start(icon)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: icon
  END SUBROUTINE sct_context_start

  SUBROUTINE sct_context_stop(icon)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: icon
  END SUBROUTINE sct_context_stop

  SUBROUTINE sct_reset_timer(itimer)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: itimer
  END SUBROUTINE sct_reset_timer

  REAL(C_DOUBLE) FUNCTION sct_val(itimer)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    sct_val = 0.0_C_DOUBLE
  END FUNCTION sct_val

  REAL(C_DOUBLE) FUNCTION sct_resolution()
    IMPLICIT NONE
    sct_resolution =  0.0_C_DOUBLE
  END FUNCTION sct_resolution

  INTEGER FUNCTION new_global_context(label)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: label
    new_global_context = 0
  END FUNCTION new_global_context

  INTEGER FUNCTION new_general_context(label, comm)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: label
    INTEGER, INTENT(in) :: comm
    new_general_context = 0
  END FUNCTION new_general_context

  INTEGER FUNCTION sct_new_timer(label)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: label
    sct_new_timer = 0
  END FUNCTION sct_new_timer

  SUBROUTINE sct_report(timer_choice, proc_choice, thread_choice, sp_merging)
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(in) :: timer_choice, proc_choice, thread_choice, sp_merging
  END SUBROUTINE sct_report

  REAL(C_DOUBLE) FUNCTION sct_last_dt(itimer)
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(in) :: itimer
    sct_last_dt = 0.0_C_DOUBLE
  END FUNCTION sct_last_dt

END MODULE sct
