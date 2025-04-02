PROGRAM test_sct
#ifdef HAVE_MPI
  USE mpi
#endif
  USE sct, ONLY: sct_init, sct_new_timer, sct_del_timer, sct_start, sct_stop, &
       &         sct_resolution, sct_val, sct_report, sct_new_context,        &
       &         sct_context_start, sct_context_stop, sct_active, SCT_GETENV, &
       &         sct_set_callstats

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  INTEGER :: timer1, timer2

  CALL init0

  timer1 = sct_new_timer('timer 1')
  timer2 = sct_new_timer('timer 2')

  call sct_start(timer2)
  CALL test1
  call sct_stop(timer2)

  !CALL sct_report(proc_choice=SCT_GETENV, thread_choice=SCT_GETENV, sp_merging=SCT_GETENV)
  CALL sct_report()
  !CALL sct_report(timer1)

  CALL exit0

CONTAINS

  SUBROUTINE test1
    REAL(dp) :: s, x, dt
    INTEGER :: i,j

    DO j = 1, 8
      s = 0.0_dp
      IF (sct_active(timer1)) CALL fail("test1:(a): (timer_active(timer1) /= 0)")
      CALL sct_start(timer1)
      IF (.NOT. sct_active(timer1)) CALL fail("test1:(b): (timer_active(timer1) /= 1)")
      DO i = 0, j*100000
        x = REAL(i,dp)
        s = s + EXP(-SQRT(x))
      ENDDO
      CALL sct_stop(timer1)
      dt = sct_val(timer1)
      IF (sct_active(timer1)) CALL fail("test1:(c): (timer_active(timer1) /= 0)")
      IF (s<2.0_dp) PRINT*,'j, s, dt=',j,s,dt ! fake interest in the result
    ENDDO

  END SUBROUTINE test1

  SUBROUTINE exit0
    INTEGER :: ierr
#ifdef HAVE_MPI
    CALL mpi_finalize(ierr)
#endif
  END SUBROUTINE exit0

  SUBROUTINE init0
    INTEGER :: ierr
#ifdef HAVE_MPI
    CALL MPI_INIT(ierr)
    CALL check_mpi(ierr)
#endif

    CALL sct_init(timer_max=3)
    CALL sct_set_callstats(1)

  END SUBROUTINE init0

  SUBROUTINE check_mpi(rc)
    INTEGER, INTENT(in) :: rc

#ifdef HAVE_MPI
    IF (rc /= MPI_SUCCESS) THEN
      WRITE(0,*) 'MPI return code  /= MPI_SUCCESS'
      STOP
    END IF
#endif

  END SUBROUTINE check_mpi

  SUBROUTINE fail(msg)
    CHARACTER(len=*), INTENT(in):: msg
    INTEGER :: ierr

    WRITE(0,*) 'fail: ',msg
#ifdef HAVE_MPI
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
    STOP

  END SUBROUTINE fail

END PROGRAM test_sct
