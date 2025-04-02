PROGRAM test_hdf5_fortran
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_long, c_float, c_double
#ifdef HAVE_MPI
  USE mpi
#endif
  USE sct, ONLY: sct_init, sct_new_timer, sct_start, sct_stop, sct_report, &
       &         sct_set_callstats, sct_add_report_attribute

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  INTEGER :: itimer, ierr
  INTEGER(c_int)  :: val_int
  INTEGER(c_long) :: val_long
  REAL(c_float)   :: val_float
  REAL(c_double)  :: val_double
  CHARACTER(len=20) :: val_string

#ifdef HAVE_MPI
  CALL MPI_INIT(ierr)
#endif

  CALL sct_init(timer_max=3)

  itimer = sct_new_timer('test timer')

  CALL sct_start(itimer)
  CALL test1
  CALL sct_stop(itimer)

  val_int  = 123
  val_long = 456
  val_float = 0.5
  val_double = 0.25
  val_string = 'abc'

  CALL sct_add_report_attribute('test_att_int',val_int)
  CALL sct_add_report_attribute('test_att_long',val_long)
  CALL sct_add_report_attribute('test_att_float',val_float)
  CALL sct_add_report_attribute('test_att_double',val_double)
  CALL sct_add_report_attribute('test_att_string',val_string)

  CALL sct_report()

#ifdef HAVE_MPI
  CALL mpi_finalize(ierr)
#endif

CONTAINS

  SUBROUTINE test1
    REAL(dp) :: s, x, dt
    INTEGER :: i,j

    DO j = 1, 8
      s = 0.0_dp
      DO i = 0, j*100000
        x = REAL(i,dp)
        s = s + EXP(-SQRT(x))
      ENDDO
      IF (s<2.0_dp) PRINT*,'j, s =',j,s ! fake interest in the result
    ENDDO

  END SUBROUTINE test1


END PROGRAM test_hdf5_fortran
