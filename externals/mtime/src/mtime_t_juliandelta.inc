!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
  ! ================================================================================
  ! juliandelta section:
  ! ================================================================================

  ! generic assignment for constructors
  !
  FUNCTION t_juliandelta_assign_raw(sign, day, ms)
    TYPE(t_juliandelta) :: t_juliandelta_assign_raw
    CHARACTER(c_char), INTENT(in) :: sign
    INTEGER(c_int64_t), INTENT(in) :: day
    INTEGER(c_int64_t), INTENT(in) :: ms
    TYPE(c_ptr) :: c_pointer
    TYPE(juliandelta), POINTER :: jd_tmp

    c_pointer = my_newjuliandelta(sign, day, ms)
    !print *,sign, c_pointer
    IF (.NOT. C_ASSOCIATED(c_pointer)) THEN
      CALL handle_errno(1*100 + 1, __FILE__, __LINE__)
      t_juliandelta_assign_raw%jd%sign = 'L'
    ELSE
      CALL C_F_POINTER(c_pointer, jd_tmp)
      t_juliandelta_assign_raw%jd = jd_tmp
      t_juliandelta_assign_raw%jd%sign = jd_tmp%sign
      CALL my_deallocatejuliandelta(c_pointer)
    END IF
  END FUNCTION t_juliandelta_assign_raw
