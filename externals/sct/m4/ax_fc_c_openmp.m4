# AX_FC_C_OPENMP
# ---------------------------------------------------------------------
# Runs AC_OPENMP for C and Fortran and checks for a linking flag for
# Fortran compiler that is required to link C code that uses OpenMP.
# The flag is set to the OPENMP_FC_C_FLAGS and cached as ax_cv_prog_fc_c_openmp.
# The result is either "none needed", "unsupported" (and OPENMP_FC_C_FLAGS is empty)
# or the actual flags (OPEN_FC_C_FLAGS is equal to the flags).
m4_define([_AX_FC_C_OPENMP],
[
      program main
      implicit none
      interface
        subroutine foo()  BIND(C, NAME='foo')
        end subroutine foo
      end interface
      call foo
      end
])

AC_DEFUN([AX_FC_C_OPENMP], [
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_FC])dnl
AC_LANG_PUSH([C])
AC_OPENMP
AC_LANG_POP([C])
AC_LANG_PUSH([Fortran])
AC_OPENMP
AC_LANG_POP([Fortran])
AC_CACHE_CHECK([for $FC option to link code that requires OpenMP and is compiled with $CC], [ax_cv_prog_fc_c_openmp],
  [OPENMP_FC_C_FLAGS=
  ax_cv_prog_fc_c_openmp=unsupported
  AS_IF([test x"$ac_cv_prog_c_openmp" != xunsupported],
    [AC_LANG_PUSH([C])
    ax_fc_c_openmp_CFLAGS_save=$CFLAGS
    CFLAGS="$OPENMP_CFLAGS $ax_fc_c_openmp_CFLAGS_save"
    ac_compile="$ac_compile && cp conftest.$ac_objext conftest_foo.$ac_objext"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#ifndef _OPENMP
 choke me
#endif
#include <omp.h>
int foo (void) { return omp_get_num_threads (); }]])],
      [AC_LANG_PUSH([Fortran])
      ax_fc_c_openmp_LIBS_save=$LIBS
      LIBS="conftest_foo.$ac_objext $ax_fc_c_openmp_LIBS_save"
      AC_LINK_IFELSE([_AX_FC_C_OPENMP],
        [ax_cv_prog_fc_c_openmp='none needed'],
        [ax_fc_c_openmp_FCFLAGS_save=$FCFLAGS
        AS_IF([test x"$ac_cv_prog_fc_openmp" != xunsupported],
          [FCFLAGS="$OPENMP_FCFLAGS $ax_fc_c_openmp_FCFLAGS_save"
          AC_LINK_IFELSE([_AX_FC_C_OPENMP],
            [ax_cv_prog_fc_c_openmp=$ac_cv_prog_fc_openmp],
            [AS_IF([test x"$ac_cv_prog_c_openmp" != "xnone needed"],
              [ax_cv_prog_fc_c_openmp=`AS_ECHO(["$OPENMP_CFLAGS"]) | sed 's%[^ ]* *%-Wl,&%g'`
              FCFLAGS="$ax_cv_prog_fc_c_openmp $ax_fc_c_openmp_FCFLAGS_save"
              AC_LINK_IFELSE([_AX_FC_C_OPENMP], [],
              [ax_cv_prog_fc_c_openmp=unsupported])])])])
        FCFLAGS=$ax_fc_c_openmp_FCFLAGS_save])
      LIBS=$ax_fc_c_openmp_LIBS_save
      AC_LANG_POP([Fortran])])
    CFLAGS=$ax_fc_c_openmp_CFLAGS_save
    AC_LANG_POP([C])])])
AS_CASE([$ax_cv_prog_fc_c_openmp],
  ["none needed" | unsupported], [],
  [OPENMP_FC_C_FLAGS=$ax_cv_prog_fc_c_openmp])
AC_SUBST([OPENMP_FC_C_FLAGS])
])
