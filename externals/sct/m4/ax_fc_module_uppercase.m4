# AX_FC_MODULE_UPPERCASE()
# -----------------------------------------------------------------------------
# Checks whether the Fortran module files have uppercase names. The result is
# "yes", "no" or "unknown".
#
# The result is stored in the FC_MODUPPERCASE output variable and is cached in
# the ax_cv_fc_module_uppercase variable.
#
AC_DEFUN([AX_FC_MODULE_UPPERCASE],
  [AC_CACHE_CHECK([whether Fortran 90 modules files get uppercase names],
     [ax_cv_fc_module_uppercase],
     [ax_cv_fc_module_uppercase=unknown
      AS_MKDIR_P([conftest.dir])
      cd conftest.dir
      AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[[      module conftest_module
      contains
      subroutine conftest_routine
      end subroutine
      end module]])],
        [AS_IF([test -f CONFTEST_MODULE.*], [ax_cv_fc_module_uppercase=yes],
           [test -f conftest_module.*], [ax_cv_fc_module_uppercase=no])])
      AC_LANG_POP([Fortran])
      cd ..
      rm -rf conftest.dir])
   AC_SUBST([FC_MODUPPERCASE], [$ax_cv_fc_module_uppercase])])
