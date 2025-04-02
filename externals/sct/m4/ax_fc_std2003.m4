# SYNOPSIS
#
#   AX_FC_STD2003([ACTION-IF-SUCCEED], [ACTION-IF-FAIL])
#
# DESCRIPTION
#
#   Try figuring out the fortran compiler flags to support 2003 standard ISO_C_BINDING
#   and set FCFLAGS if found.
#
# LICENSE
#
#   Copyright (c) 2013 Hendryk Bockelmann <bockelmann@dkrz.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_FC_STD2003],
[AC_CACHE_CHECK([for Fortran flag needed to compile std2003 ISO_C_BINDING], [ax_cv_fc_std2003],
  [ax_cv_fc_std2003="unknown"
   AC_LANG_PUSH([Fortran])
   ax_save_FCFLAGS=$FCFLAGS
   for ax_flag in none -std=f2003 -qlanglvl=2003std -std03; do
     case $ax_flag in
       none) FCFLAGS=$ax_save_FCFLAGS ;;
       *)    FCFLAGS="$ax_save_FCFLAGS $ax_flag" ;;
     esac
     AC_LINK_IFELSE([AC_LANG_PROGRAM([],[
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_FLOAT
      TYPE, BIND(C) :: MYFTYPE
          INTEGER(C_INT) :: I
          REAL(C_FLOAT) :: F
      END TYPE MYFTYPE
])],
                    [ax_cv_fc_std2003=$ax_flag; break])
   done
   FCFLAGS=$ax_save_FCFLAGS
   AC_LANG_POP([Fortran])
])
 case $ax_cv_fc_std2003 in
   error|unknown)
     ifelse([$2],,[AC_MSG_ERROR([cannot find compiler option to support Fortran 2003 standard ISO_C_BINDING])],[$2])
     ;;
   *)
     if test "x$ax_cv_fc_std2003" != xnone; then
       FCFLAGS="$FCFLAGS $ax_cv_fc_std2003"
     fi
     $1
     ;;
 esac
])
