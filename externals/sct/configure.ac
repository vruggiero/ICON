#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([SCT], [1.0], [bockelmann@dkrz.de], [], [http://dkrz-source.gitlab-pages.dkrz.de/sct])
LIB_CURRENT=1
LIB_REVISION=0
LIB_AGE=0
AC_SUBST(LIB_CURRENT)
AC_SUBST(LIB_REVISION)
AC_SUBST(LIB_AGE)

# Checks for system.
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_SYSTEM

AM_INIT_AUTOMAKE([1.14 silent-rules no-dependencies])
AM_MAINTAINER_MODE([enable])
AM_EXTRA_RECURSIVE_TARGETS([doc])
AC_CONFIG_SRCDIR([src/sct_reporter.c])

LT_INIT([pic-only])
dnl Clean up after Libtool checks for C compiler:
rm -f a.out

#
# information on the package
#

# MPI
AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi],
                              [compile with MPI support. If none is found, MPI is not used @<:@default: no@:>@])],
              [case "${enableval}" in
                 yes) mpi=true
                      AC_DEFINE([HAVE_MPI], 1, [MPI parallel enabled])
                      ;;
                 no)  mpi=false ;;
                 *) AC_MSG_ERROR([bad value '${enableval}' for --enable-mpi]) ;;
               esac],
            [mpi=false])

# for some systems (e.g. CRAY) mpirun/aprun check is not possible
AC_ARG_ENABLE([mpirun-check],
              [AS_HELP_STRING([--disable-mpirun-check],
                              [disable mpirun checks which might not work on some systems @<:@default: auto@:>@])],
              [case "${enableval}" in
                 yes) with_mpirun_check=true  ;;
                 no)  with_mpirun_check=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --disable-mpirun-check]) ;;
               esac],
              [with_mpirun_check=auto])

# OpenMP
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--enable-openmp],
                              [enable OpenMP parallel timer @<:@default: false@:>@])],
              [case "${enableval}" in
                 yes) openmp=true
                      AC_DEFINE([HAVE_OPENMP], 1, [OpenMP parallel enabled])
                      ;;
                 no)  openmp=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
               esac],
              [openmp=false])

AM_CONDITIONAL([OPENMP],   [test x"$openmp" = xtrue])
AM_CONDITIONAL([MPI],      [test x"$mpi"    = xtrue])

# DEBUG
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debug mode @<:@default: no@:>@])],
              [case "${enableval}" in
                 yes) debug=true
                      AC_DEFINE([DEBUG], 1, [DEBUG mode enabled])
                      ;;
                 no)  debug=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
               esac],
              [debug=false])

# CHECK_TIMER
AC_ARG_ENABLE([check-timer],
              [AS_HELP_STRING([--enable-check-timer],
                              [set in order to check usage of sct (spoils the accuracy of the performance measurements) @<:@default: no@:>@])],
              [case "${enableval}" in
                 yes) checktimer=true
                      AC_DEFINE([CHECK_TIMER], 1, [CHECK_TIMER mode enabled])
                      ;;
                 no)  checktimer=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-check-timer]) ;;
               esac],
              [checktimer=false])

# NESTED_TIMER
AC_ARG_ENABLE([nested-timer],
              [AS_HELP_STRING([--enable-nested-timer],
                              [enable timer nesting in reports (nested timer cannot be stopped while superordinate timer is still active) @<:@default: yes@:>@])],
              [case "${enableval}" in
                 yes) nestedtimer=true  ;;
                 no)  nestedtimer=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-nested-timer]) ;;
               esac],
              [nestedtimer=true
               AC_DEFINE([NESTED_TIMER], 1, [NESTED_TIMER mode enabled])])

# TIMESTAMP_COUNTER
AC_ARG_ENABLE([timestamp-counter],
              [AS_HELP_STRING([--enable-timestamp-counter],
                              [enable timer based on the processor’s time-stamp counter (a 64-bit MSR) if accessible @<:@default: no@:>@])],
              [case "${enableval}" in
                 yes) timestampcounter=true  ;;
                 no)  timestampcounter=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-timestamp-counter]) ;;
               esac],
              [timestampcounter=false
               AC_DEFINE([TIMESTAMP_COUNTER], 0, [TIMESTAMP_COUNTER mode enabled])])

# FAKELIB
AC_ARG_ENABLE([fakelib],
              [AS_HELP_STRING([--enable-fakelib],
                              [build additional fake libsct with same Fortran interface but no backend @<:@default: no@:>@])],
              [case "${enableval}" in
                 yes) fakelib=true
                      AC_DEFINE([FAKELIB], 1, [FAKELIB mode enabled])
                      ;;
                 no)  fakelib=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-fakelib]) ;;
               esac],
              [fakelib=false])
AM_CONDITIONAL([FAKELIB], [test x"$fakelib" = xtrue])


#
# checks for programs
#
AC_PATH_PROG([PERL], [perl])

DX_DOT_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([sct], [doc/Doxyfile], [])

AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DX_DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])

# Checks for C and Fortran compiler - which might be mpi-compiler wrappers

dnl AS_CASE([$host_os], # todo: find better check here
dnl         [aix*], [if test x"$openmp" = xtrue; then   # change to use thread safe compiler version if openmp is used
dnl          CC=${CC-xlc_r}
dnl          FC=${FC-xlf90_r}
dnl          if test x"$mpi" = xtrue; then
dnl            MPICC=${MPICC-mpcc_r}
dnl            MPIFC=${MPIFC-mpxlf90_r}
dnl            CC=${MPICC}
dnl            FC=${MPIFC}
dnl          fi
dnl        fi])

AX_PROG_CC_MPI([test x"$mpi" != xfalse],
               [],
               [if test x"$mpi" = xtrue; then
                  AC_MSG_ERROR([MPI C compiler requested, but couldn't use MPI.])
                fi])

AC_PROG_FC
dnl Clean up after Libtool checks for Fortran compiler:
rm -f a.out
AX_PROG_FC_MPI([test x"$mpi" != xfalse],
               [],
               [if test x"$mpi" = xtrue; then
                  AC_MSG_ERROR([MPI Fortran compiler requested, but couldn't use MPI.])
                fi])

AS_IF([test x"$mpi" != xfalse],
  # do not check for launch if mpirun-check is disabled
  [AS_IF([test "x$with_mpirun_check" = "xfalse"],
         [AC_MSG_WARN([check for MPI_LAUNCH disabled])],
         [AC_PATH_PROGS([MPI_LAUNCH], [poe mpirun mpiexec], [true])
       # set MP_HOSTFILE if we are using IBM PE
       # $ac_cv_path_MPI_LAUNCH can be a command with arguments:
       set dummy $ac_cv_path_MPI_LAUNCH
       if test x"$(basename $[2])" == x"poe"; then
         # create tmp hostfile and export envVar MP_HOSTFILE
         # to enable compiler check for MPI compiler in case of IBM PE
         if test -z "${MP_HOSTFILE}"; then
           hostname=$(hostname)
           if ! grep "$hostname" $HOME/.rhosts 2>/dev/null >&2 ; then
             AC_MSG_ERROR([
  -----------------------------------------------------
  It seems as if you are using IBM PE (poe) for MPI-
  parallel applications. But current hostname id not
  in ~/.rhosts! Stopping...
  Check 'config.log' for more information.
  Please contact $PACKAGE_BUGREPORT for a bugreport.
  -----------------------------------------------------])
           fi
           MP_HOSTFILE=`mktemp -t hostfile.tmpXXXXXX`
           trap '/bin/rm -f $MP_HOSTFILE' EXIT
           export MP_HOSTFILE
           echo "$hostname" >$MP_HOSTFILE
         fi
         # reset MPI_LAUNCH to wrapper script
         MPI_LAUNCH="$(pwd)/util/mpi_launch_poe"
         ac_cv_path_MPI_LAUNCH=$MPI_LAUNCH
       fi
       # check whether MPI environment is working
       AS_IF([test x"$cross_compiling" = xno],
             [AC_MSG_CHECKING([if $MPI_LAUNCH works])
              AC_LINK_IFELSE([AC_LANG_SOURCE([
@%:@include <stdio.h>
@%:@include <stdlib.h>

@%:@include <mpi.h>

@%:@define xmpi(ret)        \\
  do {                      \\
    if (ret != MPI_SUCCESS) \\
      exit(EXIT_FAILURE);   \\
  } while (0)

int
main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  char *numarg = argv@<:@1@:>@;
  int cmdnum = atoi(numarg);
  int procnum = 1;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &procnum, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD));
  xmpi(MPI_Finalize());
  return (procnum == cmdnum)?EXIT_SUCCESS:EXIT_FAILURE;
}
])],
                             [AS_IF([$MPI_LAUNCH -n 4 ./conftest$EXEEXT 4 >&AS_MESSAGE_LOG_FD 2>&1],
                                    [AC_MSG_RESULT([yes])],
                                    [AC_MSG_RESULT([no])
                                     AS_VAR_IF([with_mpirun_check], [auto],
                                               [AC_MSG_WARN([$MPI_LAUNCH doesn't work])
                                                with_mpirun_check=false],
                                               [AC_MSG_FAILURE([$MPI_LAUNCH doesn't work])])])],
                             [AC_MSG_FAILURE([Cannot compile simple MPI program])])
             ])
])])
AM_CONDITIONAL([MPICHECK], [test x"$with_mpirun_check" != xfalse])


AC_LANG_PUSH([Fortran])
ACX_FC_PP_SRCEXT([f90])
AC_LANG_POP([Fortran])
AC_FC_PP_DEFINE
# Due to a bug in AC_FC_PP_DEFINE the value of ac_fc_srcext is set to 'F',
# which means that all the following tests are run not for *.f90, which we
# need, but for *.F files. Therefore, we restore the value ourselves:
ac_fc_srcext=f90

# Set FCDEFS and adjust -D flag if needed.
FCDEFS="$DEFS -DHAVE_FCONFIG_H"
if test x"$FC_DEFINE" != x-D; then
  FCDEFS=`echo " $FCDEFS" | sed "s@ -D@ $FC_DEFINE@g"`
fi
AC_SUBST([FCDEFS])

AC_LANG_PUSH([Fortran])
ACX_FC_LINE_LENGTH([unlimited])
AC_LANG_POP([Fortran])
# Find out how to add to the Fortran include path.
AC_FC_MODULE_FLAG
# Find out the Fortran module file extension.
AC_FC_MODULE_EXTENSION
# Find out whether the names of the Fortran module files are in uppercase.
AX_FC_MODULE_UPPERCASE
AM_CONDITIONAL([FC_MODUPPERCASE], [test x"$FC_MODUPPERCASE" = xyes])
# check for Fortran 2003 standard options
AX_FC_STD2003([fortran_interface=true],
              [AC_MSG_WARN([Fortran compiler does not support ISO_C_BINDING ... suppress building Fortran interface])
               fortran_interface=false])
AM_CONDITIONAL([FORTRAN_INTERFACE], [test x"$fortran_interface" = xtrue])

#
# checks for libraries
#

# check for rtm available:

ac_save_LIBS=$LIBS

# i/ SCT_RTM_READ_REAL_TIME
AC_SEARCH_LIBS([read_real_time], [c])
# sets ac_cv_search_read_real_time to 'none required' if function is already available
#                                  to 'no' if no corresponding lib was found
#                                  '-llibrary' if found

# ii/ SCT_RTM_OMP_GET_WTIME -> omp_get_wtime
# iii/ SCT_RTM_MPI_WTIME	-> MPI_Wtime
# will be set if supported by this built and none else found

# SCT_RTM_POSIX
# iv/ SCT_RTM_CLOCK_GETTIME_MONOTONIC	clock_gettime(CLOCK_MONOTONIC, tm);
# v/  SCT_RTM_CLOCK_GETTIME_REALTIME	clock_gettime(CLOCK_REALTIME, tm);
AC_CHECK_HEADERS([unistd.h time.h])
AC_SEARCH_LIBS([clock_gettime], [rt])
if test x"$ac_cv_search_clock_gettime" != xno; then
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
    [[#include <unistd.h>]],
    [[int n[(_POSIX_TIMERS > 0) ? (+1) : (-1)];]])],
    [ac_cv_test_posix_timers=true],
    [ac_cv_test_posix_timers=false])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
    [[#include <unistd.h>
      #include <time.h>]],
    [[struct timespec tm; clock_gettime(CLOCK_MONOTONIC, &tm);]])],
    [ac_cv_test_posix_monotonic=true],
    [ac_cv_test_posix_monotonic=false])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
    [[#include <unistd.h>
      #include <time.h>]],
    [[struct timespec tm; clock_gettime(CLOCK_REALTIME, &tm);]])],
    [ac_cv_test_posix_realtime=true],
    [ac_cv_test_posix_realtime=false])
fi

# vi/ SCT_RTM_GETTIMEOFDAY
AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_FUNCS([gettimeofday])

# vii/ SCT_RTM_RDTSCP
if test x"$timestampcounter" != xfalse; then
  if test x"$ac_cv_header_inttypes_h" != xno; then
    AC_MSG_CHECKING([if RDTSCP works])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [[#include <inttypes.h>]],
      [[unsigned low, high;
        uint64_t res;
        __asm__ __volatile__ ("rdtscp": "=a" (low), "=d" (high) );
        res = ((low) | ((uint64_t)(high) << 32)) ]])],
      [AC_MSG_RESULT([yes])
       ac_cv_test_rdtscp_timers=true],
      [AC_MSG_RESULT([no])
       ac_cv_test_rdtscp_timers=false])
  else
    ac_cv_test_rdtscp_timers=false
  fi
fi

# restore LIBS and decide which rtm to take
LIBS=$ac_save_LIBS

sct_rtm="none"
if test x"$ac_cv_search_read_real_time" != xno; then
  if test x"$ac_cv_search_read_real_time" != x"none required"; then
    LIBS="$ac_save_LIBS $ac_cv_search_read_real_time"
  fi
  AC_DEFINE([SCT_RTM], [1], [use SCT_RTM_READ_REAL_TIME for time measurement])
  sct_rtm="SCT_RTM_READ_REAL_TIME"
elif test x"$ac_cv_test_rdtscp_timers" = xtrue; then
  AC_DEFINE([SCT_RTM], [64], [use SCT_RTM_RDTSCP for time measurement])
  sct_rtm="SCT_RTM_RDTSCP"
elif test x"$ac_cv_test_posix_timers" = xtrue && test x"$ac_cv_test_posix_monotonic" = xtrue; then
  if test x"$ac_cv_search_clock_gettime" != x"none required"; then
    LIBS="$ac_save_LIBS $ac_cv_search_clock_gettime"
  fi
  AC_DEFINE([SCT_RTM], [8], [use SCT_RTM_CLOCK_GETTIME_MONOTONIC for time measurement])
  sct_rtm="SCT_RTM_CLOCK_GETTIME_MONOTONIC"
elif test x"$ac_cv_test_posix_timers" = xtrue && test x"$ac_cv_test_posix_realtime" = xtrue; then
  if test x"$ac_cv_search_clock_gettime" != x"none required"; then
    LIBS="$ac_save_LIBS $ac_cv_search_clock_gettime"
  fi
  AC_DEFINE([SCT_RTM], [16], [use SCT_RTM_CLOCK_GETTIME_REALTIME for time measurement])
  sct_rtm="SCT_RTM_CLOCK_GETTIME_REALTIME"
elif test x"$openmp" = xtrue; then
  AC_DEFINE([SCT_RTM], [2], [use SCT_RTM_OMP_GET_WTIME for time measurement])
  sct_rtm="SCT_RTM_OMP_GET_WTIME"
elif test x"$mpi" = xtrue; then
  AC_DEFINE([SCT_RTM], [4], [use SCT_RTM_MPI_WTIME for time measurement])
  sct_rtm="SCT_RTM_MPI_WTIME"
elif test x"$ac_cv_header_sys_time_h" = xyes && test x"$ac_cv_func_gettimeofday" = xyes; then
  AC_DEFINE([SCT_RTM], [32], [use SCT_RTM_GETTIMEOFDAY for time measurement])
  sct_rtm="SCT_RTM_GETTIMEOFDAY"
else
  AC_MSG_ERROR([
  -----------------------------------------------------
  Could not find any usable read time method (RTM).
  Stopping...
  Check 'config.log' for more information.
  Please contact $PACKAGE_BUGREPORT for a bugreport.
  -----------------------------------------------------])
fi


# check for PAPI
dnl use gnulib module 'havelib' ... see https://www.gnu.org/software/gnulib/manual/html_node/Searching-for-Libraries.html
dnl needed to download m4 files m4/m4_lib-* and config/config.rpath
AC_ARG_ENABLE([papi],
              [AS_HELP_STRING([--enable-papi],
                              [enable search for PAPI lib @<:@default: false@:>@])],
              [case "${enableval}" in
                 yes) papi=true
                      ;;
                 no)  papi=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-papi]) ;;
               esac],
              [papi=false])

AS_IF([test x"$papi" != xfalse],
      [ AC_LIB_HAVE_LINKFLAGS(papi,,[#include <papi.h>],,
              [
  -----------------------------------------
  Using PAPI interface requires library
  and header files. Skipping ...
  Check 'config.log' for more information,
  maybe libdirstem was determined false.
  -----------------------------------------])
      ], [])
AM_CONDITIONAL([PAPI], [test x"$HAVE_LIBPAPI" = xyes])

# check for HDF5
AC_ARG_ENABLE([hdf5],
              [AS_HELP_STRING([--enable-hdf5],
                              [enable search for HDF5 lib @<:@default: false@:>@])],
              [case "${enableval}" in
                 yes) hdf5=true
                      ;;
                 no)  hdf5=false ;;
                 *) AC_MSG_ERROR([bad value ${enableval} for --enable-hdf5]) ;;
               esac],
              [hdf5=false])

AS_IF([test x"$hdf5" != xfalse],
      [AC_LIB_HAVE_LINKFLAGS(hdf5,,[#include <hdf5.h>])
       AS_VAR_IF([HAVE_LIBHDF5], [no],
                 [AC_MSG_FAILURE([HDF5 library and/or header files are not found])])])
AM_CONDITIONAL([HDF5], [test x"$HAVE_LIBHDF5" = xyes])

AM_COND_IF([HDF5],
  [AC_PATH_PROG([H5DUMP], [h5dump], [unknown])],
  [H5DUMP=unknown])
AM_CONDITIONAL([H5DUMP], [test x"$H5DUMP" != xunknown])

# checks for header files
AC_CHECK_HEADERS([math.h stdio.h limits.h stdlib.h string.h strings.h inttypes.h sys/systemcfg.h])


# checks for typedefs, structures, and compiler characteristics
AC_C_INLINE
AC_C_RESTRICT

if test x"$openmp" = xtrue; then
  if test x"$fortran_interface" = xtrue; then
    AX_FC_C_OPENMP
    if test x"$ax_cv_prog_fc_c_openmp" = xunsupported; then
      AC_MSG_WARN([Fortran compiler cannot link C code that uses OpenMP])
    elif test x"$ax_cv_prog_fc_c_openmp" != "xnone needed"; then
      OPENMP_FC_C_LT_FLAGS=`AS_ECHO(["$OPENMP_FC_C_FLAGS"]) | sed 's%-Wl,@<:@^ @:>@* *%-XCClinker &%g'`
      AC_SUBST([OPENMP_FC_C_LT_FLAGS])
    fi
  else
    AC_LANG_PUSH([C])
    AC_OPENMP
    AC_LANG_POP([C])
  fi
fi


#
# checks for library functions
#
AC_SEARCH_LIBS([sqrt], [m])
AC_CHECK_FUNCS([floor pow sqrt memset strcasecmp])


#
# Output
#
AC_CONFIG_FILES([tests/test_precision_run
                 tests/test_sct_run
                 tests/test_papi_run
                 tests/test_hdf5_run
                 tests/test_hdf5_fortran_run   
                 tests/test_serial_run
                 tests/test_omp_run
                 tests/test_mpi_run
                 tests/test_hybrid_run],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 src/Makefile
                 src/pkgconfig/sct.pc
                 src_nosct/Makefile
                 tests/Makefile])

# Create config headers: config.h comes from autoheader-generated config.h.in,
AC_CONFIG_HEADERS([config.h])

# Convert fconfig.h.in into an fconfig.h suitable for inclusion into Fortran.
# The sed script transforms C style /* comments */ into Fortran style.
AC_CONFIG_HEADERS([fconfig.h:config.h.in],
 [[sed '
    /\/\*\(.*\)\*\//{
      s,,!\1,p
      d
    }
    /\/\*/,/\*\//{
      s,/\*,,
      s,\*/,,
      s,^,!,p
      d
    }' fconfig.h > fconfig.ht
   mv -f fconfig.ht fconfig.h
]])

AC_OUTPUT

echo \
"------------------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix          : ${prefix}
C-Compiler      : ${CC} ${DEFS} ${CFLAGS} ${CPPFLAGS}
Fortran-Compiler: ${FC} ${FCDEFS} ${FCFLAGS} ${FCFLAGS_f90} ${FCCPPFLAGS}


Package features:

MPI              : ${mpi}"

echo \
"OpenMP           : ${openmp}"

if test x"$openmp" = xtrue; then
echo \
"    CFLAGS       : ${OPENMP_CFLAGS}
    FCFLAGS      : ${OPENMP_FC_C_FLAGS}"
fi

echo \
"PAPI library     : ${LTLIBPAPI}"

echo \
"HDF5 library     : ${LTLIBHDF5}"

echo \
"use timer nesting: ${nestedtimer}
use RTM          : ${sct_rtm}

Now type 'make @<:@<target>@:>@' with optional <target>:

   all           - build libsct
   check         - check some features of sct
   install       - install configured version only"

if test -n "$DX_DOXYGEN"; then
echo \
"   doc           - build documentation based on Doxygen"
fi

echo \
"
------------------------------------------------------------"
