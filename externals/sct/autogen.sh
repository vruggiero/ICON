#!/bin/bash

verlte() {
    [  "$1" = "`echo -e "$1\n$2" | sort -V | head -n1`" ]
}

autoreconf -fvi || exit $?

libtool_version=`{ libtool --version 2>/dev/null || glibtool --version 2>/dev/null; } | awk 'NR==1 {print $4}'`

verlte $libtool_version 2.4.2 && {
  # fix m4/libtool.m4 for nagfor sharedflag
  # see https://trac.mpich.org/projects/mpich/ticket/1870 for details
  patch --forward --no-backup-if-mismatch -p0 -s -i maint/libtool.m4.patch

  # The program 'patch' exits with exitcode=1 if the patch has already been applied.
  # Consider this a normal scenario:
  exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode
}

# Fix the problem when creating shared library with libtool and using MPI
# compiler wrapper:
patch --forward --no-backup-if-mismatch -p1 -s -l -i maint/libtool.m4.nag_wrapper.patch
exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode

# Fix the bug with spaces between flags and their arguments:
patch --forward --no-backup-if-mismatch -p1 -s -l -i maint/libtool.m4.arg_spaces.patch
exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode

# Rebuild configure:
autoconf -f || exit $?

# Reset libtool.m4 timestamps to avoid confusing make:
touch -r m4/ltversion.m4 m4/libtool.m4 || exit $?

# Fix the bug related to '-pthread' in the list of inherited flags, which is not
# recognized by nagfor:
patch --forward --no-backup-if-mismatch -p1 -s -l -i maint/ltmain.sh.nag_pthread.patch
exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode

# All went fine since we have not exited before:
exit 0

