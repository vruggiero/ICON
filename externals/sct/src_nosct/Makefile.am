AM_FCFLAGS = $(FC_MODINC). $(FCFLAGS_f90)

libdir = "${exec_prefix}_nosct/lib"
includedir = "${exec_prefix}_nosct/include"

lib_LTLIBRARIES = libsct.la

if FORTRAN_INTERFACE
libsct_la_SOURCES = nosct.f90
if FC_MODUPPERCASE
sct_mod = SCT.$(FC_MODEXT)
else
sct_mod = sct.$(FC_MODEXT)
endif
nodist_include_HEADERS = $(sct_mod)
else
libsct_la_SOURCES =
nodist_include_HEADERS =
endif

$(sct_mod): nosct.lo

# Clean Fortran module files and Cray's postprocessed files.
CLEANFILES = $(sct_mod) sct.i

