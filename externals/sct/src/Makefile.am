AM_CFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(FC_MODINC).. $(FCFLAGS_f90) $(FCDEFS)

lib_LTLIBRARIES = libsct.la
libsct_la_LIBADD = $(LTLIBPAPI) $(LTLIBHDF5)

libsct_la_LDFLAGS = $(OPENMP_FC_C_LT_FLAGS) -version-info @LIB_CURRENT@:@LIB_REVISION@:@LIB_AGE@
header = $(srcdir)/sct_config.h \
	$(srcdir)/sct_mach.h \
	$(srcdir)/sct_collector.h \
	$(srcdir)/sct_reduce.h \
	$(srcdir)/sct_reporter.h \
	$(srcdir)/sct_mergesort.h

nodist_include_HEADERS = sct.h

BUILT_SOURCES = sct.h

libsct_la_SOURCES = sct_collector.c sct_f2c.c sct_mach.c sct_reduce.c sct_reporter.c sct_mergesort.c
if FORTRAN_INTERFACE
libsct_la_SOURCES += sct.f90
if FC_MODUPPERCASE
sct_mod = SCT.$(FC_MODEXT)
else
sct_mod = sct.$(FC_MODEXT)
endif
nodist_include_HEADERS += $(sct_mod)
endif

CLEANFILES = sct.h

# Clean Fortran module files and Cray's postprocessed files.
CLEANFILES += $(sct_mod) sct.i

EXTRA_DIST = $(header)

# Build public header sct.h
sct.h: $(header)
if MPI
	$(top_srcdir)/util/gen_public_header.pl HAVE_MPI $(header)
else
	$(top_srcdir)/util/gen_public_header.pl N_HAVE_MPI $(header)
endif

$(sct_mod): sct.lo

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pkgconfig/sct.pc
