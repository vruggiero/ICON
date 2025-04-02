AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = -I$(top_builddir)/src
AM_FCFLAGS = $(FC_MODINC). $(FC_MODINC)$(top_builddir)/src $(FCFLAGS_f90)
if MPI
AM_FCFLAGS += $(FC_DEFINE)HAVE_MPI
endif
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_INCLUDES) $(FCCPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

LDADD = $(top_builddir)/src/libsct.la

check_PROGRAMS =

if MPI
if MPICHECK
check_PROGRAMS += test_precision
endif
else
check_PROGRAMS += test_precision
endif

if FORTRAN_INTERFACE
if MPI
if MPICHECK
check_PROGRAMS += test_sct
endif
else
check_PROGRAMS += test_sct
endif
endif

if HDF5
if MPI
if MPICHECK
check_PROGRAMS += test_hdf5
check_PROGRAMS += test_hdf5_fortran
endif
else
check_PROGRAMS += test_hdf5
check_PROGRAMS += test_hdf5_fortran
endif
LDADD += $(LTLIBHDF5)
endif

if PAPI
if MPI
if MPICHECK
check_PROGRAMS += test_papi
endif
else
check_PROGRAMS += test_papi
endif
LDADD += $(LTLIBPAPI)
endif

if !MPI
if !OPENMP
check_PROGRAMS += test_serial
endif
if OPENMP
check_PROGRAMS += test_omp
endif
endif

if MPI
if MPICHECK
if !OPENMP
check_PROGRAMS += test_mpi
endif
if OPENMP
check_PROGRAMS += test_hybrid
endif
endif
endif

test_papi_SOURCES = test_papi.c
test_hdf5_SOURCES = test_hdf5.c
test_hdf5_fortran_SOURCES = test_hdf5_fortran.f90
test_serial_SOURCES = test_serial.c
test_omp_SOURCES = test_omp.c
test_mpi_SOURCES = test_mpi.c
test_hybrid_SOURCES = test_hybrid.c
test_precision_SOURCES = test_precision.c
test_sct_SOURCES = test_sct.f90

AM_COLOR_TESTS=yes

TESTS =

if MPI
if MPICHECK
TESTS += test_precision_run
endif
else
TESTS += test_precision_run
endif

if FORTRAN_INTERFACE
if MPI
if MPICHECK
TESTS += test_sct_run
endif
else
TESTS += test_sct_run
endif
endif

if PAPI
TESTS += test_papi_run
endif

if HDF5
if MPI
if MPICHECK
TESTS += test_hdf5_run
TESTS += test_hdf5_fortran_run
endif
else
TESTS += test_hdf5_run
TESTS += test_hdf5_fortran_run
endif
endif

if !MPI
if !OPENMP
TESTS += test_serial_run
endif
if OPENMP
TESTS += test_omp_run
endif
endif

if MPI
if MPICHECK
if !OPENMP
TESTS += test_mpi_run
endif
if OPENMP
TESTS += test_hybrid_run
endif
endif
endif

# Remove module files, Cray's postprocessed files and hdf5 output upon make clean:
clean-local:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
	$(RM) *.i *.h5 *.h5_????-??-??_??:??:??
	$(RM) -rf *.dSYM
