ACLOCAL_AMFLAGS = -I m4

SUBDIRS = . src doc tests

if FAKELIB
SUBDIRS += src_nosct
endif

include aminclude.am

EXTRA_DIST = util tests/test_base.pl

doc/html: doc

doc:
if HAVE_DOXYGEN
	@echo "-----------------------------------------"
	@echo "Using doxygen to make doc."
	@echo "----------------------------------------"
	cd src && $(MAKE) sct.h
	$(MAKE) doxygen-doc
else
	@echo "-----------------------------------------"
	@echo "Using make doc requires doxygen to be"
	@echo "installed and found by configure."
	@echo "Check 'config.log' for more information."
	@echo "----------------------------------------"
	$(MKDIR_P) $(top_srcdir)/doc/html
endif
