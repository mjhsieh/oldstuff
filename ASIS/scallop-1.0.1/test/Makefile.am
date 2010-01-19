## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = poisson_test

INCLUDES = -DDIM3 -DDEFINE_TEMPLATES @FFTW_INCLUDE_DIR@ @KELP_INCLUDE_DIR@ @MPI_INCLUDE_DIR@ -I$(top_srcdir)/kAPI -I$(top_srcdir)/src

FFLAGS += @FFLAGS@ @OTHER_FFLAGS@
CXXFLAGS += @CXXFLAGS@ @OTHER_CXXFLAGS@

poisson_test_SOURCES = \
	poisson_test.h \
	poisson_test.C

poisson_test_LDADD = \
	$(top_srcdir)/src/libpoisson.a \
	$(top_srcdir)/kAPI/libkAPI.a \
	@FFTW_LD_DIR@ @FFTW_LD_LIBS@ \
	@KELP_LD_DIR@ @KELP_LD_LIBS@ \
	@MPI_LD_DIR@ @MPI_LD_LIBS@ \
	@NETWORKING_LIBS@ \
	@FLIBS@

CLEANFILES = *.rpo
