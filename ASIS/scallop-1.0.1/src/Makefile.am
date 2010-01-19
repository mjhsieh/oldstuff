## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libpoisson.a

INCLUDES = -DDIM3 @FFTW_INCLUDE_DIR@ @KELP_INCLUDE_DIR@ @MPI_INCLUDE_DIR@ -I$(top_srcdir)/kAPI

FFLAGS += @FFLAGS@ @OTHER_FFLAGS@
CXXFLAGS += @CXXFLAGS@ @OTHER_CXXFLAGS@

libpoisson_a_SOURCES = \
	exact_numerics.f \
	helpers.C \
	helpers.h \
	misc_numerics.f \
	nodal_multigrid_numerics.f\
	nodal_fftwSolver.C \
	nodal_fftwSolver.h \
	3dpoisson.F \
	3dsine.F \
	fftw_f77.i \
	nodal_infdmn.C \
	nodal_infdmn.h \
	poisson.h \
	poisson.C \
	rhs_maker_numerics.f \
	scallop_macros.h \
	stencil.h \
	stencil.C \
	timer.h

CLEANFILES = *.rpo
