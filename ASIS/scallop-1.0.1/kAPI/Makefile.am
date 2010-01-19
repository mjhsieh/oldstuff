## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libkAPI.a

noinst_SCRIPTS = dimensions.sh

INCLUDES = @KELP_INCLUDE_DIR@ @MPI_INCLUDE_DIR@

FFLAGS += @FFLAGS@ @OTHER_FFLAGS@
CXXFLAGS += @CXXFLAGS@ @OTHER_CXXFLAGS@

BUILT_SOURCES = \
	adder2.C \
	adder3.C \
	adder2.h \
	adder3.h \
	def2.h \
	def3.h \
	init2.C \
	init3.C \
	io2.C \
	io3.C \
	io2.h \
	io3.h \
	p_domain2.C \
	p_domain3.C \
	p_domain2.h \
	p_domain3.h \
	p_grid2.C \
	p_grid3.C \
	p_grid2.h \
	p_grid3.h \
	s_grid2.C \
	s_grid3.C \
	s_grid2.h \
	s_grid3.h \
	subtractor2.C \
	subtractor3.C \
	subtractor2.h \
	subtractor3.h \
	templates2.C \
	templates3.C

adder2.C adder3.C: adder.C dimensions.sh
	./dimensions.sh adder.C

adder2.h adder3.h: adder.h dimensions.sh
	./dimensions.sh adder.h

def2.h def3.h: def.h dimensions.sh
	./dimensions.sh def.h

init2.C init3.C: init.C dimensions.sh
	./dimensions.sh init.C

io2.C io3.C: io.C dimensions.sh
	./dimensions.sh io.C

io2.h io3.h: io.h dimensions.sh
	./dimensions.sh io.h

p_domain2.C p_domain3.C: p_domain.C dimensions.sh
	./dimensions.sh p_domain.C

p_domain2.h p_domain3.h: p_domain.h dimensions.sh
	./dimensions.sh p_domain.h

p_grid2.C p_grid3.C: p_grid.C dimensions.sh
	./dimensions.sh p_grid.C

p_grid2.h p_grid3.h: p_grid.h dimensions.sh
	./dimensions.sh p_grid.h

s_grid2.C s_grid3.C: s_grid.C dimensions.sh
	./dimensions.sh s_grid.C

s_grid2.h s_grid3.h: s_grid.h dimensions.sh
	./dimensions.sh s_grid.h

subtractor2.C subtractor3.C: subtractor.C dimensions.sh
	./dimensions.sh subtractor.C

subtractor2.h subtractor3.h: subtractor.h dimensions.sh
	./dimensions.sh subtractor.h

templates2.C templates3.C: templates.C dimensions.sh def2.h def3.h
	./dimensions.sh templates.C

EXTRA_DIST = \
	adder.C \
	adder.h \
	def.h \
	init.C \
	io.C \
	io.h \
	p_domain.C \
	p_domain.h \
	p_grid.C \
	p_grid.h \
	s_grid.C \
	s_grid.h \
	subtractor.C \
	subtractor.h \
	templates.C \
	dimensions.sh

CLEANFILES = \
	adder2.C \
	adder3.C \
	adder2.h \
	adder3.h \
	def2.h \
	def3.h \
	init2.C \
	init3.C \
	io2.C \
	io3.C \
	io2.h \
	io3.h \
	p_domain2.C \
	p_domain3.C \
	p_domain2.h \
	p_domain3.h \
	p_grid2.C \
	p_grid3.C \
	p_grid2.h \
	p_grid3.h \
	s_grid2.C \
	s_grid3.C \
	s_grid2.h \
	s_grid3.h \
	subtractor2.C \
	subtractor3.C \
	subtractor2.h \
	subtractor3.h \
	templates2.C \
	templates3.C \
	*.rpo

libkAPI_a_SOURCES = \
	adder2.C \
	adder3.C \
	adder2.h \
	adder3.h \
	adder_numerics.f \
	init2.C \
	init3.C \
	io2.C \
	io3.C \
	io2.h \
	io3.h \
	p_domain2.C \
	p_domain3.C \
	p_domain2.h \
	p_domain3.h \
	p_grid2.C \
	p_grid3.C \
	p_grid2.h \
	p_grid3.h \
	s_grid2.C \
	s_grid3.C \
	s_grid2.h \
	s_grid3.h \
	s_grid_numerics.f \
	subtractor2.C \
	subtractor3.C \
	subtractor2.h \
	subtractor3.h \
	subtractor_numerics.f \
	templates2.C \
	templates3.C
