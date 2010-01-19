# Multiple-machine Makefile

SHELL = /bin/sh
.IGNORE:

# Files

SRCF =		lammps.f communicate.f diagnostic.f dump_atom.f dump_vel.f \
		dump_force.f ewald.f ewald_coeff.f \
		finish.f fix.f force.f force_bond.f force_class2.f \
		force_many.f force_respa.f initialize.f input.f \
		integrate.f integrate_respa.f lapack.f mesh_3d.f \
		minimize.f min_support.f misc.f neighbor.f \
		pppm.f pppm_coeff.f pppm_remap.f \
		random.f read_data.f read_restart.f setup.f \
		setup_special.f start.f stats.f string.f thermo.f \
		velocity.f write_restart.f

SRCC = 		fft_3d.c

INC = 		lammps.h param.h

# Definitions

ROOT =		lmp
EXE =		$(ROOT)_$@
SRC = 		$(SRCF) $(SRCC)
OBJ = 		$(SRCF:.f=.o) $(SRCC:.c=.o)

# Help

help:
	@echo 'Type "make target" where target is one of:'
	@echo '      serial	(for SGI O2/10000)'
	@echo '      sgi	(for SGI O2/10000 w/ MPI)'
	@echo '      pgon	(for Intel Paragon)'
	@echo '      tflop	(for Intel Tflop)'
	@echo '      cplant	(for Computational Plant)'
	@echo '      origin	(for SGI Origin)'
	@echo '      t3e	(for Cray T3E)'
	@echo '      8400	(for Dec 8400 Cluster)'
	@echo '      geo	(for Dec geo)'
	@echo '      c90	(for Cray C90)'

# Targets

serial sgi pgon tflop cplant origin t3e 8400 geo c90:
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) $(INC) Makefile*; fi

# Cleans

clean_serial:
	rm -r Obj_serial

clean_sgi:
	rm -r Obj_sgi

clean_pgon:
	rm -r Obj_pgon

clean_cplant:
	rm -r Obj_cplant

clean_tflop:
	rm -r Obj_tflop

clean_origin:
	rm -r Obj_origin

clean_t3e:
	rm -r Obj_t3e

clean_8400:
	rm -r Obj_8400

clean_geo:
	rm -r Obj_geo

clean_c90:
	rm -r Obj_c90
