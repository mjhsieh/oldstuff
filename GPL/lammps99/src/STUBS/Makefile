# Makefile for MPI stubs

SHELL = /bin/sh
.IGNORE:

# Files

SRCF =		mpi_fortran.f
SRCC =		mpi_c.c
INC =		mpif.h mpi.h

# Definitions

EXE =		libmpi.a
OBJ = 		$(SRCF:.f=.o) $(SRCC:.c=.o)

# System-specific settings

F77 =		f77
F77FLAGS =	-O -n32
CC =		cc
CCFLAGS =	-O -n32
ARCHIVE =	ar
ARCHFLAG =	r

# Target

$(EXE):	$(OBJ)
	$(ARCHIVE) $(ARCHFLAG) $(EXE) $(OBJ)

# Clean

clean:
	rm *.o libmpi.a

# Compilation rules

.f.o:
	$(F77) $(F77FLAGS) -c $<

.c.o:
	$(CC) $(CCFLAGS) -c $<

# Individual dependencies

$(OBJ):	$(INC)

