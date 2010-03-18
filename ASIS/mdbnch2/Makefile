###################################################################
#	XLF (Mac OS X)
###################################################################
F77 = f77
#FFLAGS= -O0
#FFLAGS= -O2 -qmaxmem=-1
#FFLAGS= -O2 -qmaxmem=-1 -qhot -qextname=etime -qunroll=yes
FFLAGS= -O2 -qmaxmem=-1 -qhot -qextname=etime

##################################################################
# gfortran (Mac OS X)

mdbnch.exe: mdbnch.o second.o
	$(F77) $(FFLAGS) -lm mdbnch.o second.o \
	-o mdbnch.exe
.f.o:
	$(F77) $(FFLAGS)  -c $*.f

clean:
	rm -f *.o *.a mdbnch.exe
