objects = avere.o cpu_intel.o readmole.o i-forcee.o mdee.o mdstepe.o \
       suple.o shejke.o chenc.o restarte.o l-forcee.o unite.o 

   FC = ifc
   FFLAGS = -O3 -w  
   TARGET = mdee

default: $(objects) 
	$(FC) -o $(TARGET) $(objects)
   $(objects) : mdee.h

clean: 
	rm -f $(TARGET) *.o
