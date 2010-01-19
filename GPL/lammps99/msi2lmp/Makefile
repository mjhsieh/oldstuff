TARGET   = msi2lmp2
SRCS     = Msi2LMP2.c \
	   FrcMenu.c\
           ReadCarFile.c\
           ReadMdfFile.c\
	   LinkConnections.c\
	   MakeLists.c \
	   SearchDynamicFrc.c \
	   ReadCoefficientsF1.c ReadCoefficientsF2.c \
	   BuildLMP.c \
           ReadFrcFile.c \
           InitializeItems.c \
           SearchAndFill.c \
	   CreateOutputReport.c
OBJS     = $(SRCS:.c=.o)
HEADERS  = Msi2LMP2.h Forcefield.h
CC       = cc
CFLAGS   = -O 
FRCFILE  = cvff.frc 
FRCFILE2 = cff91.frc 
README   = README
MKFILE   = Makefile
STAMP    = `csh -c  date | sed 's, ,-,g' > stampfile`
VERSION  = `cat stampfile`
NOSTAMP  = `rm stampfile`

$(TARGET) : $(OBJS) 
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm 


$(OBJS) : $(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET) 

print:
	genscript -f Courier6 -2r $(HEADERS) $(SRCS) 

tar:
	tar cvf $(TARGET).tar $(MKFILE) $(README) $(SRCS) $(HEADERS) $(FRCFILE) $(FRCFILE2)
	gzip $(TARGET).tar
	
encode:
	uuencode $(TARGET).tar.gz $(TARGET).tar.gz > $(TARGET).tar.gz.uu

version:
	$(STAMP)
	mkdir $(VERSION)
	cp $(SRCS)     $(VERSION)
	cp $(HEADERS)  $(VERSION)
	cp $(FRCFILE)  $(VERSION)
	cp $(FRCFILE2) $(VERSION)
	cp $(README)   $(VERSION)
	cp $(MKFILE)   $(VERSION)
	$(NOSTAMP)



