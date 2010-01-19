SHELL=/bin/sh
CC=gcc
AR=ar
AR_FLAGS=rc
DEBUGFLAGS= 
CFLAGS=$(DEBUGFLAGS) 
CPPFLAGS = -I.

.PHONY: clean cleandeps

librx.a : rx.o
	-rm -f t$@
	$(AR) $(AR_FLAGS) t$@ rx.o
	mv t$@ $@

rx.o: rx.c rx.h

dist: Makefile.in rx.c rx.h regex.h
	@if [ x$(DISTDIR) = x. -o x$(DISTDIR) = x ] ; then  \
          echo "setting DISTDIR=$(DISTDIR) (current dir) too dangerous.";\
          [ "" ] ; fi
	@if [ -d $(DISTDIR) -o -f $(DISTDIR) ] ; then \
          echo $(DISTDIR) already exists ; \
          [ "" ] ; fi
	mkdir $(DISTDIR)
	ln Makefile.in rx.c rx.h regex.h $(DISTDIR)
clean:
	-rm rx.o librx.a
cleandeps:
