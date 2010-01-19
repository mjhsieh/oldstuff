dist_srcs = README assign_parse_radch.pl parse_charges parse_radii \
            assign_charmm_radch.pl sanity.pl \
            compare-multiflex-results.pl relto_with.pl tiplot.pl \
            test.pkint test.g 

dist: $(dist_srcs)
	@if [ x$(DISTDIR) = x. -o x$(DISTDIR) = x ] ; then  \
          echo "setting DISTDIR=$(DISTDIR) (current dir) too dangerous.";\
          [ "" ] ; fi
	@if [ -d $(DISTDIR) -o -f $(DISTDIR) ] ; then \
          echo $(DISTDIR) already exists ; \
          [ "" ] ; fi
	mkdir $(DISTDIR)
	ln $(dist_srcs) $(DISTDIR)
