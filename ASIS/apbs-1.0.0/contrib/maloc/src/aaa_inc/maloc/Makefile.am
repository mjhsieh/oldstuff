### ##########################################################################
### MALOC = < Minimal Abstraction Layer for Object-oriented C >
### Copyright (C) 1994--2008 Michael Holst
###
### This library is free software; you can redistribute it and/or
### modify it under the terms of the GNU Lesser General Public
### License as published by the Free Software Foundation; either
### version 2.1 of the License, or (at your option) any later version.
###
### This library is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
### Lesser General Public License for more details.
###
### You should have received a copy of the GNU Lesser General Public
### License along with this library; if not, write to the Free Software
### Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
###
### rcsid="$Id: Makefile.am,v 1.47 2008/02/08 18:00:09 fetk Exp $"
### ##########################################################################

## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will use the following:
##
##     COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
##     LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
##
## Author:  Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .

VSYS_HEADS   = vsys.h   vio.h vmem.h vnm.h vset.h
VSH_HEADS    = vsh.h
PSH_HEADS    = psh.h   vmp.h vcom.h vmpi.h

if BUILD_EFENCE
EFENCE_HEADS = 
else
EFENCE_HEADS = 
endif

MALOC_HEADS = maloc.h maloc_base.h
OTHER_HEADS = ${VSYS_HEADS} ${VSH_HEADS} ${PSH_HEADS} ${EFENCE_HEADS}

includedir = ${prefix}/include/maloc
include_HEADERS = ${MALOC_HEADS} ${OTHER_HEADS}

distclean-local: clean-local

clean-local:
	-rm -f ${OTHER_HEADS}

