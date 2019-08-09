###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.1.1.1 2007-07-31 07:00:35 rulrich Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

###########################################################################
# globals
###########################################################################

###########################################################################
# interface
###########################################################################

interface_Headers =
noinst_interface_Headers =


interface:
	mkdir -p $@

interface-stamp: $(interface_Headers) 
	for i in $^; do \
	  ln -fs @abs_srcdir@/`basename $$i` interface/`basename $$i`; \
	done
	touch $@

interface-clean:
	rm -f interface-stamp
	rm -rf interface


interface_Headers += CorsikaInterface.h

BUILT_SOURCES = interface interface-stamp
#lib_LTLIBRARIES = libCorsikaInterface.la 
#libCorsikaInterface_la_SOURCES = CorsikaROOTOut.cc
#libCorsikaInterface_la_CPPFLAGS = -fPIC 
#libCorsikaInterface_la_LDFLAGS = -fPIC --shared

EXTRA_DIST += $(interface_Headers) 

clean-local: interface-clean

#install : $(BUILT_SOURCES) install-am

nobase_dist_include_HEADERS += $(addprefix interface/, $(notdir $(interface_Headers)))
EXTRA_DIST += $(interface_Headers)
