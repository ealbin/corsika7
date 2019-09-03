###########################################################################
# Makefile template
# $Id: Makefile.am 5118 2016-01-04 19:14:11Z darko $
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

AM_CXXFLAGS = 
nobase_dist_include_HEADERS = 
EXTRA_DIST =





if COAST_INTERNAL 

# --------------------------------------   END OF INTERNAL COAST BUILD ------------------------------------------------------------------
else

endif

###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -O3

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

nobase_nodist_include_HEADERS = $(addprefix interface/, $(notdir $(interface_Headers)))
