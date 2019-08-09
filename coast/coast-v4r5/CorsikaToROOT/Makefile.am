###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.1.1.1 2007-07-31 06:58:59 rulrich Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern \
               -I$(top_builddir)/CorsikaROOT \
               -I$(top_builddir)/CorsikaFileIO

###########################################################################
# crsRead
###########################################################################

crs2r_Headers =
noinst_crs2r_Headers =

crs2r:
	mkdir -p $@

crs2r-stamp: $(crs2r_Headers) 
	for i in $^; do \
	  ln -fs @abs_srcdir@/`basename $$i` crs2r/`basename $$i`; \
	done
	touch $@

crs2r-clean:
	rm -f crs2r-stamp
	rm -rf crs2r

crs2r_Headers += MCorsikaReader.h TC2R.h

BUILT_SOURCES = crs2r crs2r-stamp
lib_LTLIBRARIES = libCorsikaToROOT.la
libCorsikaToROOT_la_SOURCES = MCorsikaReader.cc TC2R.cc
libCorsikaToROOT_la_CPPFLAGS = -fPIC @ROOT_CXXFLAGS@
libCorsikaToROOT_la_LDFLAGS = -fPIC --shared --no-undefined
libCorsikaToROOT_la_LIBADD = @ROOT_LIBS@ 

bin_PROGRAMS = corsika2root
corsika2root_SOURCES = corsika2root.cc
corsika2root_CPPFLAGS = @ROOT_CXXFLAGS@
corsika2root_LDFLAGS = -Wl,--no-as-needed \
	               -L$(top_builddir)/CorsikaIntern -lCorsikaIntern \
                       -L$(top_builddir)/CorsikaROOT -lCorsikaROOT \
                       -L$(top_builddir)/CorsikaFileIO -lCorsikaFileIO \
	               -lCorsikaToROOT
corsika2root_LDADD = @ROOT_LIBS@  

EXTRA_DIST += $(crs2r_Headers) 

clean-local: crs2r-clean

#install : $(BUILT_SOURCES) install-am

nobase_dist_include_HEADERS += $(addprefix crs2r/, $(notdir $(crs2r_Headers)))
EXTRA_DIST += $(crs2r_Headers)

