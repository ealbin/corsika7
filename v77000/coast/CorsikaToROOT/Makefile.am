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

AM_CXXFLAGS += -I$(top_builddir)/coast/CorsikaIntern \
               -I$(top_builddir)/coast/CorsikaROOT \
               -I$(top_builddir)/coast/CorsikaFileIO

# --------------------------------------   END OF INTERNAL COAST BUILD ------------------------------------------------------------------
else

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern \
               -I$(top_builddir)/CorsikaROOT \
               -I$(top_builddir)/CorsikaFileIO
endif

###########################################################################
# globals
###########################################################################


AM_CXXFLAGS += -O3

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

# right now exclude coast binaries from corsika internal build
if COAST_INTERNAL

lib_LIBRARIES = libCorsikaToROOT.a
libCorsikaToROOT_a_SOURCES = MCorsikaReader.cc TC2R.cc
libCorsikaToROOT_a_CPPFLAGS = @ROOT_CXXFLAGS@  #-fPIC 

else

lib_LTLIBRARIES = libCorsikaToROOT.la
libCorsikaToROOT_la_SOURCES = MCorsikaReader.cc TC2R.cc
libCorsikaToROOT_la_CPPFLAGS = -fPIC @ROOT_CXXFLAGS@
libCorsikaToROOT_la_LDFLAGS = -fPIC --shared --no-undefined
libCorsikaToROOT_la_LIBADD = @ROOT_LIBS@ 

bin_PROGRAMS = corsika2root
corsika2root_SOURCES = corsika2root.cc
corsika2root_CPPFLAGS = @ROOT_CXXFLAGS@
corsika2root_LDFLAGS = -Wl,--no-as-needed \
	               -L$(top_builddir)/CorsikaIntern \
                       -L$(top_builddir)/CorsikaROOT \
                       -L$(top_builddir)/CorsikaFileIO
corsika2root_LDADD = @COAST_LIBS@ @COAST_ROOT_LIBS@ @ROOT_LIBS@  

endif

EXTRA_DIST += $(crs2r_Headers) 

clean-local: crs2r-clean

#distclean-local: crs2r-clean
#	rm -f *.a ; rm -f *.la ; rm -f *.so*

#install : $(BUILT_SOURCES) install-am

nobase_nodist_include_HEADERS = $(addprefix crs2r/, $(notdir $(crs2r_Headers)))

