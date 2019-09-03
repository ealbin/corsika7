###########################################################################
# Makefile template
# $Id: Makefile.am 5118 2016-01-04 19:14:11Z darko $
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

AM_CXXFLAGS = 
nobase_dist_include_HEADERS = 
EXTRA_DIST = *.h



###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -I$(top_builddir)/coast/CorsikaIntern \
               -I$(top_builddir)/coast/CorsikaFileIO \
               -I$(top_builddir)/coast/CorsikaInterface

AM_CXXFLAGS += -O3

###########################################################################
# CoREAS
###########################################################################

lib_LIBRARIES = libInclinedPlane.a 

libInclinedPlane_a_SOURCES = InclinedPlane.cc TPlotter.cc

#libInclinedPlane_a_CPPFLAGS = -fPIC 
#libInclinedPlane_a_LDFLAGS = -fPIC --shared

