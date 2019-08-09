###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.1.1.1 2007-07-31 07:00:33 rulrich Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern \
               -I$(top_builddir)/CorsikaToROOT \
               -I$(top_builddir)/CorsikaFileIO \
               -I$(top_builddir)/CorsikaInterface

###########################################################################
# rootout
###########################################################################

rootout_Headers =
noinst_rootout_Headers =

#BUILT_SOURCES = rootout rootout-stamp
lib_LTLIBRARIES = libCorsikaROOTOut.la libMachineIndependent.la

libCorsikaROOTOut_la_SOURCES = CorsikaROOTOut.cc
libCorsikaROOTOut_la_CPPFLAGS = -fPIC 
libCorsikaROOTOut_la_LDFLAGS = -fPIC --shared

libMachineIndependent_la_SOURCES = MachineIndependent.cc
libMachineIndependent_la_CPPFLAGS = -fPIC 
libMachineIndependent_la_LDFLAGS = -fPIC --shared

#EXTRA_DIST += $(rootout_Headers) 
#clean-local: rootout-clean
#install : $(BUILT_SOURCES) install-am
#nobase_dist_include_HEADERS += $(addprefix rootout/, $(notdir $(rootout_Headers)))
#EXTRA_DIST += $(rootout_Headers)
