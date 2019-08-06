###########################################################################
# Makefile template
# $Id: Makefile.am 5118 2016-01-04 19:14:11Z darko $
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

AM_CXXFLAGS = 
nobase_dist_include_HEADERS = 
EXTRA_DIST =



###########################################################################
# globals
###########################################################################


AM_CXXFLAGS += -O3

###########################################################################
# rootout
###########################################################################

rootout_Headers =
noinst_rootout_Headers =

#BUILT_SOURCES = rootout rootout-stamp
if COAST_INTERNAL

AM_CXXFLAGS += -I$(top_builddir)/coast/CorsikaIntern \
               -I$(top_builddir)/coast/CorsikaToROOT \
               -I$(top_builddir)/coast/CorsikaFileIO \
               -I$(top_builddir)/coast/CorsikaInterface

lib_LIBRARIES = libCorsikaROOTOut.a libMachineIndependent.a

libCorsikaROOTOut_a_SOURCES = CorsikaROOTOut.cc
libCorsikaROOTOut_a_CPPFLAGS = #-fPIC 
#libCorsikaROOTOut_a_LDFLAGS = -fPIC --shared

libMachineIndependent_a_SOURCES = MachineIndependent.cc
libMachineIndependent_a_CPPFLAGS = #-fPIC 
#libMachineIndependent_a_LDFLAGS = -fPIC --shared

else

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern \
               -I$(top_builddir)/CorsikaToROOT \
               -I$(top_builddir)/CorsikaFileIO \
               -I$(top_builddir)/CorsikaInterface

lib_LTLIBRARIES = libCorsikaROOTOut.la libMachineIndependent.la

libCorsikaROOTOut_la_SOURCES = CorsikaROOTOut.cc
libCorsikaROOTOut_la_CPPFLAGS = -fPIC 
libCorsikaROOTOut_la_LDFLAGS = -fPIC --shared

libMachineIndependent_la_SOURCES = MachineIndependent.cc
libMachineIndependent_la_CPPFLAGS = -fPIC 
libMachineIndependent_la_LDFLAGS = -fPIC --shared

endif

#EXTRA_DIST += $(rootout_Headers) 
#clean-local: rootout-clean
#install : $(BUILT_SOURCES) install-am
#nobase_dist_include_HEADERS += $(addprefix rootout/, $(notdir $(rootout_Headers)))
#EXTRA_DIST += $(rootout_Headers)
