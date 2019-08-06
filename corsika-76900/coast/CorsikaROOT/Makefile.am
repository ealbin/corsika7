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

AM_CXXFLAGS += -I$(top_builddir)/coast/CorsikaIntern
# --------------------------------------   END OF INTERNAL COAST BUILD ------------------------------------------------------------------
else
AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern

endif

###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -O3

###########################################################################
# crsRead
###########################################################################

crsIO_Headers =
noinst_crsIO_Headers =

crsIO:
	mkdir -p $@

crsIO-stamp: $(crsIO_Headers) $(noinst_crsIO_Headers) Makefile
	for i in $^; do \
	  ln -fs @abs_srcdir@/`basename $$i` crsIO/`basename $$i`; \
	done
	touch $@

crsIO-clean:
	rm -f crsIO-stamp
	rm -rf crsIO
	rm -f DirectROOT_Dict.*


DirectROOT_Dict.cc: $(crsIO_Headers) DirectROOTLinkDef.h
	rootcint -f $@ -c $^


crsIO_Headers +=  TCherenkov.h TLongitudinal.h TParticle.h \
		  TRun.h TShower.h

BUILT_SOURCES = crsIO crsIO-stamp

if COAST_INTERNAL 

lib_LIBRARIES = libCorsikaROOT.a
libCorsikaROOT_a_SOURCES = TCherenkov.cc TLongitudinal.cc TParticle.cc \
			TRun.cc TShower.cc
nodist_libCorsikaROOT_a_SOURCES = DirectROOT_Dict.cc
libCorsikaROOT_a_CPPFLAGS = @ROOT_CXXFLAGS@ #-fPIC 
#libCorsikaROOT_a_LDFLAGS = --no-undefined -fPIC --shared
#libCorsikaROOT_a_LIBADD = @ROOT_LIBS@

else

lib_LTLIBRARIES = libCorsikaROOT.la
libCorsikaROOT_la_SOURCES = TCherenkov.cc TLongitudinal.cc TParticle.cc \
			TRun.cc TShower.cc
nodist_libCorsikaROOT_la_SOURCES = DirectROOT_Dict.cc
libCorsikaROOT_la_CPPFLAGS = -fPIC @ROOT_CXXFLAGS@
libCorsikaROOT_la_LDFLAGS = -fPIC --shared --no-undefined
libCorsikaROOT_la_LIBADD = @ROOT_LIBS@

endif

EXTRA_DIST += $(crsIO_Headers) DirectROOTLinkDef.h

clean-local: crsIO-clean

#distclean-local: crsIO-clean
#	rm -f *.a ; rm -f *.la ; rm -f *.so*

#install : $(BUILT_SOURCES) install-am

dist-hook :
	rm -f $(distdir)/DirectROOT_Dict.*

nobase_nodist_include_HEADERS = $(addprefix crsIO/, $(notdir $(crsIO_Headers)))

