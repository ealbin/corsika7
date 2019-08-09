###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.2 2007-08-08 15:07:31 rulrich Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern

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
lib_LTLIBRARIES = libCorsikaROOT.la
libCorsikaROOT_la_SOURCES = TCherenkov.cc TLongitudinal.cc TParticle.cc \
			TRun.cc TShower.cc DirectROOT_Dict.cc
libCorsikaROOT_la_CPPFLAGS = -fPIC @ROOT_CXXFLAGS@
libCorsikaROOT_la_LDFLAGS = -fPIC --shared --no-undefined
libCorsikaROOT_la_LIBADD = @ROOT_LIBS@

EXTRA_DIST += $(crsIO_Headers) DirectROOTLinkDef.h

clean-local: crsIO-clean

#install : $(BUILT_SOURCES) install-am

dist-hook :
	rm -f $(distdir)/DirectROOT_Dict.*

nobase_dist_include_HEADERS += $(addprefix crsIO/, $(notdir $(crsIO_Headers)))
EXTRA_DIST += $(crsIO_Headers)
