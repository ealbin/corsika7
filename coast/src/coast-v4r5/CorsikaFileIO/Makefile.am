###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.2 2007-10-19 11:13:40 rulrich Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

###########################################################################
# globals
###########################################################################

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern

###########################################################################
# crsRead
###########################################################################

crsRead_Headers =
noinst_crsRead_Headers =

crsRead:
	mkdir -p $@

crsRead-stamp: $(crsRead_Headers) $(noinst_crsRead_Headers) Makefile
	for i in $^; do \
	  ln -fs @abs_srcdir@/`basename $$i` crsRead/`basename $$i`; \
	done
	touch $@

crsRead-clean:
	rm -f crsRead-stamp
	rm -rf crsRead

crsRead_Headers += TCorsikaReader.h TSubBlockIO.h \
                   MParticleSubBlockOutput.h  \
                   MCorsikaReader.h  \
                   MLongfileReader.h

BUILT_SOURCES = crsRead crsRead-stamp
lib_LTLIBRARIES = libCorsikaFileIO.la
libCorsikaFileIO_la_SOURCES = TCorsikaReader.cc TSubBlockIO.cc \
                              MParticleSubBlockOutput.cc  \
                              MCorsikaReader.cc MLongfileReader.cc

libCorsikaFileIO_la_CPPFLAGS = -fPIC
libCorsikaFileIO_la_LDFLAGS = -fPIC --shared --no-undefined


clean-local: crsRead-clean

#install : $(BUILT_SOURCES) install-am

nobase_dist_include_HEADERS += $(addprefix crsRead/, $(notdir $(crsRead_Headers)))
EXTRA_DIST += $(crsRead_Headers)
