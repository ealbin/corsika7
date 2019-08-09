###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

###########################################################################
# globals
###########################################################################

###########################################################################
# crs
###########################################################################

crs_Headers =
noinst_crs_Headers =

crs:
	mkdir -p $@

crs-stamp: $(crs_Headers) $(noinst_crs_Headers) Makefile
	for i in $^; do \
	  ln -fs @abs_srcdir@/`basename $$i` crs/`basename $$i`; \
	done
	touch $@

crs-clean:
	rm -f crs-stamp
	rm -rf crs

crs_Headers +=  CParticle.h CInteraction.h IParticleReadout.h MCherenkov.h \
		MEventEnd.h MEventHeader.h MLongitudinalBlock.h  \
		MMuonProductionInfo.h MParticle.h MParticleBlock.h \
		MRunEnd.h MRunHeader.h TBlock.h TLongitudinal.h \
		TParticleBlockEntry.h TSubBlock.h CorsikaConsts.h \
		CorsikaTypes.h

BUILT_SOURCES = crs crs-stamp
lib_LTLIBRARIES = libCorsikaIntern.la
libCorsikaIntern_la_SOURCES = CParticle.cc CInteraction.cc IParticleReadout.cc MCherenkov.cc \
		MEventEnd.cc MEventHeader.cc MLongitudinalBlock.cc  \
		MMuonProductionInfo.cc MParticle.cc MParticleBlock.cc \
		MRunEnd.cc MRunHeader.cc TBlock.cc TLongitudinal.cc \
		TParticleBlockEntry.cc TSubBlock.cc
libCorsikaIntern_la_CPPFLAGS = -fPIC
libCorsikaIntern_la_LDFLAGS = -fPIC --shared --no-undefined

clean-local: crs-clean

#install : $(BUILT_SOURCES) install-am

nobase_dist_include_HEADERS += $(addprefix crs/, $(notdir $(crs_Headers)))
EXTRA_DIST += $(crs_Headers)

