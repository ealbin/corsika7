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
		CorsikaTypes.h PersistentParticle.h

BUILT_SOURCES = crs crs-stamp
if COAST_INTERNAL 

lib_LIBRARIES = libCorsikaIntern.a
libCorsikaIntern_a_SOURCES = CParticle.cc CInteraction.cc IParticleReadout.cc MCherenkov.cc \
		MEventEnd.cc MEventHeader.cc MLongitudinalBlock.cc  \
		MMuonProductionInfo.cc MParticle.cc MParticleBlock.cc \
		MRunEnd.cc MRunHeader.cc TBlock.cc TLongitudinal.cc \
		TParticleBlockEntry.cc TSubBlock.cc PersistentParticle.cc
libCorsikaIntern_a_CPPFLAGS = #-fPIC
#libCorsikaIntern_a_LDFLAGS = --no-undefined -fPIC --shared 

else

lib_LTLIBRARIES = libCorsikaIntern.la
libCorsikaIntern_la_SOURCES = CParticle.cc CInteraction.cc IParticleReadout.cc MCherenkov.cc \
		MEventEnd.cc MEventHeader.cc MLongitudinalBlock.cc  \
		MMuonProductionInfo.cc MParticle.cc MParticleBlock.cc \
		MRunEnd.cc MRunHeader.cc TBlock.cc TLongitudinal.cc \
		TParticleBlockEntry.cc TSubBlock.cc PersistentParticle.cc
libCorsikaIntern_la_CPPFLAGS = -fPIC
libCorsikaIntern_la_LDFLAGS = -fPIC --shared --no-undefined

endif

clean-local: crs-clean

#install : $(BUILT_SOURCES) install-am

nobase_nodist_include_HEADERS = $(addprefix crs/, $(notdir $(crs_Headers)))
EXTRA_DIST += $(crs_Headers)

