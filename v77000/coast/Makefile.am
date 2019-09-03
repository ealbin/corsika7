###########################################################################
# Makefile template
# $Id: Makefile.am 5118 2016-01-04 19:14:11Z darko $
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

AM_CXXFLAGS = 
nobase_dist_include_HEADERS = 
EXTRA_DIST =

AM_CXXFLAGS += -O3

SUBDIRS = CorsikaIntern CorsikaFileIO CorsikaInterface 
# add root stuff, if root is available
if HAVE_ROOT
  SUBDIRS += CorsikaROOT CorsikaToROOT
endif


if COAST_INTERNAL 

SUBDIRS += CorsikaOptions # this does not hurt in corsika and makes our life easier

# DIST_SUBDIRS = $(SUBDIRS) CorsikaROOT CorsikaToROOT

#distclean-local:
#	rm -f *.a ; rm -f *.la ; rm -f *.so*

#dist-hook:
#	tar --remove-file -zcf bernlohr-1.46.tar.gz [^z]*.c *.h *1.dat *2.dat *3.dat *4.dat *5.dat *6.dat *7.dat *8.dat *9.dat *.ps *.pdf Copyright GNUmakefile.org README  ; 


# --------------------------------------   END OF INTERNAL COAST BUILD ------------------------------------------------------------------
else

#nobase_dist_include_HEADERS += COASTconfig.h
EXTRA_DIST += COASTconfig.h 
EXTRA_DIST += Makefile.userlib

# docu, if doxygen, CORSIKA so far never with doxygen
if HAVE_DOXYGEN
  SUBDIRS += Documentation
endif

# ---------------------------------- END OF STANDALONE BUILD --------------------------
endif 





if !COAST_INTERNAL

install-exec-hook: userlib

clean-local: userlib-clean

userlib:
	$(MAKE) -f Makefile.userlib

userlib-clean:
	$(MAKE) -f Makefile.userlib clean

endif
