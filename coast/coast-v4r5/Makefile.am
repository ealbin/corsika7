###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.5 2007-10-22 13:46:01 pierog Exp $
###########################################################################

include $(top_srcdir)/Makefile.am.common

SUBDIRS = CorsikaIntern CorsikaFileIO CorsikaInterface

# add root stuff, if root is available
if HAVE_ROOT
  SUBDIRS += CorsikaROOT CorsikaToROOT rootout
endif


# docu, if doxygen
if HAVE_DOXYGEN
  SUBDIRS += Documentation
endif

# configure output
nobase_dist_include_HEADERS += COASTconfig.h
EXTRA_DIST += COASTconfig.h Makefile.userlib

install-exec-hook: userlib

clean-local: userlib-clean

userlib:
	$(MAKE) -f Makefile.userlib

userlib-clean:
	$(MAKE) -f Makefile.userlib clean
