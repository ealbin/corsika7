###########################################################################
# Makefile for documentation
# $Id: Makefile.am 5118 2016-01-04 19:14:11Z darko $
###########################################################################
# Global header
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

#docdir = $(pkgdatadir)/doc
nobase_doc_DATA = 
#nobase_nodist_doc_DATA =
#nobase_dist_doc_DATA =
EXTRA_DIST =

# include $(top_srcdir)/Makefile.am.common


###########################################################################
# globals
###########################################################################

# need code for installation of doxygen output
# latex requires a large poolsize to be able to process the doxygen
# output. Try to find a portable way to adjust the poolsize. In the 
# meanwhile, don't process latex sources.

if HAVE_DOXYGEN
all: COASTDocumentation 
endif



###########################################################################
# Examples
###########################################################################
Examples_DIST_SUBS = 
Examples_DIST =


Examples_NoviceSimple = \
   Examples/novice/simple/CorsikaPlotter.cc \
   Examples/novice/simple/Makefile

Examples_DIST += $(Examples_NoviceSimple)


Examples_NoviceMedium = \
   Examples/novice/medium/Makefile \
   Examples/novice/medium/ReaderExample.cc \
   Examples/novice/medium/TReaderExample.h \
   Examples/novice/medium/TReaderExample.cc 

Examples_DIST += $(Examples_NoviceMedium)


Examples_NoviceAdvanced = \
   Examples/novice/advanced/CorsikaROOTPlotter.cc \
   Examples/novice/advanced/CorsikaROOTPlotterInterface.cc \
   Examples/novice/advanced/Makefile \
   Examples/novice/advanced/TCorsikaPlotter.cc \
   Examples/novice/advanced/TCorsikaPlotter.h \
   Examples/novice/advanced/TSimpleReader.cc \
   Examples/novice/advanced/TSimpleReader.h

#Examples_NoviceAdvanced_SUBS = \
#   Examples/novice/advanced/Makefile.in

Examples_DIST += $(Examples_NoviceAdvanced)
#Examples_DIST_SUBS += $(Examples_NoviceAdvanced_SUBS)


Examples_NoviceExpert = \
   Examples/novice/expert/CompareAnim.cc \
   Examples/novice/expert/ExpertInterface.cc \
   Examples/novice/expert/Makefile \
   Examples/novice/expert/extractTH2D.C \
   Examples/novice/expert/MakeAnim2D.cc \
   Examples/novice/expert/MakeAnim.cc \
   Examples/novice/expert/TPlotter.cc \
   Examples/novice/expert/TPlotter.h

Examples_DIST += $(Examples_NoviceExpert)


Examples_Master3D = \
   Examples/master/3D/Makefile \
   Examples/master/3D/Interface.cc \
   Examples/master/3D/TPlotter.cc \
   Examples/master/3D/TPlotter.h

Examples_DIST += $(Examples_Master3D)



nobase_doc_DATA += $(Examples_DIST) $(Examples_DIST_SUBS:%.in=%) 
EXTRA_DIST += $(Examples_DIST) $(Examples_DIST_SUBS)



###########################################################################
# doxygen documentation
###########################################################################

Images = 
Documents = 
Doxy_DIST = Doxy_MainPage.dox \
            Doxy_footer.html \
            $(Images) \
            $(Documents)

Doxy_DIST_SUBS = Doxyfile.in

COASTDocumentation: Doxyfile
	doxygen $<
#	cp $(srcdir)/Images/favicon.ico $@/html/
	find $@/html/ -iname \*.html | xargs \
	  perl -pi -e 's/<\/head>/<link href="favicon.ico" rel="shortcut icon" type="image\/x-icon">\n<\/head>/'

COASTDocumentation-ps:
	cd COASTDocumentation/latex; $(MAKE) ps

clean-local:
	rm -rf COASTDocumentation

install-data-local:
	tar cf - COASTDocumentation | \
	(cd $(docdir); tar xf -)

uninstall-local:
	rm -rf $(docdir)/COASTDocumentation

nobase_doc_DATA += $(Doxy_DIST) $(Doxy_DIST_SUBS:%.in=%)
EXTRA_DIST += $(Doxy_DIST) $(Doxy_DIST_SUBS) 

