## Makefile.am -- Process this file with automake to get Makefile.in

srcdir = @srcdir@
VPATH = @srcdir@

lib_LIBRARIES = liburqmd.a
liburqmd_a_FFLAGS = @URQMD_FFLAGS@
liburqmd_a_SOURCES = addpart.f angdis.f anndec.f blockres.f boxprg.f \
			cascinit.f coload.f dectim.f delpart.f detbal.f \
	   		dwidth.f error.f getmass.f getspin.f init.f  \
	   		iso.f ityp2pdg.f jdecay2.f make22.f numrec.f  \
	   		output.f paulibl.f proppot.f saveinfo.f scatter.f \
	   		siglookup.f string.f tabinit.f urqmd.f whichres.f 

EXTRA_DIST = README Copyright GNUmakefile_corsika \
	boxinc.f colltab.f comnorm.f comres.f coms.f \
	comstr.f comwid.f freezeout.f inputs.f newpart.f \
	options.f outcom.f

README:
	@if test -f urqmd1.3_cors.tar.gz  ; then \
	tar zxf urqmd1.3_cors.tar.gz ; rm -f urqmd1.3_cors.tar.gz ;\
	touch decomp.flag ; \
	fi

dist-hook:
	if test -f decomp.flag  ; then \
	echo "Restore compression of urqmd package..." ; \
	tar --remove-file -zcf urqmd1.3_cors.tar.gz *.f Copyright \
	GNUmakefile_corsika README ; \
	rm -f decomp.flag ; \
	fi
