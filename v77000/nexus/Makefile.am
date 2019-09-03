## Makefile.am -- Process this file with automake to get Makefile.in

srcdir = @srcdir@
VPATH = @srcdir@

lib_LIBRARIES = libnexus.a
libnexus_a_SOURCES = nexus-bas-3972.f  nexus-chk-3972.f nexus-con-3972.f   \
                     nexus-dky-3972.f   nexus-dro-3972.f nexus-ems-3972.f  \
                     nexus-fra-3972.f   nexus-ico-397.f nexus-icp-397.f   \
                     nexus-ids-3972.f   nexus-int-3972.f nexus-jps-397.f   \
                     nexus-omg-3972.f   nexus-par-3972.f nexus-sem-3972.f   \
                     nexus-sha-3972.f  nexus-sto-397.f nexus-tim-397.f   \
                     nexus-uti-3972.f  nexus-xan-3972.f nexus-xpr-397.f   \
                     nexus-xsa-397.f

libnexus_a_FFLAGS = @NEXUS_FFLAGS@

EXTRA_DIST = README nexus.inc nexus.incems nexus.incico nexus.incpar nexus.incsem \
	     nexus.inics nexus.iniev nexus.inirj nexus.initl nexus.param Makefile.ka

README:
	@if test -f nexus3.97_cors.tar.gz  ; then \
	tar zxf nexus3.97_cors.tar.gz ; rm -f nexus3.97_cors.tar.gz ;\
	touch decomp.flag ; \
	fi

dist-hook:
	if test -f decomp.flag  ; then \
	echo "Restore compression of nexus package..." ; \
	tar --remove-file -zcf nexus3.97_cors.tar.gz nexus*.f nexus.inc* \
	nexus.ini* nexus.param Makefile.ka README ;  \
	rm -f decomp.flag ; \
	fi

