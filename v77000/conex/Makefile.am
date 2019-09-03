## Makefile.am -- Process this file with automake to get Makefile.in

srcdir = @srcdir@
VPATH = @srcdir@

cxcompfile.f :
	@CXCOMP_FILE@

lib_LIBRARIES = libconex.a
nodist_libconex_a_SOURCES = cxcompfile.f

libconex_a_FFLAGS = @CONEX_FFLAGS@

EXTRA_DIST = README conex.dke20  conex.dkm20  conex.dkz20  conex.incnex  \
		conex.wepo20  conex.p2epo20  conexep.h conex.p2dec20 \
		conex.wqgs20  conex.p2qgs20  conex.p2ghe20 \
		conex.dkg20  conex.dkn20  conex.h conexwei.h \
		conex.wghe20  conex.wsib20  conex.dkl20  conex.dks20  \
		conex.p2sib20  conex.p2IIq20  conex.p2flu20 \
		conex.wflu20  conex.p2urq20  conex.wurq20 \
		conex.ilo20  conex.wIIq20  ===LOG=== \
		conex_*_*.paramin conex_cors.F \
		egs4.dat egs4.out


@CONEX_PARAM@:
	@if [ ! -e $(bindir)/$@.param ] ; then  \
	echo "==conex==> generating run/$@.param from : $@.paramin  " ; \
	cat $@.paramin | sed s%CONEXROOT%$(prefix)%g > $(bindir)/$@.param ;   \
	fi

install-exec-hook: @CONEX_PARAM@


README:
	@if test -f conex-6.4000.tar.gz  ; then \
	tar zxf conex-6.4000.tar.gz ; rm -f conex-6.4000.tar.gz ;\
	touch decomp.flag ; \
	fi

dist-hook:
	@if test -f decomp.flag  ; then \
	echo "Restore compression of conex package..." ; \
	tar --remove-file -zcf conex-6.4000.tar.gz conex* egs4.* ===LOG=== README; \
	rm -f decomp.flag ; \
	fi


clean-local:
	rm -f cxcompfile.f
