## Makefile.am -- Process this file with automake to get Makefile.in

srcdir = @srcdir@
VPATH = @srcdir@

lib_LIBRARIES = libepos.a
libepos_a_SOURCES =  \
epos-bas-lhc.f  epos-con-lhc.f eposm.f epos_j.f \
epos-dky-lhc.f  epos-dro-lhc.f epos-ems-lhc.f  \
epos-fra-lhc.f  epos-ico-lhc.f epos-hnb-lhc.f \
epos-ids-lhc.f  epos-int-lhc.f epos-jps-lhc.f  \
epos-omg-lhc.f  epos-par-lhc.f epos-sem-lhc.f  \
epos-rsh-lhc.f  epos-qsh-lhc.f epos-tim-lhc.f  \
epos-uti-lhc.f  epos-xan-lhc.f epos-xpr-lhc.f  
 
libepos_a_FFLAGS = @EPOS_FFLAGS@

EXTRA_DIST = README epos.inc epos.incems epos.incico epos.incpar epos.incsem \
	     epos.iniev epos.initl epos.ini1b \
	     epos.param Makefile.ka epos.inchy epos.inirj.lhc epos.inics.lhc

README:
	@if test -f epos-lhc-v3400_cors.tar.gz  ; then \
	tar zxf epos-lhc-v3400_cors.tar.gz ; rm -f epos-lhc-v3400_cors.tar.gz ;\
	fi

dist-hook:
	tar --remove-file -zcf epos-lhc-v3400_cors.tar.gz epos*.f epos.inc* epos.ini* epos.param README ; 
