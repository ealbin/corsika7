## Makefile.am -- Process this file with automake to get Makefile.in

srcdir = @srcdir@
VPATH = @srcdir@

lib_LIBRARIES = libpythia.a

libpythia_a_SOURCES = fhhiggscorr.f fhsetflags.f fhsetpara.f \
                      pdfset.f py1ent.f  py2ent.f py2frm.f py3ent.f \
                      py4ent.f py4frm.f py4jet.f py4jts.f py4jtw.f \
                      py6frm.f pyadsh.f pyalem.f pyalps.f pyangl.f \
                      pyapps.f pybesq.f pybksb.f pyboei.f pybook.f \
                      pycba2.f pycbal.f pycdiv.f pycell.f pychge.f \
                      pycjdc.f pyckbd.f pyclus.f pycmq2.f pycmqr.f \
                      pycomp.f pycrth.f pycsrt.f pyct5l.f pyct5m.f \
                      pycteq.f pycttr.f pydata.f pydcyk.f pydecy.f \
                      pydiff.f pydisg.f pydocu.f pydump.f pyedit.f \
                      pyeevt.f pyeicg.f pyeig4.f pyerrm.f pyevnt.f \
                      pyevnw.f pyevol.f pyevwt.f pyexec.f pyfact.f \
                      pyfcmp.f pyfeyn.f pyfill.f pyfint.f pyfisb.f \
                      pyfowo.f pyfram.f pyfscr.f pygaga.f pygamm.f \
                      pygano.f pygau2.f pygaus.f pygbeh.f pygdir.f \
                      pygfxx.f pyggam.f pygive.f pyglui.f pygrvd.f \
                      pygrvl.f pygrvm.f pygrvs.f pygrvv.f pygrvw.f \
                      pygvmd.f pyh2xx.f pyhepc.f pyhext.f pyhfth.f \
                      pyhggm.f pyhist.f pyi3au.f pyinbm.f pyindf.f \
                      pyinit.f pyinki.f pyinom.f pyinpr.f pyinre.f \
                      pyjmas.f pyjoin.f pyjurf.f pykcut.f pyk.f \
                      pykfdi.f pykfin.f pyklim.f pykmap.f pylamf.f \
                      pyldcm.f pylhef.f pylist.f pylogo.f pymael.f \
                      pymass.f pymaxi.f pymemx.f pymewt.f pymign.f \
                      pymihg.f pymihk.f pymirm.f pymrun.f pymsin.f \
                      pymult.f pyname.f pynjdc.f pynmes.f pynull.f \
                      pyofsh.f pyonia.f pyoper.f pypcmp.f pypdel.f \
                      pypdfl.f pypdfu.f pypdga.f pypdpi.f pypdpo.f \
                      pypdpr.f pyp.f pypile.f pyplot.f pypole.f \
                      pyprep.f pyptdi.f pyptfs.f pyptis.f pyptmi.f \
                      pyqqbh.f pyradk.f pyrand.f pyreco.f pyremn.f \
                      pyresd.f pyrghm.f pyrnm3.f pyrnmq.f pyrobo.f \
                      pyrvch.f pyrvg1.f pyrvg2.f pyrvg3.f pyrvg4.f \
                      pyrvgl.f pyrvgw.f pyrvi1.f pyrvi2.f pyrvi3.f \
                      pyrvne.f pyrvr.f pyrvsb.f pyrvs.f pyrvsf.f \
                      pysave.f pyscat.f pysfdc.f pysgex.f pysghf.f \
                      pysghg.f pysgqc.f pysgsu.f pysgtc.f pysgwz.f \
                      pyshow.f pysigh.f pysimp.f pyslha.f pyspen.f \
                      pysphe.f pyspli.f pysspa.f pystat.f pystbh.f \
                      pystrf.f pysubh.f pysugi.f pytabu.f pytaud.f \
                      pytbbc.f pytbbn.f pytbdy.f pytbhb.f pytbhg.f \
                      pytbhq.f pytbhs.f pytecm.f pytest.f pythag.f \
                      pythrg.f pythru.f pytime.f pytune.f pyupda.f \
                      pyupev.f pyupin.f pyupre.f pyveto.f pywaux.f \
                      pywidt.f pywidx.f pyx2xg.f pyx2xh.f pyx3jt.f \
                      pyx4jt.f pyxdif.f pyxjet.f pyxkfl.f pyxtee.f \
                      pyxtot.f pyxxga.f pyxxz6.f pyzdis.f ssmssm.f \
                      structm.f structp.f sugra.f upevnt.f upinit.f \
                      upveto.f visaje.f chabapar.f chabadif.f \
                      chamepar.f chamedif.f bobapar.f bobadif.f \
                      bomepar.f bomedif.f 
# libpythia_a_FFLAGS = @PYTHIA_FFLAGS@

EXTRA_DIST = README commons60.f pythia6.4_manual.pdf update_notes6.4.10.txt

README:
	@if test -f pythia-6.4.11.tar.gz ; then \
	tar zxf pythia-6.4.11.tar.gz ; rm -f pythia-6.4.11.tar.gz ;\
	fi

dist-hook:
	tar --remove-file -zcf pythia-6.4.11.tar.gz README update* *.f *.pdf  ;

