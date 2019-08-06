#!/bin/bash
#
# sortaugerhit.sh:
# ================
# create file list and run `sortaugerhit` program;
# script must contain the current run number of a `csk??????` subpath;
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
# ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
#               (cpu time about 4300 minutes for 20 cores and 1 TBy space)
# ------------------------------------------------------------------------
# usage: ./sortaugerhit.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk002810/ 
# 
if [ -e DAT002810 ] ; then
  ls -1 DAT002810 > sortaugerhit.i002810
else
  if [ -e DAT002810-000001 ] ; then
    ls -1 DAT002810-00???? > sortaugerhit.i002810
  else
    ls -1 DAT002810-?????????-????????? > sortaugerhit.i002810
  # ls -1 DAT002810-* | grep t -v | grep l -v | grep a -v > sortaugerhit.i002810
  fi
fi
# run program sortaugerhit:
./sortaugerhit < sortaugerhit.i002810 > sortaugerhit.out002810
# rename text file `fort.20` with all particles of core position 20:
mv fort.20 DAT002810.augerhit20
