#!/bin/bash
#
# job_submit -p1 -cp -t360 -m1000 readcskralplot.sh
# 
# create file list and run `readcskralplot` program and
# sum up all particle data files to histograms (Ralphs version):
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk001234/
# 
if [ -e DAT001234-000001 ] ; then
  ls -1 DAT001234-?????? > readcskralplot.i001234
else
  ls -1 DAT001234-* | grep t -v | grep l -v | grep a -v > readcskralplot.i001234
fi
# 
# names of sub paths csk00????;
# gfortran -fbounds-check readcskralplot.f -o readcskralplot
# ifort -C -check bounds readcskralplot.f -o readcskralplot
#
  ./readcskralplot < readcskralplot.i001234 > readcskralplot.out001234
  mv fort.19 readcskralplot.fort001234
