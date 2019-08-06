#!/bin/bash
#
# sortaugerhit.sh:
# ================
# create file list and run `sortaugerhit` program;
# current path must be a `csk??????` subpath;
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
# ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
#               (cpu time about 4300 minutes for 20 cores and 1 TBy space)
# ------------------------------------------------------------------------
# usage: ./sortaugerhit.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk002810/ (example) 
# 
pwd > sortaugerhit.tmpname
workdir=`cat sortaugerhit.tmpname`
strtest="csk"
jcsk=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
if [ $jcsk -gt 0 ] ; then
  sed -i "/csk/ s/csk/   /" sortaugerhit.tmpname
  jrunnr=`cat sortaugerhit.tmpname | awk '{printf("%06d",$2)}'`
  rm -f sortaugerhit.tmpname
  if [ -e DAT$jrunnr-000001 ] ; then
    ls -1 DAT$jrunnr-0????? > sortaugerhit.i$jrunnr
  else
    ls -1 DAT$jrunnr-?????????-????????? > sortaugerhit.i$jrunnr
  fi
  if [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then 
    ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit 
  fi
  if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
    gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
  fi
  # run program sortaugerhit:
  ./sortaugerhit < sortaugerhit.i$jrunnr > sortaugerhit.out$jrunnr
  # mv fort.20 DAT$jrunnr.augerhit20
else
  echo "          ----- being not within a 'cskiiiiii' subdirectory."
fi
