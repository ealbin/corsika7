#!/bin/bash
#
# concatcorsika.sh:
# =================
# create list of (small) particle data files; then run fortran
# program `concatcorsika` to create one big particle data file
# of all small ones (of a parallel corsika simulation,
# including a lot of `zero` particles); must be csk subdirectory.
# ------------------------------------------------------------------------
#        gfortran -O0 -fbounds-check concatcorsika.f -o concatcorsika
#        ifort -C -O0 -check bounds concatcorsika.f -o concatcorsika
# ------------------------------------------------------------------------
# usage: ./concatcorsika.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk072066/
#
# - - - - test login node and concatcorsika code:
computype="ikp"
if [ ! -e concatcorsika.f ] ; then
  echo "          ----- missing source file 'concatcorsika.f'!"
  exit 1
fi
# compile and link considering the login node:
if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
  gfortran -O0 -fbounds-check concatcorsika.f -o concatcorsika
elif [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
  ifort -C -O0 -check bounds concatcorsika.f -o concatcorsika
elif [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
  ifort -C -O0 -check bounds concatcorsika.f -o concatcorsika
fi
#
# - - - - test on `csk` subdirectory of a parallel simulation:
pwd > concatcorsika.tmpname
workdir=`cat concatcorsika.tmpname`
strtest="csk"
kbeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
#
# - - - - testing current path on `csk` subdirectory:
if [ $# -eq 0 ] ; then
  # - - - - - - no execution argument given - - - - - - 
  if [ $kbeg -gt 0 ] ; then
    # being in a parallel simulation subdirectory 'cskiiiiii/':
    sed -i "/csk/ s/csk/   /" concatcorsika.tmpname
    runnrpart=`cat concatcorsika.tmpname | awk '{printf("%06d",$2)}'`
    echo "          ----- concatenate particle data in directory csk$runnrpart/"
    if [ -e DAT$runnrpart-000001 ] ; then
      ls -1 DAT$runnrpart-?????? > concatcorsika.i$runnrpart
    elif  [ -e DAT$runnrpart-?????????-000000001 ] ; then
      ls -1 DAT$runnrpart-?????????-????????? > concatcorsika.i$runnrpart
    fi
    # run program concatcorsika:
    ./concatcorsika < concatcorsika.i$runnrpart > concatcorsika.out$runnrpart
    ls -l concatcorsika.*$runnrpart
  else # not within a 'cskiiiiii/' subdirectory:
    echo "          ----- being not within a 'cskiiiiii/' subdirectory;"
  fi
fi
rm -f concatcorsika.tmpname
#
