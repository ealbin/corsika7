#!/bin/bash
#
# concatcorsika.sh:
# =================
# create list of (small) particle data files; then run fortran
# program `concatcorsika` to create one big particle data file
# of all small ones (of a parallel corsika simulation,
# including a lot of `zero` particles); must be csk subdirectory.
# ------------------------------------------------------------------------
# CompLink:
#        gfortran -O0 -fbounds-check concatcorsika.f -o concatcorsika
#        ifort -C -O0 -check bounds concatcorsika.f -o concatcorsika
# ------------------------------------------------------------------------
# usage: ./concatcorsika.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk072066/
# cp ../concatcorsika* .
#
# - - - - test on existing code `concatcorsika`:
if [ ! -e concatcorsika.f ] ; then
  echo "          ----- missing source file 'concatcorsika.f',"
  if [ -e "../concatcorsika.f" ] ; then
    cp ../concatcorsika.f .
    echo "                file copied from '../concatcorsika.f'."
  else
    echo "                not found source file '../concatcorsika.f'!"
    exit 1
  fi
fi
# - - - - compile and link code considering the login node:
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
# - - - - testing current path on a 'cskiiiiii/' subdirectory:
if [ $# -eq 0 ] ; then
  if [ $kbeg -gt 0 ] ; then
    # - - - - - being in a parallel simulation subdirectory 'cskiiiiii/':
    sed -i "/csk/ s/csk/   /" concatcorsika.tmpname
    jrunnr=`cat concatcorsika.tmpname | awk '{printf("%06d",$2)}'`
    echo "          ----- concatenate particle data files in csk$jrunnr/"
    if [ -e DAT$jrunnr-000001 ] ; then
      ls -1 DAT$jrunnr-0????? > concatcorsika.i$jrunnr
    elif [ -e DAT$jrunnr-?????????-000000001 ] ; then
      ls -1 DAT$jrunnr-[0-8]????????-????????? > concatcorsika.i$jrunnr
    fi
    # - - - - - run program concatcorsika:
    ./concatcorsika < concatcorsika.i$jrunnr > concatcorsika.out$jrunnr
    ls -l concatcorsika.*$jrunnr
  else
    # - - - - - not within a 'cskiiiiii/' subdirectory:
    echo "          ----- being not within a 'cskiiiiii/' subdirectory;"
    echo "                current path is '$workdir'."
  fi
# - - - - check given argument starting 'concatcorsika.sh':
elif [ $# -eq 1 ] ; then
  concatarg=`echo $1`
  echo $concatarg > concatcorsika.tmptext
  if [ ${#1} -le 6 ] ; then # check length of argument:
    jrunarg=`cat concatcorsika.tmptext | awk '{printf("%06d",$1)}'`
    if [ $kbeg -eq 0 ] ; then 
      echo "          ----- being not within the valid 'csk$jrunarg/' directory!"
      echo "                current path is '$workdir'."
    else
      pwd > concatcorsika.tmpname
      workdir=`cat concatcorsika.tmpname`
      strtest="csk"
      jbeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
      if [ $jbeg -gt 0 ] ; then
        sed -i "/csk/ s/csk/   /" concatcorsika.tmpname
        jrunnr=`cat concatcorsika.tmpname | awk '{printf("%06d",$2)}'`
        if [ $jrunnr -eq $jrunarg ] ; then
          # - - - - - being within a 'cskiiiiii/' subdirectory:
          echo "          ----- concatenate particle data files in csk$jrunnr/"
          if [ -e DAT$jrunnr-000001 ] ; then
            ls -1 DAT$jrunnr-0????? > concatcorsika.i$jrunnr
          elif [ -e DAT$jrunnr-?????????-000000001 ] ; then
            ls -1 DAT$jrunnr-[0-8]????????-????????? > concatcorsika.i$jrunnr
          fi
          # - - - - - run program concatcorsika:
          ./concatcorsika < concatcorsika.i$jrunnr > concatcorsika.out$jrunnr
          ls -l concatcorsika.*$jrunnr
        else
          echo "          ----- being not within the valid 'csk$jrunarg/' directory!"
          echo "                current path is '$workdir'."
        fi
      fi
    fi
  elif [ ${#1} -eq 9 ] ; then
    sed -i "/DAT/ s/DAT/   /" concatcorsika.tmptext
    sed -i "/csk/ s/csk/   /" concatcorsika.tmptext
    jrunarg=`cat concatcorsika.tmptext | awk '{printf("%06d",$1)}'`
    pwd > concatcorsika.tmpname
    workdir=`cat concatcorsika.tmpname`
    strtest="csk"
    jbeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'` 
    if [ $jbeg -gt 0 ] ; then
      sed -i "/csk/ s/csk/   /" concatcorsika.tmpname
      jrunnr=`cat concatcorsika.tmpname | awk '{printf("%06d",$2)}'`
      if [ $jrunnr -eq $jrunarg ] ; then
        # - - - - - being within a 'cskiiiiii/' subdirectory:
        echo "          ----- concatenate particle data files in csk$jrunnr/"
        if [ -e DAT$jrunnr-000001 ] ; then
          ls -1 DAT$jrunnr-0????? > concatcorsika.i$jrunnr
        elif [ -e DAT$jrunnr-?????????-000000001 ] ; then
          ls -1 DAT$jrunnr-[0-8]????????-????????? > concatcorsika.i$jrunnr
        fi
        # - - - - - run program concatcorsika:
        ./concatcorsika < concatcorsika.i$jrunnr > concatcorsika.out$jrunnr
        ls -l concatcorsika.*$jrunnr
      else
        echo "          ----- being not within the valid 'csk$jrunarg/' directory!"
        echo "                current path is '$workdir'."
      fi
    fi
  fi
  rm -f concatcorsika.tmptext
fi
rm -f concatcorsika.tmpname
#
