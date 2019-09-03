#!/bin/bash
#
# readcsk2prtcls.sh:
# ==================
# reading a corsika particle data file and write three smaller files
# of selected particles as DATnnnnnn.muons, .elect, .gamma;
# ------------------------------------------------------------------------
# check:
#        gfortran -O0 -fbounds-check readcsk2prtcls.f -o readcsk2prtcls
#        ifort -C -O0 -check bounds readcsk2prtcls.f -o readcsk2prtcls
# ------------------------------------------------------------------------
# usage: ./readcsk2prtcls.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
# - - parallel path or simply one particle data file:
if [ ! -e readcsk2prtcls.f ] ; then
  echo "          ----- missing source file 'readcsk2prtcls.f'!"
  exit 1
fi
# - - compile and link fortran code considering the login node:
if [ -e "readcsk2prtcls.f" ] ; then
  if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
    gfortran -O0 -fbounds-check readcsk2prtcls.f -o readcsk2prtcls
  elif [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
    ifort -C -O0 -check bounds readcsk2prtcls.f -o readcsk2prtcls
  elif [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
    ifort -C -O0 -check bounds readcsk2prtcls.f -o readcsk2prtcls
  fi
fi
if [ $# -eq 1 ]; then
  # - - - - - - - - - execution argument given - - - - - - - - - - - - - -
  partfile=`echo $1`
  if [ -e "$partfile" ] ; then
    # - - - confirm name of existing particle data file:
    echo $partfile > readcsk2prtcls.inp
    # run program with particle data file name:
    ./readcsk2prtcls < readcsk2prtcls.inp 
    # display names of sorted particle data files:
    ls -l $partfile*
    echo "      "
    rm -f readcsk2prtcls.tmp
  else
    # - - - particle data file not found via or in this path:
    echo "         ... given particle data file not found in this path!" 
  fi
  #
else # elif [ $# -eq 0 ]; then # - - - no execution argument given - - - -
  pwd > readcsk2prtcls.tmp
  cskdirect=`cat readcsk2prtcls.tmp`
  chartest="csk"
  # - - test on a parallel simulation subdirectory cskiiiiii:
  kbeg=`awk -v a="$cskdirect" -v b="$chartest" 'BEGIN{print index(a,b)}'`
  if [ $kbeg -gt 0 ] ; then
    echo "      "
    # - - - extract corsika run number from path name:
    sed -i "/csk/ s/\/csk/\/   /" readcsk2prtcls.tmp
    jrunnr=`cat readcsk2prtcls.tmp | awk '{ printf("%06d",$2) }'`
    if [ -e "DAT$jrunnr-000001" ] ; then
       # newer naming of 'parallel' files:
       ls -1 DAT$jrunnr-0????? > readcsk2prtcls.i$jrunnr
    elif  [ -e DAT$jrunnr-?????????-000000001 ] ; then
       # older naming of 'parallel' files:  
       ls -1 DAT$jrunnr-?????????-????????? > readcsk2prtcls.i$jrunnr
    fi
    jfiles=`wc readcsk2prtcls.i$jrunnr | awk '{ printf("%d",$1) }'`
    echo "          ... sorting electron-, gamma-, muon-particles by"
    echo "              reading $jfiles particle data files in directory"
    echo "              $cskdirect"
    rm -f readcsk2prtcls.tmp
    # - - - run fortran program readcsk2prtcls and display names:
    let "jfiles = $jfiles + 1"
    jrunew=`echo $jfiles | awk '{ printf("%06d",$1) }'`
    ./readcsk2prtcls < readcsk2prtcls.i$jrunnr  > readcsk2prtcls.out$jrunnr
    ls -l readcsk2prtcls.*$jrunnr
    ls -l DAT$jrunnr-$jrunew*
    echo "      "
  else
    # - - - not within a parallel simulation subdirectory 'cskiiiiii/':
    echo "         ... being not within a 'cskiiiiii/' subdirectory,"
    echo "         ... or missing a name of a particle data file!"
    echo "             ./readcsk2prtcls.sh DAT987654"
  fi
  #
fi
