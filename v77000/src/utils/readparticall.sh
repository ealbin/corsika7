#!/bin/bash
#
# readparticall.sh:
# =================
# run particle analysis to count all particles of all data files. 
# ------------------------------------------------------------------------
# check:
# gfortran -O0 -fbounds-check readparticall.f -o readparticall
# ifort -C -O0 -check bounds readparticall.f -o readparticall
# ------------------------------------------------------------------------
# usage: ./readparticall.sh [particledatafilename]
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd corsikapath/
#
# - - parallel path or simply one particle data file:
if [ ! -e readparticall.f ] ; then
  echo "          ----- missing source file 'readparticall.f'!"
  exit 1
fi
# - - compile and link fortran code considering the login node:
if [ -e "readparticall.f" ] ; then
  if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
    gfortran -O0 -fbounds-check readparticall.f -o readparticall
  elif [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
    ifort -C -O0 -check bounds readparticall.f -o readparticall
  elif [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
    ifort -C -O0 -check bounds readparticall.f -o readparticall
  fi
fi
if [ $# -eq 1 ] ; then
  # - - - - - - - - - execution argument given - - - - - - - - - - - - - -
  dirpartfile=`echo $1`
  partfile=`basename $dirpartfile | awk '{printf("%s",$1)}'`
  # extract run number and confirm name of particle data file:
  workfile=`echo $partfile`
  strtest="_"
  ibeg=`awk -v a="$workfile" -v b="$strtest" 'BEGIN{print index(a,b)}'`
  if [ $ibeg -gt 0 ] ; then
    jrunnr=987654
    echo $dirpartfile > readparticall.i$jrunnr
    # take care on files like `prE16M562_41` to run;
    # other primaries of this kind were co, gm, he, fe, si.
  else
    strtest=".part"
    jbeg=`awk -v a="$workfile" -v b="$strtest" 'BEGIN{print index(a,b)}'`
    if [ $jbeg -gt 0 ] ; then
      echo $partfile > readparticall.tmp
      sed -i "/DAT/ s/\./ /" readparticall.tmp
      sed -i "/DAT/ s/DAT/   /" readparticall.tmp
      jrunnr=`cat readparticall.tmp | awk '{printf("%06d",$1)}'`
      echo $dirpartfile > readparticall.i$jrunnr
    else
      strtest="DAT"
      lbeg=`awk -v a="$workfile" -v b="$strtest" 'BEGIN{print index(a,b)}'`
      if [ $lbeg -gt 0 ] ; then
        echo $partfile > readparticall.tmp
        sed -i "/DAT/ s/DAT/   /" readparticall.tmp
        jrunnr=`cat readparticall.tmp | awk '{printf("%06d",$1)}'`
        echo $dirpartfile > readparticall.i$jrunnr
      fi
    fi
  fi
  if [ -e "$dirpartfile" ] ; then
    # run program with particle data file name:
    ./readparticall < readparticall.i$jrunnr > readparticall.out$jrunnr
    cat readparticall.out$jrunnr
    ls -l readparticall.*$jrunnr
    echo "      "
  else
    # - - - particle data file not found via or in this path:
    echo "         ... given particle data file not found in this path!" 
  fi
  #
else # elif [ $# -eq 0 ]; then # - - - no execution argument given:
  pwd > readparticall.tmp
  cskdirect=`cat readparticall.tmp`
  chartest="csk"
  # - - test on a parallel simulation subdirectory cskiiiiii:
  kbeg=`awk -v a="$cskdirect" -v b="$chartest" 'BEGIN{print index(a,b)}'`
  if [ $kbeg -gt 0 ] ; then
    echo "      "
    # - - - extract corsika run number from path name:
    sed -i "/csk/ s/\/csk/\/   /" readparticall.tmp
    jrunnr=`cat readparticall.tmp | awk '{ printf("%06d",$2) }'`
    if [ -e "DAT$jrunnr-000001" ] ; then
       # newer naming of 'parallel' files:
       ls -1 DAT$jrunnr-0????? > readparticall.i$jrunnr
    elif  [ -e DAT$jrunnr-?????????-000000001 ] ; then
       # older naming of 'parallel' files:  
       ls -1 DAT$jrunnr-?????????-????????? > readparticall.i$jrunnr
    fi
    jfiles=`wc readparticall.i$jrunnr | awk '{ printf("%d",$1) }'`
    echo "          ... reading $jfiles particle data files in directory"
    echo "              $cskdirect"
    # - - - run fortran program readparticall and display names:
    let "jfiles = $jfiles + 1"
    jrunew=`echo $jfiles | awk '{ printf("%06d",$1) }'`
    ./readparticall < readparticall.i$jrunnr > readparticall.out$jrunnr
    ls -l readparticall.*$jrunnr
    ls -l DAT$jrunnr-$jrunew*
    echo "      "
  else
    # - - - not within a parallel simulation subdirectory 'cskiiiiii/':
    echo "         ... being not within a 'cskiiiiii/' subdirectory,"
    echo "         ... or missing a name of a particle data file!"
    echo "             ./readparticall.sh DAT987654"
  fi
  #
fi
rm -f readparticall.tmp
