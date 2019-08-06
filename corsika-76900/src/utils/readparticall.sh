#!/bin/bash
#
# readparticall.sh:
# =================
# run particle analysis to count all particles of all data files. 
# ------------------------------------------------------------------------
# check:
#        gfortran -O0 -fbounds-check readparticall.f -o readparticall
#        ifort -O0 -C -check bounds readparticall.f -o readparticall
# ------------------------------------------------------------------------
# usage: ./readparticall.sh [particledatafilename]
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd corsikapath/
#
rm -f readparticall.inpfile
if [ $# -eq 1 ] ; then
# - - - sequential corsika simulation (one file as argument):
  echo $1 > readparticall.inpfile
  echo $1 > jobfilename.tmp
  sed -i "/DAT/ s/\./ /" jobfilename.tmp
  sed -i "/DAT/ s/DAT/   /" jobfilename.tmp
  jrunnr=`cat jobfilename.tmp | awk '{printf("%06d\n",$1)}'`
  # ls -l readparticall.inpfile # already created.
  rm -f jobfilename.tmp
else
# - - - parallel corsika simulation (is current path):
  pwd > jobfilepwd.tmp
  # check position of string `csk`:
  workdir=`cat "jobfilepwd.tmp"`
  strtest="csk"
  ibeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
  if [ $ibeg -gt 0 ] ; then
    let "iend = $ibeg + 2"
    # but use awk to separate run number:
    sed -i "/csk/ s/csk/   /" jobfilepwd.tmp
    jrunnr=`cat jobfilepwd.tmp | awk '{printf("%s",$2)}'`
    if [ -e DAT$jrunnr-000001 ] ; then
      ls -1 DAT$jrunnr-\[0-8\]????? > readparticall.i$jrunnr
    else
      ls -1 DAT$jrunnr-?????????-????????? > readparticall.i$jrunnr
    fi
    cp readparticall.i$jrunnr readparticall.inpfile
    rm -f jobfilepwd.tmp
  fi
fi
# - - - run fortran program readparticall:
# gfortran -O0 -fbounds-check readparticall.f -o readparticall
# ifort -O0 -C -check bounds readparticall.f -o readparticall
if [ -e "readparticall.inpfile" ]; then
  echo "             ... writing protocol file readparticall.out$jrunnr"
  ./readparticall < readparticall.inpfile > readparticall.out$jrunnr
  rm -f readparticall.inpfile
else 
  echo "             ... No corsika particle data file given as argument."
fi
