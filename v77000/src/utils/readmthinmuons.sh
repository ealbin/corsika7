#!/bin/bash
#
# readmthinmuons.sh:
# ==================
# run multi-thinning analysis including all particle counting. 
# ------------------------------------------------------------------------
# check:
#        gfortran -O0 -fbounds-check readmthinmuons.f -o readmthinmuons
#        ifort -O0 -C -check bounds readmthinmuons.f -o readmthinmuons
# ------------------------------------------------------------------------
# usage: ./readmthinmuons.sh [DATnnnnnn]
#        ./readmthinmuons < readmthinmuons.innnnnn
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd corsikapath/
#
rm -f readmthinmuons.inpmthi
if [ $# -eq 1 ] ; then
# - - - sequential corsika simulation (one file as argument):
  echo $1 > readmthinmuons.inpmthi
  echo $1 > jobfilename.tmp
  sed -i "/DAT/ s/\./ /" jobfilename.tmp
  sed -i "/DAT/ s/DAT/   /" jobfilename.tmp
  jrunnr=`cat jobfilename.tmp | awk '{printf("%06d\n",$1)}'`
  # ls -l readmthinmuons.inpmthi # already created.
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
      ls -1 DAT$jrunnr-\[0-8\]????? > readmthinmuons.i$jrunnr
    else
      ls -1 DAT$jrunnr-?????????-????????? > readmthinmuons.i$jrunnr
    fi
    cp readmthinmuons.i$jrunnr readmthinmuons.inpmthi
    rm -f jobfilepwd.tmp
  fi
fi
# - - - run fortran program readmthinmuons:
# gfortran -O0 -fbounds-check readmthinmuons.f -o readmthinmuons
# ifort -O0 -C -check bounds readmthinmuons.f -o readmthinmuons
if [ -e "readmthinmuons.inpmthi" ]; then
  echo "             ... writing protocol file readmthinmuons.out$jrunnr"
  ./readmthinmuons < readmthinmuons.inpmthi > readmthinmuons.out$jrunnr
  rm -f readmthinmuons.inpmthi
else 
  echo "             ... No corsika particle data file given as argument."
fi
