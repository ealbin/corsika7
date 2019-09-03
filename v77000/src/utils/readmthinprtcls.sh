#!/bin/bash
#
# readmthinprtcls.sh:
# ==================
# run multi-thinning analysis including all particle counting. 
# ------------------------------------------------------------------------
# check:
#        gfortran -O0 -fbounds-check readmthinprtcls.f -o readmthinprtcls
#        ifort -O0 -C -check bounds readmthinprtcls.f -o readmthinprtcls
# ------------------------------------------------------------------------
# usage: ./readmthinprtcls.sh [particledatafilename]
#        ls -1 DAT987654 > readmthinprtcls.i987654
#        ls -1 DAT987654-?????? > readmthinprtcls.i987654
#        ./readmthinprtcls < readmthinprtcls.i987654
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd corsikapath/
#
if [ $# -eq 1 ] ; then
pwd > jobfilepwd.tmp
# - - - sequential corsika simulation (one file as argument):
  echo $1 > readmthinprtcls.inpll2asci
  echo $1 > jobfilename.tmp
  sed -i "/DAT/ s/\./ /" jobfilename.tmp
  sed -i "/DAT/ s/DAT/   /" jobfilename.tmp
  jrunnr=`cat jobfilename.tmp | awk '{printf("%06d\n",$1)}'`
else
# - - - parallel corsika simulation (current path):
  # check position of string `csk`:
  workdir=`cat "jobfilepwd.tmp"`
  strtest="csk"
  ibeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
  let "iend = $ibeg + 2"
  # but use awk to separate run number:
  sed -i "/csk/ s/csk/   /" jobfilepwd.tmp
  jrunnr=`cat jobfilepwd.tmp | awk '{printf("%s",$2)}'`
  if [ -e DAT$jrunnr-000001 ] ; then
    ls -1 DAT$jrunnr-?????? > readmthinprtcls.i$jrunnr
  else
    ls -1 DAT$jrunnr-?????????-????????? > readmthinprtcls.i$jrunnr
  fi
  cp readmthinprtcls.i$jrunnr readmthinprtcls.inpll2asci
fi
rm -f jobfilepwd.tmp
# check fortran executable:
if [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
  ifort -C -O0 -check bounds readmthinprtcls.f -o readmthinprtcls
fi
if [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
  ifort -C -O0 -check bounds readmthinprtcls.f -o readmthinprtcls
fi
if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
  gfortran -O0 -fbounds-check readmthinprtcls.f -o readmthinprtcls
fi
# - - - run fortran program readmthinprtcls:
if [ -e "readmthinprtcls.inpll2asci" ]; then
  echo "             ... writing protocol file readmthinprtcls.out$jrunnr"
  ./readmthinprtcls < readmthinprtcls.inpll2asci > readmthinprtcls.out$jrunnr
  rm -rf readmthinprtcls.inpll2asci
else 
  echo "             ... No corsika particle data file given as argument."
fi
