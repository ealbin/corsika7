#!/bin/bash
# 
# sumrawcoreas.sh:
# ================
# delete some subsidiary coreas protocol files and prepare
# to sum up antenna data files of all subdirectories `*_coreas`
# to one new file for each antenna (originally performed 
# by script parallelization).
# ------------------------------------------------------------------------
# CompLink:
#     gfortran -O0 -fbounds-check sumrawcoreas.f -o sumrawcoreas
# Running:
#     ./sumrawcoreas > sumrawcoreas.outrunnr
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
# - - - - get run number from current directory name:
pwd > simcoreaspath.tmp
sed -i "/csk/ s/csk/   /" simcoreaspath.tmp
jrunnr=`cat simcoreaspath.tmp | awk '{printf("%06d\n",$2)}'`
nrunnr=`echo $jrunnr | awk '{printf("%d\n",$1)}'`
rm -f simcoreaspath.tmp
#
# - - - - writing all names of `SIM*_coreas/raw*.dat` files:
ls -1 SIM$jrunnr-*_coreas/raw*.dat > SIM$jrunnr.raw.lst
rawfirst=`head -1 "SIM$jrunnr.raw.lst" | awk '{printf("%s",$1)}'`
nrawdat=`wc $rawfirst | awk '{printf("%d",$1)}'`
#
# - - - - keep only last `.reas` file in this path:
reasfiles=`echo $jrunnr | awk '{printf("SIM%06d-*.reas",$1)}'`
ls -1 $reasfiles 2> simreasfiles.err > simreasfiles.txt
rm -f simreasfiles.err
nreas=`wc "simreasfiles.txt" | awk '{printf("%7d\n",$1)}'`
if [ $nreas -gt 0 ] ; then
  let "ndel = $nreas - 1"
  if [ $ndel -gt 0 ] ; then
    cat simreasfiles.txt | head -$ndel > simreasfiles-del.tmp
    for reasname in $( cat "simreasfiles-del.tmp" ) ; do
      # echo $reasname ;  
      rm -f $reasname ;
    done
    rm -f simreasfiles-del.tmp
    echo "      ... deleting almost all SIM$jrunnr-*.reas ..."
  fi
  rm -f simreasfiles.txt
fi
#
# - - - - keep only last `.bins` file in this path:
binsfiles=`echo $jrunnr | awk '{printf("SIM%06d-*.bins",$1)}'`
ls -1 $binsfiles 2> simbinsfiles.err > simbinsfiles.txt
rm -f simbinsfiles.err
nbins=`wc "simbinsfiles.txt" | awk '{printf("%7d\n",$1)}'`
if [ $nbins -gt 0 ] ; then
  let "ndel = $nbins - 1"
  if [ $ndel -gt 0 ] ; then
    cat simbinsfiles.txt | head -$ndel > simbinsfiles-del.tmp
    for binsname in $( cat "simbinsfiles-del.tmp" ) ; do
      # echo $binsname ;  
      rm -f $binsname ;
    done
    rm -f simbinsfiles-del.tmp
    echo "      ... deleting almost all SIM$jrunnr-*.bins ..."
  fi
  rm -f simbinsfiles.txt
fi
ls -1 $binsfiles > simbinslast.tmp
simbinslast=`cat "simbinslast.tmp" | awk '{printf("%s",$1)}'`
#
# - - - - create new directory for summed antenna data:
newnumbers=`echo ${simbinslast:10:9}`
simnewpath=`echo "${simbinslast:0:19}-${newnumbers}_coreas"`
mkdir -p "$simnewpath"
#
# - - - - writing txt file to be read by fortran code:
echo "      ... creating parameter file 'SIMcoreas.rawinfo':" 
echo "          $jrunnr   corsika runnr" > SIM$jrunnr.rawinfo
nrawtotal=`wc "SIM$jrunnr.raw.lst" | awk '{printf("%d",$1)}'`
echo $nrawtotal | awk '{printf("%7d   total number of raw files\n",$1)}' >> SIM$jrunnr.rawinfo
nrawsubdi=`wc "$simbinslast" | awk '{printf("%d",$1)}'`
echo $nrawsubdi | awk '{printf("%7d   raw files per subdirectory (nAnts)\n",$1)}' >> SIM$jrunnr.rawinfo
echo $nrawdat | awk '{printf("%7d   entries per raw file\n",$1)}' >> SIM$jrunnr.rawinfo
echo $PWD >> SIM$jrunnr.rawinfo
cat simbinslast.tmp >> SIM$jrunnr.rawinfo
cp -p SIM$jrunnr.rawinfo SIMcoreas.rawinfo
rm -f simbinslast.tmp
cat SIM$jrunnr.rawinfo
# - - - - compile and link source code and run summation:
gfortran -O0 -fbounds-check sumrawcoreas.f -o sumrawcoreas
echo " ===> run the following command line (reading 'SIMcoreas.rawinfo.')"
echo "          ./sumrawcoreas > sumrawcoreas.out$jrunnr"
# rm -f SIMcoreas.rawinfo
exit
