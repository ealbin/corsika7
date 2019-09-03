#!/bin/bash
#
# showversion-corsika.sh
# ======================
#        script to display numbers of corsika versions used in
#        (parallel) simulations written to this path of DAT*.lst
#        and to display names of executables used in parallel
#        corsika simulations on the HP XC3000 at KIT_CN.
# ------------------------------------------------------------------------
# usage: ./showversion-corsika.sh
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# - - - - check regular DAT*.lst files:
grep " OF VERSION" DAT??????.lst > version-corsika.tmp
mlst=`wc version-corsika.tmp | awk '{ printf("%d",$1) }'`
let "nlst = mlst / 2"
#
# - - - - work with regular DAT*.lst files:
if [ $mlst -gt 0 ]; then
  echo "         Number of DAT??????.lst files found: " $nlst
  chabeg=`head -1 "version-corsika.tmp" | cut -b 1-3`  
  info=0
  for i in $( cat "version-corsika.tmp" ) ; do
    let "info = info + 1"
    if [ $info -eq 1 ]; then
       lstfile=`echo $i`
    fi
    if [ $info -eq 6 ]; then
       numb=`echo $i`
    fi
    if [ $info -eq 12 ]; then
       month=`echo $i`
    fi
    if [ $info -eq 13 ]; then
       day=`echo $i`
    fi
    if [ $info -eq 14 ]; then
       year=`echo $i`
       echo "     " $lstfile " " $numb "   " $day\. "" $month "" $year
       info=0
    fi
  done
  /bin/rm version-corsika.tmp
fi
#
# - - - - use grep command and awk print out:
grep " < aug" saug?????? | awk '{ printf("         %s   %s \n",$1,$2) }'
#
