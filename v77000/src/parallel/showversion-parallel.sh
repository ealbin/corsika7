#!/bin/bash
#
# showversion-parallel.sh
# =======================
#        script to display numbers of corsika versions used in
#        (parallel) simulations written to this path of DAT*.lst
#        and to display names of executables used in parallel
#        corsika simulations on the HP XC3000 at KIT_CN.
# ------------------------------------------------------------------------
# usage: ./showversion-parallel.sh
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# - - - - check parallel subdirectories csk00????/:
grep " mpirun " csk00????/jobhc3-* > version-hc3cors.tmp
ncsk=`wc version-hc3cors.tmp | awk '{ printf("%d",$1) }'`
if [ $ncsk -gt 0 ]; then
  echo "         Number of csk00???? subpaths found: " $ncsk
  chabeg=`head -1 "version-hc3cors.tmp" | cut -b 1-3`  
  info=0
  for i in $( cat "version-hc3cors.tmp" ) ; do
    let "info = info + 1"
    if [ ${#i} -eq 1 ]; then
      let "info = info - 1"
    else
      if [ $info -eq 1 ]; then
        cskversion=`echo $i`
        echo $cskversion > version-corstext.tmp
        sed -i "/csk/ s/jobhc3-/   /" version-corstext.tmp
        cskdirect=`cat "version-corstext.tmp" | awk '{ printf("%s",$1) }'`
        /bin/rm version-corstext.tmp
      else
        if [ $info -eq 9 ]; then
          versioninp=`echo $i`
          echo $versioninp > version-infotext.tmp
          sed -i "/corsika/ s/runner//" version-infotext.tmp 
          sed -i "/corsika/ s/mpi_/_/" version-infotext.tmp
          versiontmp=`cat version-infotext.tmp | awk '{ printf("%s",$1) }'`
          echo "        " $cskdirect "  " $versiontmp
          /bin/rm version-infotext.tmp
          let "info = info - 10"
        fi
      fi
    fi
  done
  /bin/rm version-hc3cors.tmp
fi
#
