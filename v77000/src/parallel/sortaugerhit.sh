#!/bin/bash
#
# sortaugerhit.sh:
# ================
# create file list and run `sortaugerhit` program;
# run script with a run number or file name as argument;
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
# ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
#               (cpu time about 4300 minutes for 20 cores and 1 TBy space)
# ------------------------------------------------------------------------
# usage: ./sortaugerhit.sh 123456 
#        ./sortaugerhit.sh DAT123456 
#        ./sortaugerhit.sh DAT123456.part 
#        ./sortaugerhit.sh                       # being in cskiiiiii/ #
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# 
runinpart=`echo $1` 
computype="ikp"
pwd > sortaugerhit.tmpname
workdir=`cat sortaugerhit.tmpname`
strtest="csk"
kbeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
# test on source code sortaugerhit.f:
if [ ! -e sortaugerhit.f ] ; then  
  echo "          ----- missing source file 'sortaugerhit.f'!"
  exit 1
fi
# compile and link considering the login node:
if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
  gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
elif [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
  ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
elif [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
  ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $# -eq 0 ] ; then
  # - - - - - - no execution argument given - - - - - - 
  if [ $kbeg -gt 0 ] ; then
    # being in a parallel simulation subdirectory 'cskiiiiii/':
    sed -i "/csk/ s/csk/   /" sortaugerhit.tmpname
    runnrpart=`cat sortaugerhit.tmpname | awk '{printf("%06d",$2)}'`  
    echo "          being in subdirectory csk$runnrpart/" 
    rm -f sortaugerhit.tmpname
    if [ -e DAT$runnrpart-000001 ] ; then
      ls -1 DAT$runnrpart-0????? > sortaugerhit.i$runnrpart
    else
      ls -1 DAT$runnrpart-?????????-????????? > sortaugerhit.i$runnrpart
    fi
    # run program sortaugerhit:
    echo "./sortaugerhit < sortaugerhit.i$runnrpart > sortaugerhit.out$runnrpart"
    # mv fort.20 DAT$runnrpart.augerhit20.txt 
    # ls -l DAT$runnrpart-augerhit* 
    # ls -l sortaugerhit.*$runnrpart
  else # not within a 'cskiiiiii/' subdirectory:
    echo "          ----- being not within a 'cskiiiiii/' subdirectory;"
    echo "          ----- or missing runnr or particle data file name."
  fi
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $# -eq 1 ] ; then
  # - - - - - - name or number as execution argument given - - - - - - 
  if [ -e "$runinpart" ] ; then
    ls -l $runinpart
    datnrpart=`echo $runinpart`
    workdir=`echo $datnrpart`
    strtest="_"
    ibeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
    if [ $ibeg -gt 0 ] ; then
      runnrpart=987654
      echo $datnrpart > sortaugerhit.i$runnrpart
    else
      strtest=".part"
      jbeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
      if [ $jbeg -gt 0 ] ; then 
        echo $datnrpart > sortaugerhit.tmpname 
        sed -i "/DAT/ s/\./ /" sortaugerhit.tmpname
        sed -i "/DAT/ s/DAT/   /" sortaugerhit.tmpname
        runnrpart=`cat sortaugerhit.tmpname | awk '{printf("%06d",$1)}'`
        rm -f sortaugerhit.tmpname
        echo $datnrpart > sortaugerhit.i$runnrpart 
        # files like `prE16M562_41` not yet supported to run.
        # other primaries of this kind were co, gm, he, fe, si.
      else
        strtest="DAT"
        lbeg=`awk -v a="$workdir" -v b="$strtest" 'BEGIN{print index(a,b)}'`
        if [ $lbeg -gt 0 ] ; then
          echo $datnrpart > sortaugerhit.tmpname
          sed -i "/DAT/ s/DAT/   /" sortaugerhit.tmpname
          runnrpart=`cat sortaugerhit.tmpname | awk '{printf("%06d",$1)}'`
          rm -f sortaugerhit.tmpname
          echo $datnrpart > sortaugerhit.i$runnrpart
        fi
      fi
    fi
  else
   runnrpart=`echo $runinpart | awk '{printf("%06d",$1)}'`
   if [ -e "DAT$runnrpart" ] ; then
     ls -l DAT$runnrpart
     ls -1 "DAT$runnrpart" > sortaugerhit.i$runnrpart 
   else
     if [ -e "DAT$runnrpart-999999" ] ; then
       ls -l DAT$runnrpart-999999
       ls -1 "DAT$runnrpart-999999" > sortaugerhit.i$runnrpart
     else
       if [ -e "DAT$runnrpart-000001" ] ; then
         ls -l DAT$runnrpart-000001
         ls -1 DAT$runnrpart-0????? > sortaugerhit.i$runnrpart
       else
         if [ -e "DAT$runnrpart.part" ] ; then 
           ls -l DAT$runnrpart.part
           ls -1 "DAT$runnrpart.part" > sortaugerhit.i$runnrpart
         else
           echo "           case 999999"
           runnrpart=999999
           echo "DAT$runnrpart" > sortaugerhit.i$runnrpart
         fi
       fi
     fi
   fi
  fi
  rm -f sortaugerhit.tmpname
  # run program sortaugerhit:
  echo "./sortaugerhit < sortaugerhit.i$runnrpart > sortaugerhit.out$runnrpart"
  # mv fort.20 DAT$runnrpart.augerhit20.txt
  if [ $runnrpart -lt 987654 ] ; then
  # ls -l DAT$runnrpart[-,.]augerhit* sortaugerhit.*$runnrpart
  else
  # runinpcut=`echo $runinpart | cut -b 1-9 | awk '{printf("%s_augerhit*",$1)}'`
  # ls -l $runinpcut sortaugerhit.*$runnrpart 
  fi
fi
