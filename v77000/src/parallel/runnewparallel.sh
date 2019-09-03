#!/bin/bash
#
# runnewparallel.sh:
# ==================
# (a) runnr >= 11000, runnr < 100:
#         copy set of 4 steering files to new parallel coreas simulation.
# (b) runnr < 11000:
#         copy steering file and runscript of a new simulation.
# ------------------------------------------------------------------------
# usage: ./runnewparallel.sh <existnumb> <newnumb>
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
computype="ik3"
if [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
  computype="bwc"
else
  if [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
    computype="fh1"
  fi
fi
#
userhomepath=`echo $USER | awk '{printf("/cr/users/%s",$1)}'`
if [ $computype == "bwc" ] ; then
  userhomepath=`echo $USER | awk '{printf("/work/kit/ikp/%s",$1)}'`
else
  if [ $computype == "fh1" ] ; then
    userhomepath=`echo $USER | awk '{printf("/home/kit/ikp/%s",$1)}'`
  fi
fi
# 
# ------------------------------------------------------------------------
#
if [ $# -eq 2 ]; then
  # - - - - - - - two arguments given:
  jrunarg1=`echo $1`
  jrunarg2=`echo $2`
  #
  jrunshort=`echo $jrunarg1 | awk '{ printf("%d",$1)}'`
  jnewshort=`echo $jrunarg2 | awk '{ printf("%d",$1)}'`
  #
  numberold=`echo $jrunarg1 | awk '{ printf("%06d",$1)}'`
  numbernew=`echo $jrunarg2 | awk '{ printf("%06d",$1)}'`
  #
  parallold=`echo $numberold | awk '{ printf("parallel-%06d",$1)}'`
  parallnew=`echo $numbernew | awk '{ printf("parallel-%06d",$1)}'`
  #
  if [ -e "csk$numbernew/$parallnew" ] ; then
    echo "         ... simulation path already available csk$numbernew/"
    exit
  fi
  if [ $jrunarg1 -ge 11000 ] || [ $jrunarg1 -lt 100 ] ; then
    cp -p "SIM$numberold.inp" "$parallold"
    cp -p "SIM$numberold.sh" jobwcl-$numberold
    if [ $computype == "ik3" ] ; then
       cp -p "SIM$numberold.sh" jobik3-$numberold
    fi
  fi
  if [ ! -e $parallold ]; then
    echo "         ... first runnumb parallel-.. and jobwcl-.. must exist!"
    exit
  fi
  #
  # - - - - - edit parallel steering file:
  cp -p "$parallold" "$parallnew"
  sed -i "/RUNNR/ s/$jrunshort/$jnewshort/" $parallnew 
  sed -i "/$numberold/ s/$numberold/$numbernew/g" $parallnew
  #
  # - - - - - edit jobwcl file:
  cp -p jobwcl-$numberold jobwcl-$numbernew
  sed -i "/$numberold/ s/$numberold/$numbernew/g" jobwcl-$numbernew
  #
  # - - - - - if coreas simul. copy SIM files:
  if [ -e "SIM$numberold.reas" ]; then
    cp -p "SIM$numberold.reas" "SIM$numbernew.reas" 
    cp -p "SIM$numberold.list" "SIM$numbernew.list" 
    cp -p jobwcl-$numbernew "SIM$numbernew.sh"
    cp $parallnew "SIM$numbernew.inp"
  fi
  ls -l *$numbernew*
# 
# ------------------------------------------------------------------------
#
else
  # - - - - - - - less arguments given: 
  if [ $# -eq 1 ]; then
    jrunarg1=`echo $1`
    jrunshort=`echo $jrunarg1 | awk '{ printf("%d",$1)}'`
    numberold=`echo $jrunarg1 | awk '{ printf("%06d",$1)}'`
    parallold=`echo $numberold | awk '{ printf("parallel-%06d",$1)}'`
    
    if [ $jrunarg1 -ge 11000 ] || [ $jrunarg1 -lt 100 ] ; then
      cp -p "SIM$numberold.inp" "$parallold"
      cp -p "SIM$numberold.sh" jobwcl-$numberold
      if [ $computype == "ik3" ] ; then
        cp -p "SIM$numberold.sh" jobik3-$numberold
      fi
    else

      if [ -e "jobwcl-$numberold" ] ; then
        cp -p "jobwcl-$numberold" jobik3-$numberold
      else
         if [ -e "jobik3-$numberold" ] ; then
           cp -p "jobik3-$numberold" jobwcl-$numberold  
         fi
      fi
    fi

    if [ -e $parallold ]; then
      echo "         ... need second number, first does exist!"    
    fi
  else
    echo "         ... need two numbers, first must already exist! Stop!"
  fi
  #
fi
