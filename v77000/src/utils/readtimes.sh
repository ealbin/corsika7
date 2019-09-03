#!/bin/bash
#
# readtimes.sh:
# =============
# print start times and end times of a set of corsika simulations
# using `.lst` files to get the sum of all run times (hours, days).
# Command `head -108 ` must be checked for your simulation group,
# the given line must contain the protocol line like
#                 PRESENT TIME : 08.10.2012  14:57:41
# exactly two lines after the print out of `== START OF RUN ==`. 
# ------------------------------------------------------------------------
#              gfortran -fbounds-check readtimes.f -o readtimes
#              f77 -fbounds-check readtimes.f -o readtimes
#              ifort -C -check bounds readtimes.f -o readtimes
# ------------------------------------------------------------------------
# usage: ./readtimes.sh DAT987 [ cd9-e19m100as ]
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
listdats=`echo $1` 
if [ -n "$listdats" ]; then
  # - - - - - - length of argument is > 0:
  listpath=`echo $2`
  if [ ! -n "$listpath" ]; then
     listpath=`echo "lst" | awk '{ printf("%s",$1) }'`
  fi
  # - - - one or two arguments defined: 
  timeinpt=`echo $listpath | awk '{ printf("readtimes.input-%s",$1) }'`
  echo "        ./readtimes <" $timeinpt
  if [ -e $timeinpt ] ; then
     rm $timeinpt ;
  fi
  touch $timeinpt
  listfiles=`echo $listdats | awk '{ printf("%s*.lst",$1) }'`
  for i in $( ls -1 $listfiles ) ; do 
    ls -1 $i >> $timeinpt
    head -104 $i | tail -1 >> $timeinpt # check 108 for simulation group!
    tail -3 $i >> $timeinpt
  done
  # - - - readtimes.f must be compiled and linked completely.
  ./readtimes < $timeinpt
#
# - - - - - - length of 1st argument = 0 (missing arguments):
else
  echo "           ... missing simulation group like \"DAT987\" as 1st arg."
  echo "           ... (optionally path argument like \"cd9-e19m100as\")"
fi
#
