#!/bin/bash
# 
# create the tabular of available parallel corsika simulations:
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check showparallel.f -o showparallel
# ifort -C -O0 -check bounds showparallel.f -o showparallel
# ------------------------------------------------------------------------
# Primary   lg(E)  theta    phi  runtsk  sizeGBy  procs
#         T(days)  ecutmax  t(min)  files  RATIO  obslev
#                  Xmagn  Zmagn  _corsika_executable_
#                  Antns  gamma  _corsika_executable_
#         < ecutha  ecutmu  ecutel  ecutga  thilev  wmax  lg(thirad) >
# names of subdirectories csk[0-7]?????;
# job protocol files Job??????_%jobid.err, Job??????_%jobid.out;
# ------------------------------------------------------------------------
# usage: ./showparallel.sh [reas]
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# 
plltype="cors"
plltext=""
# - - - - - check argument of call `showparallel.sh` and run numbers:
if [ $# -eq 1 ] ; then
  plltype=`echo $1`
  plltext="-coreas"
fi
# - - - - - distinguish argument of showparallel.sh:
if [ $plltype == "cors" ] ; then 
  # - - - - - - no script argument given or arg is `cors`:
  ls -1 csk[0-7]?????/Job*.out > showparallel.jobinfos
  ./showparallel < showparallel.jobinfos > showparallel.bwc-work-jobinfos
else
  # - - - - - - argument `reas` (including coreas simulations):
  if [ $plltype == "reas" ] ; then
    ls -1 csk[0-7]?????/Job*.out > showparallel.jobinfos
    ./showparallel < showparallel.jobinfos > showparallel.bwc-work-jobinfos$plltext
  fi
fi
#
