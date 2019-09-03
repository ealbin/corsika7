#!/bin/bash
#
# showpllparams.sh:
# =================
# create a tabular of all available parallel corsika steering files:
# ------------------------------------------------------------------------
# check:
#        gfortran -O0 -fbounds-check showpllparams.f -o showpllparams
#        ifort -C -O0 -check bounds showpllparams.f -o showpllparams
# ------------------------------------------------------------------------
# usage: ./showpllparams.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# - - - - write file names via `ls -l` or `ls -1` to text file:
#
# - - - - test login node and showpllparams code:
computype="ikp"
if [ ! -e showpllparams.f ] ; then
  echo "          ----- missing source file 'showpllparams.f'!"
  exit 1
fi
# compile and link considering the login node:
if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
  gfortran -O0 -fbounds-check showpllparams.f -o showpllparams
elif [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
  ifort -C -O0 -check bounds showpllparams.f -o showpllparams
elif [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
  ifort -C -O0 -check bounds showpllparams.f -o showpllparams
fi
#
# ls -1 parallel-* > showpllparams.inpll
  ls -l parallel-* > showpllparams.inpll
# - - - - compile and link fortran source code:
  ifort -C -O0 -check bounds showpllparams.f -o showpllparams
# - - - - run program and display tabular of parameters:
  ./showpllparams < showpllparams.inpll > showpllparams.printout
  rm -rf showpllparams.inpll
  cat showpllparams.printout
