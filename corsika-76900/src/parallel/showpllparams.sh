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
# ls -1 parallel-* > showpllparams.inpll
  ls -l parallel-* > showpllparams.inpll
# - - - - compile and link fortran source code:
  ifort -C -O0 -check bounds showpllparams.f -o showpllparams
# - - - - run program and display tabular of parameters:
  ./showpllparams < showpllparams.inpll > showpllparams.printout
  cat showpllparams.printout
