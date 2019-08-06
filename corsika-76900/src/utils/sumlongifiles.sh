#!/bin/bash
#
# job_submit -p1 -cp -t300 -m1000 sumlongifiles.sh002543
# 
# create file list DAT*.long and run `sumlongifiles` program
# to get the corresponding total `DAT*.long` file:  
# ----------------------------------------------------------------------
#                                     juergen.oehlschlaeger@kit.edu
#-----------------------------------------------------------------------
# gfortran -O0 -fbounds-check sumlongifiles.f -o sumlongifiles
# ifort -C -O0 sumlongifiles.f -o sumlongifiles
#
  ls -1 DAT002543-*.long > sumlongifiles.i002543
#
  ./sumlongifiles < sumlongifiles.i002543 > sumlongifiles.out002543
#
  /bin/cp -p sumlongifiles.sum002543 DAT002543-999989999.long
