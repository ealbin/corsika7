#!/bin/bash
# 
# sum up all *.long files of a parallel simulation by scripts to one file
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check sumlongifiles.f -o sumlongifiles
# ifort -C -O0 -check bounds sumlongifiles.f -o sumlongifiles
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# 
  ls -1 DAT002604-*.long > sumlongifiles.i002604
# ./sumlongifiles < sumlongifiles.i002604 > sumlongifiles.out002604
# cp -p sumlongifiles.sum002604 DAT002604-999989999.long
