#!/bin/bash
#
# job_submit -p1 -cp -t3600 -m1000 shdatatanks.sh002084
# 
# create file list and run `shdatatanks` program:
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk002084/
# 
ls -1 DAT002084-* | grep t -v | grep l -v | grep a -v > shdatatanks.i002084
# 
# names of sub paths csk00????;
# gfortran -fbounds-check shdatatanks.f -o shdatatanks
# f77 -fbounds-check shdatatanks.f -o shdatatanks
# ifort -C -check bounds shdatatanks.f -o shdatatanks
#
  ./shdatatanks < shdatatanks.i002084 > shdatatanks.out002084
