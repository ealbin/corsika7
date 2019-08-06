#!/bin/bash
#
# job_submit -p1 -cp -t300 -m1000 shdatareduction.sh000351
# 
# create file list and run `shdatareduction` program:
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk000351/
# 
ls -1 DAT000351* | grep t -v | grep l -v | grep a -v > shdatareduction.i000351
# 
# names of sub paths csk??????;
# gfortran -fbounds-check shdatareduction.f -o shdatareduction
# f77 -fbounds-check shdatareduction.f -o shdatareduction
# ifort -C -check bounds shdatareduction.f -o shdatareduction
#
  ./shdatareduction < shdatareduction.i000351 > shdatareduction.out000351
