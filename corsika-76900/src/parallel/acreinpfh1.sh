#!/bin/bash
#
# acreinpfh1.sh:
# ==============
# Automatic creation of single or successive steering files named
# `parallel-0iiiii` and shell script files named `jobfh1-0iiiii`
# to run parallel corsika simulations on the KIT-CS fh1 processors
# (i.e. fh1.scc.kit.edu) with MPI parallelization system;
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# compilation:
#     gfortran -O0 -fbounds-check acreinpfh1.f -o acreinpfh1
#     ifort -C -O0 -check bounds acreinpfh1.f -o acreinpfh1
#
# execution:
  ifort -C -O0 -check bounds acreinpfh1.f -o acreinpfh1
  ./acreinpfh1 > acreinpfh1.tabout
  cat acreinpfh1.tabout
  chmod +x jobfh1-*
#
# now jobfh1-iiiiii, parallel-iiiiii (and SIMiiiiii.*) exist. 
# 
# submit of script SIMiiiiii.sh for radio coreas run: 
#      msub SIMiiiiii.sh
# 
# submit of script jobfh1-iiiiii for non-radio simulation:
#      msub jobfh1-iiiiii    
#
