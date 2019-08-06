#!/bin/bash
#
# acreinpbwc.sh:
# ==============
# Automatic creation of single or successive steering files named
# `parallel-iiiiii` and shell script files named `jobwcl-iiiiii`
# to run parallel corsika simulations on the KIT-CS bwunicluster 
# (i.e. bwunicluster.scc.kit.edu) with MPI parallelization system;
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# compilation:
#     gfortran -O0 -fbounds-check acreinpbwc.f -o acreinpbwc
#     ifort -C -O0 -check bounds acreinpbwc.f -o acreinpbwc
#
# execution:
  ifort -C -O0 -check bounds acreinpbwc.f -o acreinpbwc
  ./acreinpbwc > acreinpbwc.tabout
  cat acreinpbwc.tabout
  chmod +x jobwcl-*
#
# now jobwcl-iiiiii, parallel-iiiiiii (and all SIMiiiiii.*) exist. 
# 
# submit of script SIMiiiiii.sh for radio coreas run: 
#      msub SIMiiiiii.sh
# 
# submit of script jobwcl-iiiiii for non-radio simulation:
#      msub jobwcl-iiiiii    
#
