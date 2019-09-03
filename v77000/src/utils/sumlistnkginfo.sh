#!/bin/bash
#
# sum up all NKG averages of all `DAT000345-*.lst` files:
# --------------------------------------------------------------------
# gfortran -O0 -fbounds-check sumlistnkginfo.f -o sumlistnkginfo
# ifort -C -O0 -check bounds sumlistnkginfo.f -o sumlistnkginfo
# --------------------------------------------------------------------
  ls -1 DAT000345-*.lst > sumlistnkginfo.i000345
# ./sumlistnkginfo < sumlistnkginfo.i000345 > sumlistnkginfo.out000345
