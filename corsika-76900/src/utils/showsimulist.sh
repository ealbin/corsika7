#!/bin/bash
#
# showsimulist.sh:
# ================
#           create tabular of air shower simulations of corsika
#           structure by reading first record of all particle
#           data files DATiiiiii and additional reading of the
#           corresponding protocol file DATiiiiii.lst .
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
# prim, lg(E), theta, phi, nsh, task, size, obslev, h1stme, thilev,
#       thiwmax, lg(thirad), verspgm, models, rundate, Xmagn, Zmagn.
#
  ls -l DAT?????? | ./showsimulist > showsimulist.cdj-corsika-trunk-run
#
# usage: ls -l [c,f,p]?e??m???_* ./showsimulist # simprod simulations. 
#
#       gfortran -fbounds-check showsimulist.f -o showsimulist
#       ifort -O0 -C -check bounds showsimulist.f -o showsimulist
#
