#!/bin/bash
#
# showsimulist.sh:
# ================
#     create tabular of air shower simulations of corsika
#     by reading first record of all particle data files DATiiiiii
#     and possibly additional reading of the corresponding protocol
#     (or log-) file DATiiiiii.lst; a tabular of one line per file
#     will be created with the following quantities 
# prim, lg(E), theta, phi, nsh, runnr, size, obslevme, h1stme, thilev,
#       thiwmax, lg(thirad), verspgm, models, rundate, Xmagn, Zmagn.
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
  ls -l DAT?????? | ./showsimulist > showsimulist.cdj-corsika-trunk-run
#
# gfortran -O0 -fbounds-check showsimulist.f -o showsimulist
# ifort -C -O0 -check bounds showsimulist.f -o showsimulist
#
