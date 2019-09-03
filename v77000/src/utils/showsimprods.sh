#!/bin/bash
#
# showsimprods.sh:
# ================
#     reading corsika particle data files of special user `simprod`
#     with file names like `pre16m316_17` or `sie17m178_09` a.s.o. 
#     a one line per file tabular will be created with the following
#     quantities 
#         primary, lg(E), theta, phi, nsh, runnr, size,
#             obslvme, h1stme, thilev, thiwmax, lg(thirad),
#                 verspgm, models, rundate, Xmagn, Zmagn;
# ----------------------------------------------------------------------
#                                 juergen.oehlschlaeger@kit.edu
# ----------------------------------------------------------------------
#
  ls -l [c,f,g,h,m,p,s]?[e,E]* | ./showsimprods > showsimprods.d11lx24-qgs_cont-2z
#
#     gfortran -fbounds-check showsimprods.f -o showsimprods
#     f77 -fbounds-check showsimprods.f -o showsimprods
#
