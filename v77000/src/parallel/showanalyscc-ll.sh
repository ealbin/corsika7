#!/bin/bash
#
# showanalyscc.sh002139:
# ======================
# - - - - create histograms as ascii tables of all data files 
#         for a parallel corsika simulation on campus grid (opus cluster). 
#
  ls -1 DAT002139-0* | grep t -v | grep n -v > showanalyscc.i002139-0
  ls -1 DAT002139-1* | grep t -v | grep n -v > showanalyscc.i002139-1
  ls -1 DAT002139-2* | grep t -v | grep n -v > showanalyscc.i002139-2
  ls -1 DAT002139-3* | grep t -v | grep n -v > showanalyscc.i002139-3
  ls -1 DAT002139-4* | grep t -v | grep n -v > showanalyscc.i002139-4
  ls -1 DAT002139-5* | grep t -v | grep n -v > showanalyscc.i002139-5
  ls -1 DAT002139-6* | grep t -v | grep n -v > showanalyscc.i002139-6
  ls -1 DAT002139-7* | grep t -v | grep n -v > showanalyscc.i002139-7
  ls -1 DAT002139-8* | grep t -v | grep n -v > showanalyscc.i002139-8
  cat showanalyscc.i002139-* > showanalyscc.i002139
  ./showanalyscc < showanalyscc.i002139 > showanalyscc.out002139
  mv fort.9 showanalyscc.fort002139
