#!/bin/bash
#
# showanalyrun.sh:
# ================
#     to run showanalyhc3 analysis program on campus grid
#     for corsika run number 2083 (parallel simulation):
# @ job_type = serial
# @ initialdir = /fzk/cgwork/joe/corsika.conex/run/csk002083
# @ environment = COPY_ALL;TMPDIR=/tmp
# @ restart = no
# @ requirements = (Arch == "x86_64") && (Feature == "penryn")
# @ class = medium
# @ wall_clock_limit = 140:00:00
# @ output = showanalyrun.$(cluster).out
# @ error = showanalyrun.$(cluster).err
# @ queue
  cd /fzk/cgwork/joe/corsika.conex/run/csk002083
  ./showanalyscc.sh002083
