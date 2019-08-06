#!/bin/bash
#
# job_submit -p1 -cp -t35 -m1000 sumlistnkginfo.sh000345
#
# sum up all NKG averages of all `DAT000345-*.lst` files:
# --------------------------------------------------------------------
  ls -1 DAT000345-*.lst > sumlistnkginfo.i000345
  ./sumlistnkginfo < sumlistnkginfo.i000345 > sumlistnkginfo.out000345
