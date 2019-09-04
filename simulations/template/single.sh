#!/usr/bin/env bash

# Example single-core job script

here=$(pwd)
infile=$here/input.txt                          # <-- Edit this if you need to
outfile=$here/output.txt                        # <-- Edit this if you want
executable="corsika77000Linux_OPTIONS"          # <-- !!! Definitely Edit This !!!
dir="your/path/to/v77000/run"                   # <-- !!! Definitely Edit This !!!

# to run:
# $ ./single.sh

# Optional, when you're done, if you've compiled corsika2root (see ../../coast/README)
# and didn't "coconut" ROOTOUT, do this from the command line:
# $ ../your/path/coast/v4r5/bin/coast2root ./DATxxxxxxx  <-- whatever your DAT file is
# (note: be sure to have sourced your/path/to/root/thisroot.sh version < 6 before-hand)
#
# Check out some PyROOT tools in ../tools

#-----------------------------------------------------------------------------

errfile=$here/.err
perfile=$here/.per
uname -a > $outfile
echo >> $outfile
date >> $outfile
cd $dir
/usr/bin/time -v -o $perfile ./$executable < $infile >> $outfile 2> $errfile
printf '\n\nErrors (if any):\n' >> $outfile
printf '================\n\n'   >> $outfile
cat $errfile >> $outfile
printf '\n\nPerformance:\n' >> $outfile
printf '============\n\n'   >> $outfile
cat $perfile >> $outfile
rm $errfile $perfile

