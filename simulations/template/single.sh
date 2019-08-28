#!/usr/bin/env bash

# Example single-core job script

here=$(pwd)
infile=$here/input.txt                          # <-- Edit this if you need to
outfile=$here/single_output.txt                 # <-- Edit this if you want
executable="corsika76900Linux_QGSJET_gheisha"   # <-- !!! Definitely Edit This !!!

cd /your/path/to/corsika-76900/run              # <-- Edit this

/usr/bin/time -v -a -o $outfile ./$executable < $infile > $outfile 

# run:
# $ ./single.sh

# Optional, when you're done, if you've compiled corsika2root (see ../../coast/README)
# Do this from the command line:
# $ ../../coast/v4r5/bin/coast2root DATxxxxxxx  <-- whatever your DAT file is
# (note: be sure to have sourced your/path/to/root/thisroot.sh version < 6 before-hand)
#
# Check out some PyROOT tools in ../tools
