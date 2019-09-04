#!/bin/env bash

# Example parallel job script for UCI greenplanet cluster

#SBATCH --job-name=my_job               # <-- Edit this
#SBATCH -p a_partition_i_can_use        # <-- Edit this (e.g. nes2.8)
#SBATCH -c 1                            # cpus-per-task, 1 is correct
#SBATCH -n 123                          # <-- Edit this (how many cores do you want)
#SBATCH -t 2-8:00                       # <-- Edit this (days-hours:minutes)
#SBATCH -o /path/slurm_output.txt       # <-- Edit this

here=$(pwd)
infile=$here/input.txt                  # <-- Edit this if you need to (!! don't forget to uncomment PARALLEL in input.txt !!)
outfile=$here/output.txt                # <-- Edit this if you want
executable="corsika77000Linux_OPTIONS_parallel"  # <-- !!! Definitely Edit This !!!
dir="your/path/to/v77000/run"                    # <-- !!! Definitely Edit This !!!

# submit this to slurm with:
# $ sbatch parallel.sh

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
/usr/bin/time -v -o $perfile srun ./$executable < $infile >> $outfile 2> $errfile
printf '\n\nErrors (if any):\n' >> $outfile
printf '================\n\n'   >> $outfile
cat $errfile >> $outfile
printf '\n\nPerformance:\n' >> $outfile
printf '============\n\n'   >> $outfile
cat $perfile >> $outfile
rm $errfile $perfile

