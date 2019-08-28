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
outfile=$here/parallel_output.txt                       # <-- Edit this if you want
executable="corsika76900Linux_QGSJET_gheisha_parallel"  # <-- !!! Definitely Edit This !!!

cd /your/path/to/corsika-76900/run                      # <-- Edit this

/usr/bin/time -v -a -o $outfile srun ./$executable < $infile > $outfile 

# submit this to slurm with:
# $ sbatch parallel.sh

# Optional, when you're done, if you've compiled corsika2root (see ../../coast/README)
# Do this from the command line:
# $ ../../coast/v4r5/bin/coast2root DATxxxxxxxx  <-- whatever your DAT file is
# (note: be sure to have sourced your/path/to/root/thisroot.sh version < 6 before-hand)
#
# Check out some PyROOT tools in ../tools
