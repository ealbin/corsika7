#!/usr/bin/env bash

# >>>>> EDIT THESE <<<<<
name="proton-1e12"
executable="corsika76900Linux_QGSJET_gheisha"
#===========================================

in_file="$(pwd)/input-$name.txt"
out_file="$(pwd)/output-$name.txt"

cd ../../corsika-76900/run

/usr/bin/time -v -a -o $out_file ./$executable < $in_file > $out_file && \
printf "\nFinished!\n" || \
printf "\nFAIL!\n"

