#!/usr/bin/env bash

# >>>>> EDIT THESE <<<<<
name="proton-1e12"
executable="corsika76900Linux_QGSJET_gheisha"
corsika2root="../../coast/v4r5/bin/corsika2root"
#===============================================

here=$(pwd)
in_file="$here/input-$name.txt"
out_file="$here/output-$name.txt"

cd ../../corsika-76900/run

/usr/bin/time -v -a -o $out_file ./$executable < $in_file > $out_file && \
cd $here && \
for f in $(find . -name 'DAT*' ! -name '*.*')
do
    $corsika2root $f
done && \
printf "\nFinished!\n" || \
printf "\nFAIL!\n"

