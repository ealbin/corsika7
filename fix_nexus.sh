#!/usr/bin/env bash

file=./corsika-76900/nexus/nexus-uti-3972.f 

[ ! -f $file ] && \
echo "$file hasn't been extracted yet, try compiling NEXUS with ./corsika-76900/coconut first, and run this if it fails" && \
exit 1

echo "Commenting out the \"sbet\" function in $file, lines 654 to 658"
sed -i '654 {s/^/c/}' $file
sed -i '656,658 {s/^/c/}' $file
sed -i '652 {s/^/\nc !!! Commented out by fix_nexus.sh !!!/}' $file
echo "Finished, try running ./corsika-76900/coconut again"
