#!/usr/bin/env bash

file=../v77000/nexus/nexus-uti-3972.f 

if [ ! -d patch ] 
then
   echo "Please re-run this script from the main directory:"
   echo "  $ ./patch/sibyll.sh"
   exit 1
fi 

[ ! -f $file ] && \
echo "$file hasn't been extracted yet, try compiling NEXUS with ./v77000/coconut first, and run this if it fails" && \
exit 1

echo "Commenting out the \"sbet\" function in $file, lines 654 to 658"
sed -i '654 {s/^/c/}' $file
sed -i '656,658 {s/^/c/}' $file
sed -i '652 {s/^/\nc !!! Commented out by fix_nexus.sh !!!/}' $file
echo "Finished, try running ./v77000/coconut again"

