#!/usr/bin/env bash

file=./v77000/src/corsika.F 

if [ ! -d patch ] 
then
   echo "Please re-run this script from the main directory:"
   echo "  $ ./patch/sibyll.sh"
   exit 1
fi 

echo "Adding '|| __SIBYLL__' to $file, line 30446"
sed -i '30446 {s/$/ || __SIBYLL__/}' $file
echo "Finished, try running ./77000/coconut again"
