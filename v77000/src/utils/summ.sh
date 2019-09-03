#!/bin/bash
#
# summ.sh:
# ========
#          display sum of bytes of files of names beginning with ...
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check summe.f -o summe
# ifort -C -O0 -check bounds summe.f -o summe
# ------------------------------------------------------------------------
# usage:    ./summ <datfiles>
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
datfiles=`echo "\"$1\"" | awk '{ printf("%s",$1) }'`
echo "$datfiles"
if [ -n "$datfiles" ]; then
  echo $datfiles 
  # regular beginning of files given:
  echo $datfiles > testdatfiles.tmp 
  sed -i "/\*/ s/\*//" testdatfiles.tmp
  sed -i "/\?/ s/\?//g" testdatfiles.tmp 
  cat testdatfiles.tmp
  begfiles=`cat testdatfiles.tmp | awk '{ printf("%s",$1) }'`
  echo $begfiles
  du -k $begfiles* | ./summe
  rm testdatfiles.tmp
else
  # missing beginning of files:
  echo "sum of bytes: missing beginning of files;"
fi
#
