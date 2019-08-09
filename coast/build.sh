#!/usr/bin/env bash

# >>> Edit this if needed <<<
COAST_VERSION='v4r5'
# optional feature packages:
#export COAST_USER_LIB="$(pwd)/src/coast-interfaces-v4r1p3/Histogram"

if [ -z $ROOTSYS ]
then
    echo '$ROOTSYS not defined.'
    exit 1
fi

ROOT_VERSION=$(root -q | grep 'Version' | awk '{print $3}' | awk -F '.' '{print $1}') || exit 1
if (( ROOT_VERSION > 5 ))
then
    echo 'Only ROOT versions <= 5 are supported'
    exit 1
fi    

export COAST_DIR="$(pwd)/$COAST_VERSION"
[ ! -d $COAST_DIR ] && mkdir $COAST_DIR
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$COAST_DIR/lib"

cd "./src/coast-$COAST_VERSION"
./configure --with-root=yes && \
make install && \
printf "\n\nFinished successfully!\n"

