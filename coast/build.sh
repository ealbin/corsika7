#!/usr/bin/env bash

COAST_VERSION='v4r5'

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
export COAST_USER_LIB="$(pwd)/src/coast-interfaces-v4r1p3/Histogram"

cd "./src/coast-$COAST_VERSION"
./configure --with-root=yes && \
make install && \
printf "\n\nFinished successfully!\n"

