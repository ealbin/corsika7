#!/usr/bin/env bash

COAST_VERSION='v4r5'

export COAST_DIR="$(pwd)/$COAST_VERSION"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$COAST_DIR/lib"

# optional COAST tools:
#export COAST_USER_LIB="$(pwd)/coast/src/coast-interfaces-v4r1p3/Histogram"

