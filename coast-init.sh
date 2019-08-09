#!/usr/bin/env bash

COAST_VERSION='v4r5'

export COAST_DIR="$(pwd)/coast/$COAST_VERSION"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$COAST_DIR/lib"
export COAST_USER_LIB="$(pwd)/coast/src/coast-interfaces-v4r1p3/Histogram"

