#!/bin/bash

if [[ -z "$1" ]]; then
    echo "You need to supply a script name when starting the container."
    exit -1
fi

if [[ -z "$2" ]]; then
    echo "You need to supply a script an init file when starting the container."
    exit -1
fi

export DUNE_CONTROL_PATH=/dune/amdis:/dune
/dune/dune-common/bin/dunecontrol --current all

build-cmake/src/$1 $2
