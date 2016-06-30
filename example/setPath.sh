#!/bin/bash

ALG=alldisks
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ppath=`realpath ${DIR}/../build/debug/pyWrapper`

export LD_LIBRARY_PATH=`realpath ${DIR}/../../../../lib`:${LD_LIBRARY_PATH}
export PYTHONPATH=${ppath}:${PYTHONPATH}
