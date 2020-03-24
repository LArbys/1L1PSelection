#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup numpy v1_14_3 -q e17:p2714b:openblas:prof
setup xgboost  v0_82   -q e17:prof
setup ubutil   v08_00_00_27 -q e17:prof
#setup dllee_unified develop -q e17:prof # comment out this line if using your own copy
setup ubMPIDnet v1_0_0 -q NULL
unsetup dllee_unified

export UBDLANA_DIR=$PWD
export PYTHONPATH=$UBDLANA_DIR/ubdlana:$UBDLANA_DIR:${UBMPIDNET_DIR}/uboone:$PYTHONPATH
export FHICL_FILE_PATH=${UBDLANA_DIR}/cfg:$FHICL_FILE_PATH
