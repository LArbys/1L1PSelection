#!/bin/bash

source setup_tufts_container.sh
source /cluster/tufts/wongjiradlab/twongj01/production/dllee_unified_v1_0_5/configure.sh

export UBDLANA_DIR=$PWD
echo "UBDLANA_DIR: $UBDLANA_DIR"

export MPIDMODEL_DIR=/cluster/tufts/wongjiradlab/twongj01/production/MPID_pytorch
export UBMPIDNET_DIR=/cluster/tufts/wongjiradlab/twongj01/production/MPID_pytorch
export UBMPIDNET_WEIGHT_DIR=/cluster/tufts/wongjiradlab/twongj01/production/MPID_pytorch/weights/

export PYTHONPATH=$UBDLANA_DIR:$PYTHONPATH
export PYTHONPATH=$UBMPIDNET_DIR/:$PYTHONPATH
export PYTHONPATH=$MPIDMODEL_DIR/uboone:$PYTHONPATH
