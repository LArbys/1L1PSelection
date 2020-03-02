#!/bin/bash

if [ -z ${DLLEE_UNIFIED_BASEDIR+x} ]; then
    export DLLEE_UNIFIED_BASEDIR=/cluster/tufts/wongjiradlab/dcianc01/dllee_unified
fi

# setup environment variables
#source $DLLEE_UNIFIED_BASEDIR/setup.sh

# setup larlite
source $DLLEE_UNIFIED_BASEDIR/larlite/config/setup.sh

# setup laropencv
source $DLLEE_UNIFIED_BASEDIR/LArOpenCV/setup_laropencv.sh

# setup Geo2D
source $DLLEE_UNIFIED_BASEDIR/Geo2D/config/setup.sh

# setup Cilantro
source $DLLEE_UNIFIED_BASEDIR/cilantro/dllee_setup.sh

# setup LArCV
source $DLLEE_UNIFIED_BASEDIR/LArCV/configure.sh

# setup larlitecv
source $DLLEE_UNIFIED_BASEDIR/larlitecv/configure.sh

# ENV VARIABLES TO SETUP BINARIES/RECO SCRIPTS

# Tagger
export TAGGER_BIN_DIR=$DLLEE_UNIFIED_BASEDIR/larlitecv/app/TaggerCROI/bin
[[ ":${PATH}:" != *":${TAGGER_BIN_DIR}:"* ]] && export PATH="${TAGGER_BIN_DIR}:${PATH}"

# Vertexer
export VERTEXER_SCRIPT_DIR=${LARCV_BASEDIR}/app/LArOpenCVHandle/cfg/mac/
[[ ":${PATH}:" != *":${VERTEXER_SCRIPT_DIR}:"* ]] && export PATH="${VERTEXER_SCRIPT_DIR}:${PATH}"
