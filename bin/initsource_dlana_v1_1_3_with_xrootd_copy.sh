#!/bin/bash

unsetup ubdl
ups list -aK+ dllee_unified
ups list -aK+ ubdlana
ups list -aK+ ubMPIDnet
setup dllee_unified v1_0_6 -q e17:prof
source $DLLEE_UNIFIED_BASEDIR/configure.sh
setup ubdlana v1_1_3 -q e17:prof

export IFDH_COPY_XROOTD=1
if [ x$TMPDIR = x ]; then
 export TMPDIR=.
fi

export


