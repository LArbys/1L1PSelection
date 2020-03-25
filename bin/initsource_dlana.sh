#!/bin/bash

unsetup ubdl
ups list -aK+ dllee_unified
ups list -aK+ ubdlana
ups list -aK+ ubMPIDnet
setup dllee_unified v1_0_5 -q e17:prof
source $DLLEE_UNIFIED_BASEDIR/configure.sh
setup ubdlana develop -q e17:prof

export


