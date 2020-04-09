#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

WORKDIR=$1
INFILE_NUMBER=$2
INFILE_LIST=$3
IN_ROOT_DIR=$4



echo "Here"
echo $WORKDIR
echo $INFILE_NUMBER
echo $INFILE_LIST
echo $IN_ROOT_DIR
# go to the working directory
cd $WORKDIR
source setup_tufts_container.sh
source configure.sh
CALIB=/cluster/tufts/wongjiradlab/jmills09/CalibrationMaps_MCC9.root
OUT_ROOT_DIR=/cluster/tufts/wongjiradlab/jmills09/files_merged
lines=`wc -l $INFILE_LIST`
lines=$(< $INFILE_LIST wc -l)
# lines=1
line=0
while [ $line -lt $lines ]
do
  echo
  let line+=1
  THIS_FILE_ID=`sed -n ${line}p ${INFILE_LIST}`
  RELATIVE_PATH=${THIS_FILE_ID%/*}
  FILE_NAME=${THIS_FILE_ID##*/}
  echo Relativepath
  echo ${RELATIVE_PATH}
  echo Name
  echo ${FILE_NAME}
  python MakeFinalVertexFiles-prime-dlmerged_NUMUONLY_mcprecuts.py ${IN_ROOT_DIR}/${THIS_FILE_ID} ${CALIB} ${OUT_ROOT_DIR}/${RELATIVE_PATH}/ > test_outs/log_nue_numu.log
  echo ${RELATIVE_PATH} >> test_outs/successes_nue_numu.txt
done
