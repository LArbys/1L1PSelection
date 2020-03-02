#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

WORKDIR=$1
INFILE_NUMBER=$2
INFILE_LIST=$3



echo "Here"
echo $WORKDIR
echo $INFILE_NUMBER
echo $INFILE_LIST
# go to the working directory
cd $WORKDIR
source setup_tufts_container.sh
source configure.sh
CALIB=/cluster/tufts/wongjiradlab/jmills09/CalibrationMaps_MCC9.root
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
  python MakeFinalVertexFiles-prime-dlmerged_NUMUONLY_mcprecuts.py ${THIS_FILE_ID} ${CALIB} ${RELATIVE_PATH}/ > test_outs/log_ext_dirt.log
  echo ${RELATIVE_PATH} >> test_outs/successes_ext_dirt.txt
done
