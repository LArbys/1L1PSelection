#!/bin/bash

# sample name
SAMPLE=mcc9jan_bnb5e19
showerreco=/cluster/tufts/wongjiradlab/twongj01/1L1PSelection/SSNetShowerRecoOut_bnb5e19.txt

# where the selection scripts live
SELECT_DIR=/cluster/tufts/wongjiradlab/twongj01/1L1PSelection/

# files to run
FILELIST=$SELECT_DIR/filelists/dlreco_$SAMPLE.txt

# files to run per job
NFILES=100
CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/singularity_ubdl_deps_py2_10022019.simg

# DLLEE UNIFIED REPO location
DLLEE_DIR=/cluster/tufts/wongjiradlab/twongj01/production/dllee_unified_opencv3.1_copy2

# final output dir
FINAL_OUT_DIR=$SELECT_DIR/out_$SAMPLE
mkdir -p $FINAL_OUT_DIR

# for debug (slurm will set this up for us automatically if we use the array argument)
#SLURM_ARRAY_TASK_ID=10

# define start and end lines in filelist to run
let start=$NFILES*$SLURM_ARRAY_TASK_ID+1
let end=$NFILES*${SLURM_ARRAY_TASK_ID}+$NFILES

# working directory
workdir=/tmp

# need to have unique job directory so jobs dont collide
jobdir=$workdir/selection_1l1p_$SAMPLE_$SLURM_ARRAY_TASK_ID
mkdir -p $jobdir
cd $jobdir

# get list of files
sed -n ${start},${end}p $FILELIST > runlist.txt

# define parameters for running script
outdir=$jobdir/out
calibmap=/cluster/tufts/wongjiradlab/twongj01/1L1PSelection/CalibrationMaps_MCC9.root

# make output directory and clear it out
mkdir -p $outdir
rm -f $outdir/*

# SETUP COMMANDS for container
JOB_SETUP="export PYTHONPATH=${SELECT_DIR}:$PYTHONPATH && cd ${DLLEE_DIR} && source setup_tufts_container.sh && source configure.sh && cd $jobdir"

# loop over files to run
for f in `cat runlist.txt`
do
    dlreco=$f;
    mpid=$(echo ${dlreco%/*}/multipid_out* | sed 's/\/dlreco/\/stage2\/test21_mpid_pytorch_0121_dlreco/g');
    #echo BOTH ${dlreco} ${mpid};
    COMMAND="python ${SELECT_DIR}/MakeFinalVertexFiles-prime-dlmerged.py ${dlreco} ${calibmap} ${mpid} ${showerreco} ${outdir}"
    echo $COMMAND
    singularity exec $CONTAINER bash -c "$JOB_SETUP && $COMMAND"
done

cp $outdir/*.root $FINAL_OUT_DIR/

echo "FIN"