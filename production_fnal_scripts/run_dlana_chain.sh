#!/bin/bash


##################################################################
# INTERACTIVE OPTIONS
# variables script assumes has been set by CONDOR

_CONDOR_SCRATCH_DIR=/tmp/
CONDOR_DIR_INPUT=/pnfs/uboone/resilient/users/tmw/mpid_tests
CLUSTER=1000
PROCESS=0

###################
#### ARGS #########

INFILE_LIST=$CONDOR_DIR_INPUT/merged_flist_test.dat
LOGDIR=/pnfs/uboone/persistent/users/tmw/dlana_test/log
OUTDIR=/pnfs/uboone/persistent/users/tmw/dlana_test/out
IFDH_OPT=""


##################################################################
### SETUP WORKING ENVIRONMENT, COPIED FROM WRAP script by Herb ###
# Create the scratch directory in the condor scratch diretory.
# Copied from condor_lBdetMC.sh.
# Scratch directory path is stored in $TMP.
# Scratch directory is automatically deleted when shell exits.

# Make sure scratch directory is defined.
# For batch, the scratch directory is always $_CONDOR_SCRATCH_DIR
# For interactive, the scratch directory is specified by option 
# --scratch or --outdir.

## FROM CONDOR DOCUMENTATION
## _CONDOR_SCRATCH_DIR gives the directory where the job may place temporary data files. 
## This directory is unique for every job that is run, 
## and it's contents are deleted by Condor when the job stops running on a machine, 
## no matter how the job completes.

SCRATCH=$_CONDOR_SCRATCH_DIR


# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
TMP=`mktemp -d ${SCRATCH}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${SCRATCH}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { rm -rf \"$TMP\"; }" 0
chmod 755 $TMP
cd $TMP
# End of the section you should not change.

echo "Scratch directory: $TMP"

# Copy files from work directory to scratch directory.

echo "No longer fetching files from work directory."
echo "that's now done with using jobsub -f commands"
mkdir work
cp ${CONDOR_DIR_INPUT}/* ./work/
cd work
mkdir out
mkdir log
find . -name \*.tar -exec tar xf {} \;
find . -name \*.py -exec chmod +x {} \;
find . -name \*.sh -exec chmod +x {} \;
echo "Local working directoroy:"
pwd
ls
echo

# Save the hostname and condor job id.
hostname > hostname.txt
echo ${CLUSTER}.${PROCESS} > jobid.txt

# Set default CLUSTER and PROCESS environment variables for interactive jobs.
echo "Cluster: $CLUSTER"
echo "Process: $PROCESS"

# Construct name of output subdirectory.

OUTPUT_SUBDIR=${CLUSTER}_${PROCESS}
echo "Output subdirectory: $OUTPUT_SUBDIR"
### END OF SETUP
##########################################################################

#######################
##### DO WORK #########

# SETUP UBOONE CVMFS
#source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh

# SETUP DLLEE_UNIFIED,MPID UPS
#setup ifdhc v2_4_1 -q e17:p2714b:prof
#setup ubMPIDnet v1_0_0 -q NULL
#setup dllee_unified develop -q e17:prof
#setup ubdlana develop -q e17:prof

# DETERMINE INPUT DLMERGED FILE
let lineno=${PROCESS}+1
inputfile=`sed -n ${lineno}p $INFILE_LIST`
echo $inputfile > log/infile.txt
inputfile_xrootd=`sed 's|/pnfs|xroot://fndca1.fnal.gov/pnfs/fnal.gov/usr|g' log/infile.txt`

echo "INFILE_LIST: "$INFILE_LIST
echo "OUTDIR: "$OUTDIR
echo "INPUTFILE="${inputfile}
echo "INPUTFILE (XROOTD): "${inputfile_xrootd}

### MPID ###
echo "<<<<<<< RUN MPID >>>>>>>"

# copy script to log folder
cp ${MPIDMODEL_DIR}/production_cfg/inference_config_tufts_WC.cfg log/

# execute
echo "run inference_pid_torch_dlmerger_WC.py"
inference_pid_torch_dlmerger_WC.py ${inputfile_xrootd} out log/inference_config_tufts_WC.cfg >& log/mpid.out

### SHOWER ###
echo "<<<<<<< RUN SHOWER RECO >>>>>>>"
run_ssnetshowerreco.py --input-larcv ${inputfile_xrootd} --input-larlite ${inputfile_xrootd} -f both -o out_showerreco -sec >& log/showerreco.out
cp out_showerreco* out/

### SELECTION VARIABLE TREES ###
#make_

### END OF WORK #####################
#####################################

#####################################
#### COPY OUTPUT ####################
# Make a tarball of the log directory contents, and save the tarball in the log directory.

rm -f log.tar
tar -cjf log.tar -C log .
mv log.tar log

# Create remote output and log directories.
export IFDH_CP_MAXRETRIES=5

echo "Make directory ${LOGDIR}/${OUTPUT_SUBDIR}."
date
ifdh mkdir $IFDH_OPT ${LOGDIR}/$OUTPUT_SUBDIR
echo "Done making directory ${LOGDIR}/${OUTPUT_SUBDIR}."
date

if [ ${OUTDIR} != ${LOGDIR} ]; then
  echo "Make directory ${OUTDIR}/${OUTPUT_SUBDIR}."
  date
  ifdh mkdir $IFDH_OPT ${OUTDIR}/$OUTPUT_SUBDIR
  echo "Done making directory ${OUTDIR}/${OUTPUT_SUBDIR}."
  date
fi

# Transfer tarball in log subdirectory.

statout=0
echo "ls log"
ls log
echo "ifdh cp -D $IFDH_OPT log/log.tar ${LOGDIR}/$OUTPUT_SUBDIR"
ifdh cp -D $IFDH_OPT log/log.tar ${LOGDIR}/$OUTPUT_SUBDIR
date
stat=$?
if [ $stat -ne 0 ]; then
  statout=1
  echo "ifdh cp failed with status ${stat}."
fi

# Transfer data files in out subdirectory.
if [ "$( ls -A out )" ]; then
  echo "ifdh cp -D $IFDH_OPT out/* ${OUTDIR}/$OUTPUT_SUBDIR"
  ifdh cp -D $IFDH_OPT out/* ${OUTDIR}/$OUTPUT_SUBDIR
  stat=$?
  if [ $stat -ne 0 ]; then
      statout=1
      echo "ifdh cp failed with status ${stat}."
  fi
fi
