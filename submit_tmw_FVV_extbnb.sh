#!/bin/bash
#
#SBATCH --job-name=FVV-tmw-extbnb
#SBATCH --output=FVV-tmw-extbnb.txt
#SBATCH --time=120:00
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-132

source /cluster/tufts/wongjiradlab/twongj01/1L1PSelection/run_tmw_FVV_extbnb.sh