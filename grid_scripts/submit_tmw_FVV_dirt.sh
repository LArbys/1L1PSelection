#!/bin/bash
#
#SBATCH --job-name=FVV-tmw-dirt
#SBATCH --output=FVV-tmw-dirt.txt
#SBATCH --time=240:00
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-25

source /cluster/tufts/wongjiradlab/twongj01/1L1PSelection/run_tmw_FVV_dirt.sh