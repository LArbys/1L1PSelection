#!/bin/bash
#
#SBATCH --job-name=FVV-tmw-bnb5e19
#SBATCH --output=FVV-tmw-bnb5e19.txt
#SBATCH --time=120:00
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-137

source /cluster/tufts/wongjiradlab/twongj01/1L1PSelection/run_tmw_FVV_bnb5e19.sh