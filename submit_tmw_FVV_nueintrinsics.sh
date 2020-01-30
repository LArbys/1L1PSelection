#!/bin/bash
#
#SBATCH --job-name=FVV-tmw
#SBATCH --output=FVV-tmw.txt
#SBATCH --time=120:00
#SBATCH --mem-per-cpu=2000
#SBATCH --array=104

source /cluster/tufts/wongjiradlab/twongj01/1L1PSelection/run_tmw_FVV_nueintrinsics.sh