#!/bin/bash
#
#SBATCH --job-name=FVV-tmw-overlay
#SBATCH --output=FVV-tmw-overlay.txt
#SBATCH --time=240:00
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-252

source /cluster/tufts/wongjiradlab/twongj01/1L1PSelection/run_tmw_FVV_overlay.sh