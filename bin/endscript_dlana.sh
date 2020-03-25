#!/bin/bash

ls -lh
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup ubdlana develop -q e17:prof
rootcp out_showerrecov2.root:*ssnetshowerrecov2* temp.root || { echo "error copying larlite shower trees to temp"; exit 1; };
echo temp.root >> dlanalyze_input_list.txt
mv merged_dlana0.root hist.root
mv merged_dlana0.root.json merged_dlana.root.json
ls hist*.root >> dlanalyze_input_list.txt
hadd -f merged_dlana.root @dlanalyze_input_list.txt || { echo "error hadding temp and input files"; exit 1; }
rm temp.root
rm hist.root
rm out_showerrecov2.root
