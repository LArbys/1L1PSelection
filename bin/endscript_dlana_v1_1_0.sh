#!/bin/bash


echo "start end-script for dl-ana"
echo "SAM SCHEMA: "$SAM_SCHEMA
echo "local files:"
ls -lh
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup ubdlana v1_1_0 -q e17:prof
rootcp out_showerrecov2.root:*ssnetshowerrecov2* temp1.root || { echo "error copying larlite shower trees to temp"; exit 1; };
rootcp out_showerrecov2.root:*dqdx_* temp2.root || { echo "error copying larlite dqdx trees to temp"; exit 1; };

input=transferred_uris.list

if [ x$SAM_SCHEMA = xroot ]; then
    # if using xrootd
    echo "use xrootd for merging"
else
    # not using xrootd, retransfer input files
    while IFS= read -r line
    do
	echo "copy $line"
	fname=`basename $line`
	ifdh cp $line $fname
    done < "$input"
    echo "after (re-)transfer"
    ls -lh
fi

# add extracted shower code
echo temp1.root >> dlanalyze_input_list.txt
echo temp2.root >> dlanalyze_input_list.txt
# change name of output ana file and its json
mv merged_dlana0.root hist.root
mv merged_dlana0.root.json merged_dlana.root.json
# add ana file to hadd list
ls hist*.root >> dlanalyze_input_list.txt

echo "files to merge: "
cat dlanalyze_input_list.txt
echo ""

echo "attempt hadd"
hadd -f merged_dlana.root @dlanalyze_input_list.txt || { echo "error hadding temp and input files"; exit 1; };

echo "clean up"
rm temp1.root
rm temp2.root
rm hist.root
rm out_showerrecov2.root

if [ x$SAM_SCHEMA = xroot ]; then
    # if using xrootd
    echo "use xrootd for merging, no files to clean."
else
    while IFS= read -r line
    do
	fname=`basename $line`
	echo "clean-up $fname"
	rm $fname
    done < "$input"
    echo "after cleaning"
    ls -lh
fi

