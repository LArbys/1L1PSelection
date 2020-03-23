#!/bin/bash

setup ubdlana develop -q e17:prof
rootcp out_showerrecov2.root:*ssnetshowerrecov2* temp.root
echo temp.root >> dlanalyze_input_list.txt
echo hist.root >> dlanalyze_input_list.txt
hadd -f merged_dlana.root @dlanalyze_input_list.txt
rm temp.root
