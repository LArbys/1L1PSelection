#!/usr/bin/env python
import os,sys,argparse,json
from sys import argv

# ------------------------------------------------------------------------------- #
# PARSE ARGUMENTS
# need to parse arguments before loading ROOT, so args not passed to it's parser
# ------------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(description="Run 1L1P FVV Compiler")
parser.add_argument('-dl', '--dlmerged',required=True,type=str,help="dlmerged file [REQUIRED]")
parser.add_argument('-c','--calibmap',required=True,type=str,help="calibration map file [REQUIRED]")
parser.add_argument('-t','--sample-type',required=True,type=str,help="Sample type. choices: BNB,EXT,Overlay")
parser.add_argument('-f','--filter-type',required=True,type=str,help="Filter type")
parser.add_argument('-o','--outfile',default=".",type=str,help="Output file name")
parser.add_argument('--ismc',help='are we looking at montecarlo?',action="store_true")
parser.add_argument('-pmt','--run-precuts',action='store_true',default=False,help="if true, will run precut code on ophits in file")
parser.add_argument('-oh','--ophits',type=str,default="ophitBeamCalib",help="tree name to use if running PMT precuts. [default: ophitBeam]")
#parser.add_argument('-se','--start-entry',type=int,default=0,help="starting entry")

args = parser.parse_args()

# ------------------------------------------------------------------------------- #
# IMPORT MODULES
#from root_analyze import RootAnalyze
import ubdlana
from ubdlana.dlfilter import DLFilter
import ROOT as rt

ismc = args.ismc

# need to make the configuration dictionary
dlfilter_cfg = {"sample_type":args.sample_type,
                "rerun_precuts":args.run_precuts,
                "precut_ophits":args.ophits,
                "filter_type":args.filter_type}

config = {"modules":{"dlfilter":dlfilter_cfg}}


print "Load module"
dlfilter = DLFilter(config)

# Load input file
tfinput = rt.TFile( args.dlmerged, "open" )

# make output file
tfout = rt.TFile(args.outfile, "recreate" )
dlfilter.open_output( tfout )

# get tree
tree = tfinput.Get( "dlana/FinalVertexVariables" )
dlfilter.open_input(tfinput)

#ientry = 0
#bytes_read = tree.GetEntry(ientry)
#while bytes_read>0:
#    dlfilter.analyze_entry(tree)
#    ientry+=1 
#    bytes_read = tree.GetEntry(ientry)

dlfilter.end_job()

tfout.Write()

print "Finished"

