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
parser.add_argument("--rerun-1mu1p-bdt",action='store_true',default=False,help="if true, will rerun 1mu1p BDT")
parser.add_argument("--bdt-weightfile-1mu1p",type=str,default=None,help="Must be set if rerunning 1mu1p BDT")
parser.add_argument("--use-ubdlana-eventtree",action='store_true',default=False,help="if true, will use ubdlana event tree instead of defaul larlite_id_tree")
parser.add_argument("--rse-list",type=str,default=None,help="Must provide RSE list path if filter-type set to rse-list")
#parser.add_argument('-se','--start-entry',type=int,default=0,help="starting entry")

args = parser.parse_args()

# ------------------------------------------------------------------------------- #
# IMPORT MODULES
#from root_analyze import RootAnalyze
import ubdlana
from ubdlana.dlfilter import DLFilter
import ROOT as rt

ismc = args.ismc

if args.rerun_1mu1p_bdt and args.bdt_weightfile_1mu1p is None:
    raise ValueError("1mu1p BDT is flagged to be rerun, but weight file is not specified")
    

# need to make the configuration dictionary
dlfilter_cfg = {"sample_type":args.sample_type,
                "event_tree":"larlite_id_tree",
                "rerun_precuts":args.run_precuts,
                "precut_ophits":args.ophits,
                "filter_type":args.filter_type}

if args.rerun_1mu1p_bdt:
    print "RERUN 1mu1p BDT"
    dlfilter_cfg["rerun_1mu1p_bdt"] = True
    dlfilter_cfg["bdt_weights_1mu1p"] = args.bdt_weightfile_1mu1p

if args.use_ubdlana_eventtree:
    print "USE UBDLANA ID TREE"
    dlfilter_cfg["event_tree"] = "dlana/ubdlana_id_tree"

if args.filter_type=="rse-list":
    if args.rse_list is None or not os.path.exists(args.rse_list):
        raise ValueError("Must provide valid path to RSE list if running rse-list filter")
    dlfilter_cfg['rse-list'] = args.rse_list

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

