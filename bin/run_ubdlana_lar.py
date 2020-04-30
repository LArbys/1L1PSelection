#!/usr/bin/env python
import os,sys,argparse,json
from sys import argv

# ------------------------------------------------------------------------------- #
# PARSE ARGUMENTS
# need to parse arguments before loading ROOT, so args not passed to it's parser
# ------------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(description="Run 1L1P FVV Compiler")
parser.add_argument('-d','--dlmerged',required=True,type=str,help="dlmerged file [REQUIRED]")
parser.add_argument('-c','--calibmap',required=True,type=str,help="calibration map file [REQUIRED]")
#parser.add_argument('-m','--mpid',required=True,type=str,help="mpid file [REQUIRED]")
parser.add_argument('-t','--sample-type',required=True,type=str,help="Sample type. choices: BNB,EXT,Overlay")
parser.add_argument('-o','--outfile',default=".",type=str,help="Output file name")
parser.add_argument('-g','--tag',default="test",type=str,help="Tag: used to name files. [default: 'test']")
parser.add_argument('--ismc',help='are we looking at montecarlo?',action="store_true")
parser.add_argument('-pmt','--run-precuts',action='store_true',default=False,help="if true, will run precut code on ophits in file")
parser.add_argument('-oh','--ophits',type=str,default="ophitBeamCalib",help="tree name to use if running PMT precuts. [default: ophitBeamCalib]")
parser.add_argument('-se','--start-entry',type=int,default=0,help="starting entry")
parser.add_argument("-cnn","--run-cnn",action='store_true',default=False,help="if true, run shower cnn")
args = parser.parse_args()

# ------------------------------------------------------------------------------- #
# IMPORT MODULES
#from root_analyze import RootAnalyze
import ubdlana
import ubdlana.dlanalyze
from ubdlana.dlanalyze import DLAnalyze
import ROOT as rt

if args.ismc:
    ismc = True
else:
    ismc = False

# need to make the configuration dictionary
dlanalyze_cfg = {"tracker_tree":"_recoTree",
                 "ismc":ismc,
                 "sample_type":args.sample_type,
                 "another_tree":"DLAnaTree",
                 "precut_ophits":args.ophits,
                 "crtveto":{"opflash_producer":"simpleFlashBeam",
                            "crthit_producer":"crthitcorr"},
                 "bdt_1e1p_weights":"bdtweights_1e1p_MissingLowEPatch_4-28-20.pickle",
                 #"bdt_1e1p_weights":"Selection_Weights_1e1p_3-24-20.pickle",
                 "bdt_1mu1p_cosmic_weights":"bdtweights_1mu1p_WC_apr1.pickle",
                 "bdt_1mu1p_nu_weights":"x",
                 "showercnn_weights":"DoNotRun",
                 "showerreco": { "out_larlite_tree": "ssnetshowerrecov2",
                                 "out_ana_tree": "ssnetshowerrecov2ana",
                                 "adctree": "wire",
                                 "second_shower": True,
                                 "use_calib": False }
}

if args.run_cnn:
    dlanalyze_cfg['showercnn_weights'] = "ResNet18_nueintrinsic_JoshImages_epoch120.pytorch_weights"

config = {"modules":{"dlanalyze":dlanalyze_cfg}}


print "Load module"
dlanalyze = DLAnalyze(config)
dlanalyze.set_start_entry(args.start_entry)

# Load input file
tfinput = rt.TFile( args.dlmerged, "open" )

# make output file
tfout = rt.TFile(args.outfile, "recreate" )
dlanalyze.open_output( tfout )

# get tree
tree = tfinput.Get( dlanalyze_cfg['tracker_tree'] )
dlanalyze.open_input(tfinput)


ientry = args.start_entry
bytes_read = tree.GetEntry(ientry)

while bytes_read>0:
    dlanalyze.analyze_entry(tree)
    ientry+=1 
    bytes_read = tree.GetEntry(ientry)

dlanalyze.end_job()

tfout.Write()

print "Finished"

