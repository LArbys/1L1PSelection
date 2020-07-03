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
parser.add_argument('-t','--sample-type',required=True,type=str,help="Sample type. choices: BNB,EXT,Overlay")
parser.add_argument('-o','--outfile',default=".",type=str,help="Output file name")
parser.add_argument('-g','--tag',default="test",type=str,help="Tag: used to name files. [default: 'test']")
parser.add_argument('--ismc',default=False,help='are we looking at montecarlo?',action="store_true")
parser.add_argument('-pmt','--run-precuts',action='store_true',default=False,help="if true, will run precut code on ophits in file")
parser.add_argument('-oh','--ophits',type=str,default="ophitBeamCalib",help="tree name to use if running PMT precuts. [default: ophitBeamCalib]")
parser.add_argument('-se','--start-entry',type=int,default=0,help="starting entry")
parser.add_argument("-cnn","--run-cnn",action='store_true',default=False,help="if true, run shower cnn")
parser.add_argument("--bdt-1m1p",default=None,type=str,help="Specify 1m1p BDT weight file [if not using, default]")
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

shower_pars = {("BNB",False):[-69.049,0.112],# needs update 6/22/20
               ("EXT",False):[-69.049,0.112],# needs udpate 6/22/20
               ("Overlay",True):[-69.049,0.112]}

if (args.sample_type,ismc) not in shower_pars:
    raise ValueError("Sample+MC flag combination not valid: sample={} ismc={}".format( args.sampletype,ismc ) )


# need to make the configuration dictionary
dlanalyze_cfg = {"tracker_tree":"_recoTree",
                 "ismc":ismc,
                 "sample_type":args.sample_type,
                 "another_tree":"DLAnaTree",
                 "precut_ophits":args.ophits,
                 "crtveto":{"opflash_producer":"simpleFlashBeam",
                            "crthit_producer":"crthitcorr"},
                 "bdt_1e1p_weights":"BDT_1e1p_Run1_6-2-20.pickle",
                 "bdt_1mu1p_weights":"bdtweight_series2_june1_run1.pickle",
                 "showercnn_weights":"DoNotRun",
                 "showerreco": { "out_larlite_tree": "ssnetshowerrecov2",
                                 "out_ana_tree": "ssnetshowerrecov2ana",
                                 "adctree": "wire",
                                 "second_shower": True,
                                 "use_calib": False,
                                 "pix2energy_params":shower_pars[(args.sample_type,ismc)]}
}

if args.run_cnn:
    dlanalyze_cfg['showercnn_weights'] = "ResNet18_nueintrinsic_JoshImages_epoch120.pytorch_weights"

if args.bdt_1m1p is not None:
    if not os.path.exists(args.bdt_1m1p):
        raise ValueError("could not find 1m1p BDT weightfile: ",args.bdt_1m1p)
    dlanalyze_cfg['bdt_1mu1p_weights'] = args.bdt_1m1p

config = {"modules":{"dlanalyze":dlanalyze_cfg}}


print "Load module"
dlanalyze = DLAnalyze(config)
dlanalyze.set_start_entry(args.start_entry)

# Load input file
tfinput = rt.TFile( args.dlmerged, "open" )

# make output file
tfout = rt.TFile(args.outfile, "recreate" )
dlanalyze.open_output( tfout )

# get tracker tree
tracker_tree = tfinput.Get( dlanalyze_cfg['tracker_tree'] )
dlanalyze.open_input(tfinput)

# get the event tree
event_tree = tfinput.Get( 'larlite_id_tree' )

ientry = args.start_entry
bytes_read = event_tree.GetEntry(ientry)

while bytes_read>0:
    dlanalyze.analyze_entry(event_tree)
    ientry+=1 
    bytes_read = event_tree.GetEntry(ientry)

dlanalyze.close_input(tfinput)

tfout.Write()

print "Finished"

