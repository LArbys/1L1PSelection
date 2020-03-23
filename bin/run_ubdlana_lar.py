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
parser.add_argument('-m','--mpid',required=True,type=str,help="mpid file [REQUIRED]")
parser.add_argument('-s','--showerreco',required=True,type=str,help="shower reco file [REQUIRED]")
parser.add_argument('-t','--sample-type',required=True,type=str,help="Sample type. choices: BNB,EXT,Overlay")
parser.add_argument('-o','--outdir',default=".",type=str,help="Output directory")
parser.add_argument('-g','--tag',default="test",type=str,help="Tag: used to name files. [default: 'test']")
parser.add_argument('--ismc',help='are we looking at montecarlo?',action="store_true")
parser.add_argument('-pmt','--run-precuts',action='store_true',default=False,help="if true, will run precut code on ophits in file")
parser.add_argument('-oh','--ophits',type=str,default="ophitBeam",help="tree name to use if running PMT precuts. [default: ophitBeam]")
args = parser.parse_args()

# ------------------------------------------------------------------------------- #
# IMPORT MODULES
from root_analyze import RootAnalyze
import ubdlana
import ubdlana.dlanalyze
from ubdlana.dlanalyze import DLAnalyze

print DLAnalyze

