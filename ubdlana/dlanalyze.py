#! /usr/bin/env python
###############################################################################
#
# Name: dlanalyze.py
# 
# Purpose: Analyze dlreco files.
#
# Created: 19-Mar-2020, H. Greenlee
#
###############################################################################
import sys, os
from root_analyze import RootAnalyze

# Prevent root from printing garbage on initialization.
if os.environ.has_key('TERM'):
    del os.environ['TERM']

# Hide command line arguments from ROOT module.
myargv = sys.argv
sys.argv = myargv[0:1]

import ROOT
#ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# MPID


# DL Final Vertex Variables
from dlanatree import DLanaTree
import mpidutil
from larlite import larlite
from larcv import larcv
from larlitecv import larlitecv

def make(config):
    #----------------------------------------------------------------------
    #
    # Purpose: Factory function.
    #
    # Arguments: config - FCL configuration.
    #
    # Returns: Instance of class AnalyzeHits.
    #
    #----------------------------------------------------------------------

    obj = DLAnalyze(config)
    return obj

# Analyze hit class

class DLAnalyze(RootAnalyze):

    def __init__(self, config):
        #----------------------------------------------------------------------
        #
        # Purpose: Constructor.
        #
        #----------------------------------------------------------------------

        another_tree = config['modules']['dlanalyze']['another_tree']
        print 'DLAnalyze constructed with second tree = %s' % another_tree
        self.tree_name = "larlite_id_tree"
        self.tree_obj = None

        self.llout_name = config['modules']['dlanalyze']['showerreco']['out_larlite_tree']
        shr_ana         = config['modules']['dlanalyze']['showerreco']['out_ana_tree']
        self.adc_tree   = config['modules']['dlanalyze']['showerreco']['adctree']
        use_calib       = config['modules']['dlanalyze']['showerreco']['use_calib']
        uselarlite      = True
        self.showerreco = larlitecv.ssnetshowerreco.SSNetShowerReco(uselarlite,self.llout_name)
        self.showerreco.set_adc_treename( self.adc_tree )
        if use_calib:
            self.showerreco.use_calibrated_pixsum2mev( True )
        if True:
            self.showerreco.use_second_shower( True )
        if False:
            self.showerreco.use_ncpi0( True )
        if False:
            self.showerreco.use_nueint( True )
        if False:
            self.showerreco.use_bnb( True )
        
        self.showerreco.set_output_treename( shr_ana )

        # list larlite data types
        #self.larlite_types = []
        #for i in xrange(larlite.kDATA_TYPE_MAX):
        #    self.larlite_types.append( larlite.kDATA_TREE_NAME[i] )

        return


    def branches(self):
        #----------------------------------------------------------------------
        #
        # Purpose: Return list of branches we want read for this tree.
        #
        # Returns: List of hit-related branches.
        #
        #----------------------------------------------------------------------

        return ['*']


    def open_output(self, output_file):
        #----------------------------------------------------------------------
        #
        # Purpose: Add content to output file.  Add hit-related histograms.
        #          Called by framework.
        #
        # Arguments: output_file - Open TFile.
        #
        #----------------------------------------------------------------------

        # Make output directory.
        dir = output_file.mkdir('dlana')
        dir.cd()

        # Done.
        self.anatreeclass = DLanaTree()
        self.output_file = output_file

        # Make MPID tree
        self.mpid_data, self.mpid_anatree = mpidutil.make_mpid_anatree(output_file)

        # Make Shower Reco Ana tree
        shr_dir = output_file.mkdir("ssnetshowerreco")
        shr_dir.cd()
        self.showerreco.setupAnaTree()
        self.shr_anatree = self.showerreco.getAnaTree()

        return


    def getLeaf(self, tree, branch_name):
        #----------------------------------------------------------------------
        #
        # Purpose: Utility function to return leaf information for a particular
        #          branch.
        #
        # Arguments: branch name.
        #
        # Returns: Leaf (TLeaf).
        #
        # This function assumes one leaf/branch (true in case of analysis tree).
        # The returned value is an instance of class TLeaf.  To get numeric 
        # values, call TLeaf function GetValue(i), where i=array index or 0 for 
        # scalar leaf.
        #
        #----------------------------------------------------------------------

        result = None
        br = tree.GetBranch(branch_name)
        leaves = br.GetListOfLeaves()
        if len(leaves) > 0:
            result = leaves[0]
        return result


    def analyze_entry(self, tree):
        #----------------------------------------------------------------------
        #
        # Purpose: Analyze loaded tree (fill histograms).  Called by framework.
        #
        # Arguments: tree - Loaded tree.
        #
        #----------------------------------------------------------------------

        entry = tree.GetReadEntry()
        print 'DLReco::analyze_entry called for entry %d.' % entry

        # Get a leaf from the main tree.

        rse_br = []
        for brname in ['_run_id','_subrun_id','_event_id']:
            br = tree.GetBranch(brname)
            leaves = br.GetListOfLeaves()
            print '%s = %d' % (brname,int(leaves[0].GetValue()))
            rse_br.append( int(leaves[0].GetValue()) )

        # load larcv iomanager entry (larlite already loaded)
        self.in_lcv.read_entry( entry )

        # get larcv data
        ev_adc    = self.in_lcv.get_data(larcv.kProductImage2D,self.adc_tree)
        lcv_rse = [ev_adc.run(),ev_adc.subrun(),ev_adc.event()]
        ll_rse  = [self.io_ll.run_id(),self.io_ll.subrun_id(),self.io_ll.event_id()]
        if rse_br != lcv_rse or rse_br != ll_rse:
            raise RuntimeError("(run,subrun,event) for event loop tree[{}], larlite tree[{}], larcv tree[{}] do not match".format(rse_br,ll_rse,lcv_rse))
        
        # check run, subrun, event
        print "process: ",rse_br

        # run shower reco
        self.showerreco.process( self.in_lcv, self.io_ll, entry )
        self.showerreco.store_in_larlite(self.io_ll)
        self.io_ll.next_event()


    def open_input(self, input_file):
        #----------------------------------------------------------------------
        #
        # Purpose: Called when a new input file is opened.
        #
        # Arguments: input_file - Recently opened TFile.
        #
        #----------------------------------------------------------------------

        print 'DLReco::open_input called.'
        print 'Input file: %s' % input_file.GetName()
        input_file.ls()

        # we open the larcv and larlite iomanagers
        self.io_ll  = larlite.storage_manager(larlite.storage_manager.kBOTH)
        self.io_ll.add_in_filename( input_file.GetName() )
        self.io_ll.set_out_filename( "out_showerrecov2.root" )
        self.io_ll.set_data_to_read( larlite.data.kTrack, "trackReco" )
        #self.io_ll.set_data_to_write( larlite.data.kShower, self.llout_name )
        #self.io_ll.set_data_to_write( larlite.data.kShower, self.llout_name+"_sec" )
        #self.io_ll.set_data_to_write( larlite.data.kLArFlowCluster, self.llout_name )
        self.io_ll.open()
        self.io_ll.next_event() # go to first entry

        self.in_lcv = larcv.IOManager(larcv.IOManager.kREAD,"input_larcv")
        self.in_lcv.add_in_file( input_file.GetName() )
        self.in_lcv.initialize()

        # Use the larlite index tree as the index tree
        print 'Looking for TTree named %s' % self.tree_name
        obj = input_file.Get(self.tree_name)
        if obj.InheritsFrom('TTree'):
            print 'Found %s with %d entries.' % (self.tree_name, obj.GetEntriesFast())

            # Activate all branches.

            obj.SetBranchStatus('*', 1)
            print 'List of activated branches of tree %s' % self.tree_name
            for branch in obj.GetListOfBranches():
                if obj.GetBranchStatus(branch.GetName()):
                    print '  %s' % branch.GetName()

            self.tree_obj = obj
            print




    def end_job(self):
        """ close larcv and larlite files. larlite output file will write."""
        self.in_lcv.finalize()
        self.io_ll.close()

