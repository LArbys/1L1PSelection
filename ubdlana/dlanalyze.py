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

        dir = output_file.mkdir('dlreco')
        dir.cd()

        # Done.

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

        # Load entry for second tree.
        if True:
            return

        self.tree_obj.GetEntry(entry)

        # Get leaf from second tree.

        br = self.tree_obj.GetBranch('selected1L1P')
        leaves = br.GetListOfLeaves()
        print 'selected1L1P = %d' % int(leaves[0].GetValue())



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






