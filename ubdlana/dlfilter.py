#! /usr/bin/env python
###############################################################################
#
# Name: dlfilter.py
#
# Purpose: Filter Events from DLANA files
#
# Created: 15-Apr-2020, T. Wongjirad
#
###############################################################################
import sys, os, array
from root_analyze import RootAnalyze
from ROOT import TFile

# Prevent root from printing garbage on initialization.
if os.environ.has_key('TERM'):
    del os.environ['TERM']

# Hide command line arguments from ROOT module.
myargv = sys.argv
sys.argv = myargv[0:1]

# IMPORT ROOT
import ROOT
#ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# IMPORT DL STACK
from larlite import larlite,larutil
from larcv import larcv
from larlitecv import larlitecv

from LEEPreCuts_Functions import makePMTpars,performPMTPrecuts,getPMTPrecutDict

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

    obj = DLFilter(config)
    return obj

# Analyze hit class

class DLFilter(RootAnalyze):

    def __init__(self, config):
        #----------------------------------------------------------------------
        #
        # Purpose: Constructor.
        #
        #----------------------------------------------------------------------

        self.input_file_list = [] # gets append to during open_input function
        self.filter_pars = config['modules']['dlfilter']
        self.filter_type = self.filter_pars['filter_type']
        self.sample_type = self.filter_pars['sample_type']

        DEFINED_FILTERS = ["numu-sideband",
                           "1e1p-highE-sideband",
                           "1e1p-lowBDT-sideband",
                           "1e1p-signal"]

        if self.fitler_type not in DEFINED_FILTERS:
            raise ValueError("Invalid filter type specified [{}]. Defined: {}".format(self.filter_type,DEFINED_FILTERS))

        if bool(self.filter_pars['rerun_precuts']):
            self.rerun_pmtprecuts = True
            self.precutpars = makePMTpars( self.sample_type )
            self.precutpars['ophittree'] = self.filter_pars['precut_ophits']
        else:
            self.rerun_pmtprecuts = False

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

    def event_info(self, tree):
        """
        tree that is passed is suppose to be 'larlite_id_tree'
        """
        return tree._run_id,tree._subrun_id,tree._event_id

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
        print "open_output"
        self.output_file = output_file

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
        """ this is a loop over vertex entries """
        entry = tree.GetReadEntry()

        print "----- ANALYZE ENTRY [%d] ---------------------------------"%(entry)
        #for branch in self.tree_obj.GetListOfBranches():
        #    if self.tree_obj.GetBranchStatus(branch.GetName()):
        #        print '  %s' % branch.GetName()


    def open_input(self, input_file):
        #----------------------------------------------------------------------
        #
        # Purpose: Called when a new input file is opened.
        #
        # Arguments: input_file - Recently opened TFile.
        #
        #----------------------------------------------------------------------

        print 'DLFilter::open_input called.'
        print 'Input file: %s' % input_file.GetName()
        #input_file.ls()
        self.input_file_list.append(input_file.GetName())

        # rerun precuts
        if self.rerun_pmtprecuts:
            self.PMTPrecut_Dict = performPMTPrecuts( input_file.GetName(), **self.precutpars )

        # we need to get the event and vertex trees
        larlite_id_tree = input_file.Get("larlite_id_tree") # provides event-indexing tree
        finalvertextree = input_file.Get("dlana/FinalVertexVariables")
        nevent_entries  = larlite_id_tree.GetEntries()
        nvertex_entries = finalvertextree.GetEntries()

        self.vertex_indexed_trees = []
        self.event_indexed_trees  = []
        self.other_trees = []

        fvv_tree = None

        dirlist = [None]

        while len(dirlist)>0:
            dirname = dirlist.pop(0)
            tdir = input_file
            if dirname is not None:
                tdir = input_file.Get(dirname)
            tdir.cd()
            tree_keys = tdir.GetListOfKeys()
            nkeys = tree_keys.GetEntries()
            for ikey in xrange(nkeys):
                treename = tree_keys.At(ikey).GetName()
                if dirname is not None:
                    treename = dirname+"/"+treename
                atree = input_file.Get(treename)
                if atree is not None and atree.ClassName()=="TTree":
                    print "Tree: ",treename," ",atree.ClassName()
                    nentries = atree.GetEntries()
                    if nentries==nevent_entries:
                        self.event_indexed_trees.append(atree)
                    elif nentries==nvertex_entries:
                        self.vertex_indexed_trees.append(atree)
                        if treename=="dlana/FinalVertexVariables":
                            # save pointer to this tree in case we want to modify it
                            fvv_tree = atree
                    else:
                        print "A tree that does not match either event or vertex indexed entries: ",treename
                        self.other_trees.append(atree)
                elif atree is not None and atree.ClassName()=="TDirectoryFile":
                    dirlist.append(treename)
                else:
                    print "unrecognized: ",atree.ClassName()

            print "directories remaining: ",len(dirlist)
            #raw_input()


        if self.filter_type=="numu-sideband":
            self.run_numu_filter(finalvertextree)
        elif self.filter_type=="1e1p-highE-sideband":
            self.run_1e1p_highE_filter(finalvertextree)
        elif self.filter_type=="1e1p-lowBDT-sideband":
            self.run_1e1p_lowBDT_filder(finalvertextree)
        elif self.filter_type=="1e1p-signal":
            self.run_1e1p_signal_filter(finalvertextree)
        else:
            raise ValueError("unrecognized filter type: ",self.filter_type)

        print "Cloning output trees"

        self.out_event_indexed_trees  = []
        self.out_vertex_indexed_trees = []

        self.output_file.cd()
        outfvv_tree = None
        for tree in self.event_indexed_trees:
            self.out_event_indexed_trees.append( tree.CloneTree(0) )
        for tree in self.vertex_indexed_trees:
            self.out_vertex_indexed_trees.append( tree.CloneTree(0) )
            if tree==fvv_tree:
                outfvv_tree = self.out_vertex_indexed_trees[-1]
        

        neventsout = 0
        nverticesout = 0
        for ientry in xrange( larlite_id_tree.GetEntries() ):
            larlite_id_tree.GetEntry(ientry)
            rse = ( larlite_id_tree._run_id, larlite_id_tree._subrun_id, larlite_id_tree._event_id )
            if rse in self.rse_dict and self.rse_dict[rse]:
                neventsout+=1
                for tree in self.event_indexed_trees:
                    tree.GetEntry(ientry)
                for tree in self.out_event_indexed_trees:
                    tree.Fill()

        # if rerunning steps, we have to replace the branch addresses with new ones
        if self.rerun_pmtprecuts:
            rerun_totpe      = array.array('f',[0.0])
            rerun_maxpefrac  = array.array('f',[0.0])
            rerun_porchtotpe = array.array('f',[0.0])
            rerun_pass       = array.array('i',[0])
            outfvv_tree.SetBranchAddress( "TotPE", rerun_totpe )
            outfvv_tree.SetBranchAddress( "PorchTotPE", rerun_porchtotpe )
            outfvv_tree.SetBranchAddress( "MaxPEFrac", rerun_maxpefrac )
            outfvv_tree.SetBranchAddress( "PassPMTPrecut", rerun_pass )


        for ientry in xrange( finalvertextree.GetEntries() ):
            finalvertextree.GetEntry(ientry)
            rse  = ( finalvertextree.run, finalvertextree.subrun, finalvertextree.event)            
            if rse in self.rse_dict and self.rse_dict[rse]:
                nverticesout+=1

                for tree in self.vertex_indexed_trees:
                    tree.GetEntry(ientry)

                if self.rerun_pmtprecuts:
                    print 'replacing precut results with those from rerun of rse[',rse,']: '
                    print ' totpe old=',fvv_tree.TotPE,' new=',self.PMTPrecut_Dict[rse]['_totpe']
                    print ' maxpefrac old=',fvv_tree.MaxPEFrac,' new=',self.PMTPrecut_Dict[rse]['_maxpefrac']
                    rerun_totpe[0]      = self.PMTPrecut_Dict[rse]['_totpe']
                    rerun_porchtotpe[0] = self.PMTPrecut_Dict[rse]['_porchtotpe']
                    rerun_maxpefrac[0]  = self.PMTPrecut_Dict[rse]['_maxpefrac']
                    rerun_pass[0]       = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0

                for tree in self.out_vertex_indexed_trees:

                    tree.Fill()

        print "Num of event-indexed trees: ",len(self.event_indexed_trees)
        print "Num of vertex-indexed trees: ",len(self.vertex_indexed_trees)
        print "Num of events saved: ",neventsout
        print "Num of vertices saved: ",nverticesout
        print "[ End input tree prep ]"
        print "================================"



    def end_job(self):
        """ end of job tasks """
        print "[End of Job]"
        #self.output_file.cd()
        self.output_file.ls()

        #for tree in self.out_event_indexed_trees:
        #    if tree is not None:
        #        tree.Write()
        #for tree in self.out_vertex_indexed_trees:
        #    if tree is not None:
        #        tree.Write()

    def run_numu_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "run numu filter"
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if ( passprecuts==1
                 and dlanatree.PassSimpleCuts==1
                 and dlanatree.MaxShrFrac<0.2
                 and dlanatree.OpenAng>0.5
                 and dlanatree.ChargeNearTrunk>0
                 and dlanatree.FailedBoost!=1
                 and dlanatree.Lepton_EdgeDist>15.0
                 and dlanatree.Proton_EdgeDist>15.0
                 and dlanatree.BDTscore_1mu1p_cosmic>0.0
                 and dlanatree.BDTscore_1mu1p_nu>0.0 ):
                passes = True
            
            # for debug: make something pass in order to check
            if True:
                passes = True # for debug
                
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes
            self.rsev_dict[rsev] = passes

        # for debug only
        #rsekeys = self.rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",self.rsev_dict[k]

        return

    def run_1e1p_highE_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "run numu filter"
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if ( passprecuts==1
                 and dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 ,
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2
                 and dlanatree.BDTscore_1e1p>0.8 
                 and dlanatree.Enu_1e1p>700.0 ):
                passes = True
            
            # for debug: make something pass in order to check
            if False:
                passes = True # for debug
                
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes
            self.rsev_dict[rsev] = passes

        return

    def run_1e1p_lowBDT_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "run numu filter"
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if ( passprecuts==1
                 and dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 ,
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2
                 and dlanatree.BDTscore_1e1p<=0.8 ):
                passes = True
            
            # for debug: make something pass in order to check
            if True:
                passes = True # for debug
                
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes
            self.rsev_dict[rsev] = passes

        # for debug only
        #rsekeys = self.rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",self.rsev_dict[k]

        return

    def run_1e1p_signal_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "run numu filter"
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if ( passprecuts==1
                 and dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 ,
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2
                 and dlanatree.BDTscore_1e1p>0.8 ):
                passes = True
            
            # for debug: make something pass in order to check
            if True:
                passes = True # for debug
                
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes
            self.rsev_dict[rsev] = passes

        # for debug only
        #rsekeys = self.rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",self.rsev_dict[k]

        return
