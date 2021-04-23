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
try:    
    from root_analyze import RootAnalyze
except:
    """ dumpy class """
    print "No root_analyze. Making a dummy class"
    class RootAnalyze:
        def __init__(self):
            pass
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
from event_indexed_trees_util import is_tree_event_indexed
import bdtutil

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
        self.DATARUN     = self.filter_pars['data_run']
        
        if 'event_tree' in self.filter_pars:
            self.event_tree  = self.filter_pars['event_tree']
        else:
            self.event_tree  = "dlana/ubdlana_id_tree"

        DEFINED_FILTERS = ["numu-sideband",
                           "pi0-lowBDT-sideband",
                           "1e1p-highE-sideband",
                           "1e1p-nearE-sideband",
                           "1e1p-lowBDT-sideband",
                           "1e1p-midBDT-sideband",                           
                           "1e1p-signal",
                           "rse-list"]

        if self.filter_type not in DEFINED_FILTERS:
            raise ValueError("Invalid filter type specified [{}]. Defined: {}".format(self.filter_type,DEFINED_FILTERS))

        if bool(self.filter_pars['rerun_precuts']):
            self.rerun_pmtprecuts = True
            self.precutpars = makePMTpars( self.sample_type )
            self.precutpars['ophittree'] = self.filter_pars['precut_ophits']
        else:
            self.rerun_pmtprecuts = False

        if "rerun_1mu1p_bdt" in self.filter_pars and bool(self.filter_pars["rerun_1mu1p_bdt"]):
            self.rerun_1mu1p_bdt = True
            self.bdt_1mu1p_weightfile = self.filter_pars["bdt_weights_1mu1p"]
            if self.filter_pars["bdt_model_1mu1p"]=="single":
                self.bdt_model_1mu1p = bdtutil.load_BDT_model( self.bdt_1mu1p_weightfile )
            else if self.filter_pars["bdt_model_1mu1p"]=="ensemble":
                self.bdt_model_1mu1p = bdtutil.load_BDT_ensemble( "1m1p", self.bdt_1mu1p_weightfile, nbdts=10, runs=[self.DATARUN] )
            print "DLFilter: RERUN 1MU1P BDT"
        else:
            self.rerun_1mu1p_bdt = False

        if "rerun_1e1p_bdt" in self.filter_pars and bool(self.filter_pars["rerun_1e1p_bdt"]):
            self.rerun_1e1p_bdt = True            
            self.bdt_1e1p_weightfile = self.filter_pars["bdt_weights_1e1p"]
            if self.filter_pars["bdt_model_1mu1p"]=="single":
                self.bdt_model_1e1p = bdtutil.load_BDT_model( self.bdt_1e1p_weightfile )
            else if self.filter_pars["bdt_model_1mu1p"]=="ensemble":                
                self.bdt_model_1e1p = bdtutil.load_BDT_ensemble( "1e1p", self.bdt_1e1p_weightfile, nbdts=10, runs=[self.DATARUN] )
            print "DLFilter: RERUN 1e1P BDT"
        else:
            self.rerun_1e1p_bdt = False
            

        self._DEBUG_MODE_ = False
            
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
        if not self._use_ubdlana_idtree:
            rse = (tree._run_id,tree._subrun_id,tree._event_id)
        else:
            rse = (tree.run,tree.subrun,tree.event)
        print "DLFilter::event_info = ",rse
        return rse

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

        # we need to get the event and vertex trees
        larlite_id_tree = input_file.Get(self.event_tree)
        finalvertextree = input_file.Get("dlana/FinalVertexVariables")
        self._use_ubdlana_idtree = False
        if self.event_tree=="dlana/ubdlana_id_tree":
            try:
                print "DLFilter::Using dlana/ubdlana_id_tree"
                nevent_entries  = larlite_id_tree.GetEntries()
                self._use_ubdlana_idtree = True                
            except:
                print "DLFilter::Could not load ubdlana. try larlite_id_tree"
                larlite_id_tree = None
                nevent_entries = 0
                self._use_ubdlana_idtree = False
        if self.event_tree=="larlite_id_tree" or larlite_id_tree is None:
            larlite_id_tree = input_file.Get("larlite_id_tree")
            try:
                nevent_entries = larlite_id_tree.GetEntries()
                print "DLFilter: using larlite_id_tree"
                self._use_ubdlana_idtree = False
            except:
                print "DLFilter::Could not load ubdlana_id_tree either"
                larlite_id_tree = None
                nevent_entries = 0
        print "DLFilter::event id tree: ",self.event_tree
        print "DLFilter::_use_ubdlana_idtree: ",self._use_ubdlana_idtree

        nevent_entries  = larlite_id_tree.GetEntries()
        nvertex_entries = finalvertextree.GetEntries()
            
        self.vertex_indexed_trees = []
        self.event_indexed_trees  = []
        self.other_trees = []

        fvv_tree = None

        dirlist = [None]
        dirdict = {}
        
        self.pot_sum_tree = None
        
        while len(dirlist)>0:
            dirname = dirlist.pop(0)
            tdir = input_file
            if dirname is not None:
                tdir = input_file.Get(dirname)
            tdir.cd()
            tree_keys = tdir.GetListOfKeys()
            nkeys = tree_keys.GetEntries()
            for ikey in xrange(nkeys):
                basetreename = tree_keys.At(ikey).GetName()
                if dirname is not None:
                    treename = dirname+"/"+basetreename
                    dirdict[basetreename] = dirname
                else:
                    treename = basetreename
                    
                atree = input_file.Get(treename)
                if atree is not None and atree.ClassName()=="TTree":
                    print "Tree: ",treename," ",atree.ClassName()
                    nentries = atree.GetEntries()
                    if treename=="potsummary_generator_tree":
                        print "Found potsummary_generator_tree. special case."
                        self.pot_sum_tree = atree
                        continue

                    if is_tree_event_indexed( basetreename ):
                        self.event_indexed_trees.append(atree)
                        if nevent_entries>0 and nentries!=nevent_entries and basetreename!="ubdlana_id_tree":
                            raise RuntimeError("Event-indexed tree({}) nentries ({}) does not seem to match known event-indexed tree entries ({})".format(basetreename,nentries,nevent_entries))
                    else:
                        self.vertex_indexed_trees.append(atree)
                        if nentries!=nvertex_entries:
                            raise RuntimeError("Vertex-indexed tree({}) nentries ({}) does not seem to match known vertex-indexed tree entries ({})".format(basetreename,nentries,nvertex_entries))
                        if treename=="dlana/FinalVertexVariables":
                            # save pointer to this tree in case we want to modify it
                            fvv_tree = atree

                elif atree is not None and atree.ClassName()=="TDirectoryFile":
                    dirlist.append(treename)
                else:
                    print "unrecognized: ",atree.ClassName()

            print "directories remaining: ",len(dirlist)
            #raw_input()

        # FOR DEBUG
        if False:
            print "[EVENT-INDEXED TREES]"
            for atree in self.event_indexed_trees:
                print atree.GetName()
        if True:
            print "[VERTEX-INDEXED TREES]"
            for atree in self.vertex_indexed_trees:
                print atree.GetName()

        # rerun precuts
        if self.rerun_pmtprecuts and nevent_entries>0:
            self.PMTPrecut_Dict = performPMTPrecuts( input_file.GetName(), **self.precutpars )

        # rerun the 1mu1p
        # we get back a dictionary indexed by (run,subrun,event,vertex)
        # which stores replacement values for the BDT variables and BDT score
        if self.rerun_1mu1p_bdt:
            if self.filter_pars["bdt_model_1mu1p"]=="single":            
                self.bdtoutput_1mu1p = bdtutil.rerun_1mu1p_models( self.bdt_model_1mu1p, finalvertextree )
            else:
                self.bdtoutput_1mu1p = bdtutil.rerun_1mu1p_ensemble( self.bdt_model_1mu1p, finalvertextree, self.DATARUN, nbdts=10 )
        else:
            self.bdtoutput_1mu1p = {}

        # rerun the 1e1p
        # we get back a dictionary indexed by (run,subrun,event,vertex)
        # which stores replacement values for the BDT variables and BDT score
        if self.rerun_1e1p_bdt:
            if self.filter_pars["bdt_model_1e1p"]=="single":                        
                self.bdtoutput_1e1p = bdtutil.rerun_1e1p_models( self.bdt_model_1e1p, finalvertextree )
            else:
                self.bdtoutput_1e1p = bdtutil.rerun_1e1p_ensemble( self.bdt_model_1e1p, finalvertextree, self.DATARUN, nbdts=10 )
        else:
            self.bdtoutput_1e1p = {}

        # We run the filter on the finalvariablevertex trees
        # Produces a list of run,subrun,event,vertexid indicating if events pass
        if self.filter_type=="numu-sideband":
            self.run_numu_filter(finalvertextree)
        elif self.filter_type=="1e1p-highE-sideband":
            self.run_1e1p_highE_filter(finalvertextree)
        elif self.filter_type=="1e1p-nearE-sideband":
            self.run_1e1p_nearE_filter(finalvertextree)
        elif self.filter_type=="1e1p-lowBDT-sideband":
            self.run_1e1p_lowBDT_filter(finalvertextree)
        elif self.filter_type=="1e1p-midBDT-sideband":
            self.run_1e1p_midBDT_filter(finalvertextree)
        elif self.filter_type=="1e1p-signal":
            self.run_1e1p_signal_filter(finalvertextree)
        elif self.filter_type=="pi0-lowBDT-sideband":
            self.run_lowBDT_pi0_filter(finalvertextree)
        elif self.filter_type=="rse-list":
            self.run_rse_filter(finalvertextree)
        else:
            raise ValueError("unrecognized filter type: ",self.filter_type)

        print "Cloning output trees"

        self.out_event_indexed_trees  = []
        self.out_vertex_indexed_trees = []

        self.output_file.cd()
        outfvv_tree = None
        for tree in self.event_indexed_trees:
            if str(tree.GetName()) in dirdict:
                try:
                    rootdir = self.output_file.mkdir( dirdict[str(tree.GetName())] )
                except:
                    rootdir = self.output_file.Get( dirdict[str(tree.GetName())] )
                rootdir.cd()
            self.out_event_indexed_trees.append( tree.CloneTree(0) )
            self.output_file.cd()
        for tree in self.vertex_indexed_trees:
            if str(tree.GetName()) in dirdict:
                rootdir = self.output_file.Get( dirdict[str(tree.GetName())] )
                try:
                    rootdir.cd()
                except:
                    rootdir = self.output_file.mkdir( dirdict[str(tree.GetName())] )
                    rootdir.cd()

            self.out_vertex_indexed_trees.append( tree.CloneTree(0) )
            if tree==fvv_tree:
                outfvv_tree = self.out_vertex_indexed_trees[-1]
            self.output_file.cd()
        

        neventsout = 0
        nverticesout = 0
        for ientry in xrange( nevent_entries ):
            larlite_id_tree.GetEntry(ientry)
            if not self._use_ubdlana_idtree:
                rse = ( larlite_id_tree._run_id, larlite_id_tree._subrun_id, larlite_id_tree._event_id )
            else:
                rse = (larlite_id_tree.run,larlite_id_tree.subrun,larlite_id_tree.event)
            if rse in self.rse_dict and self.rse_dict[rse]:
                neventsout+=1
                for tree in self.event_indexed_trees:
                    tree.GetEntry(ientry)
                for tree in self.out_event_indexed_trees:
                    tree.Fill()
            else:
                if rse in self.rse_dict:
                    print "[do not save event] ",rse," ",self.rse_dict[rse]
                else:
                    print "[do not save event] ",rse," not in RSE dict"

        # if rerunning steps, we have to replace the branch addresses with new ones
        if self.rerun_pmtprecuts and nevent_entries>0:
            rerun_totpe      = array.array('f',[0.0])
            rerun_maxpefrac  = array.array('f',[0.0])
            rerun_porchtotpe = array.array('f',[0.0])
            rerun_pass       = array.array('i',[0])
            outfvv_tree.SetBranchAddress( "TotPE", rerun_totpe )
            outfvv_tree.SetBranchAddress( "PorchTotPE", rerun_porchtotpe )
            outfvv_tree.SetBranchAddress( "MaxPEFrac", rerun_maxpefrac )
            outfvv_tree.SetBranchAddress( "PassPMTPrecut", rerun_pass )
        if self.rerun_1mu1p_bdt:
            rerun_1mu1p_cosmic = array.array('f',[0.0])
            rerun_1mu1p_nu     = array.array('f',[0.0])
            outfvv_tree.SetBranchAddress("BDTscore_1mu1p_cosmic",rerun_1mu1p_cosmic)
            outfvv_tree.SetBranchAddress("BDTscore_1mu1p_nu",rerun_1mu1p_nu)
        if self.rerun_1e1p_bdt:
            rerun_1e1p = array.array('f',[0.0])
            outfvv_tree.SetBranchAddress("BDTscore_1e1p",rerun_1e1p)


        for ientry in xrange( finalvertextree.GetEntries() ):
            finalvertextree.GetEntry(ientry)
            rse  = ( finalvertextree.run, finalvertextree.subrun, finalvertextree.event)
            rsev = ( finalvertextree.run, finalvertextree.subrun, finalvertextree.event, finalvertextree.vtxid)
            if rse in self.rse_dict and self.rse_dict[rse]:
                nverticesout+=1
                print "[Vertex passes]",rse," vtxid=",finalvertextree.vtxid
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
                if self.rerun_1mu1p_bdt:
                    print 'replacing 1mu1p results from rerun of rsev[',rsev,']'
                    rerun_1mu1p_cosmic[0] = 0.0
                    rerun_1mu1p_nu[0]     = self.bdtoutput_1mu1p[rsev]

                if self.rerun_1e1p_bdt:
                    print 'replacing 1e1p results from rerun of rsev[',rsev,']'
                    rerun_1e1p[0] = self.bdtoutput_1e1p[rsev]

                for tree in self.out_vertex_indexed_trees:
                    tree.Fill()

        if self.pot_sum_tree is not None:
            print "Copy POT summary tree"
            outpot = self.pot_sum_tree.CloneTree(0)
            abytes = self.pot_sum_tree.GetEntry(0)
            ii = 0
            while abytes>0:
                outpot.Fill()
                ii += 1
                abytes = self.pot_sum_tree.GetEntry(ii)

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

	    # get rid of pmtprecuts to rely entirely on common optical filter (basically only getting rid of maxfrac requirement)

            if self.rerun_1mu1p_bdt:
                print "replaced bdt scores with recalculated ones"
                print "  nu: old=",dlanatree.BDTscore_1mu1p_nu," new=",self.bdtoutput_1mu1p[rsev]
                bdtscore_1mu1p_cosmic = 0.0
                bdtscore_1mu1p_nu     = self.bdtoutput_1mu1p[rsev]
            else:
                bdtscore_1mu1p_cosmic = dlanatree.BDTscore_1mu1p_cosmic
                bdtscore_1mu1p_nu = dlanatree.BDTscore_1mu1p_nu

            if ( dlanatree.PassSimpleCuts==1
                 and dlanatree.MaxShrFrac<0.2
                 and dlanatree.OpenAng>0.5
                 and dlanatree.ChargeNearTrunk>0
                 and dlanatree.FailedBoost_1m1p!=1
                 and bdtscore_1mu1p_nu>0.5 ):
                # note this cut is not compatible with single bdt!!
                passes = True
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            # update RSE dictionary if
            #  no previous entry
            #  previous entry and this passes
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes
            # add to RSEV dictionary
            self.rsev_dict[rsev] = passes

            print "RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  maxshrfrac: ",dlanatree.MaxShrFrac<0.2," (",dlanatree.MaxShrFrac,")"
            print "  opening angle: ",dlanatree.OpenAng>0.5," (",dlanatree.OpenAng,")"
            print "  chargeneartrunk: ",dlanatree.ChargeNearTrunk>0," (",dlanatree.ChargeNearTrunk,")"
            print "  failedboost_1m1p: ",dlanatree.FailedBoost_1m1p!=1," (",dlanatree.FailedBoost_1m1p,")"
            print "  bdt score: ",bdtscore_1mu1p_nu>0.5," (",bdtscore_1mu1p_nu,")"

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
        print "[ dlfilter::run_1e1p_highE_filter ] first pass to get max 1e1p BDT score for passing vtx per event"
        max_rse = {}
        rse_vtxid = {}

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

            if self.rerun_1e1p_bdt:
                print "[highE filter] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",dlanatree.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                bdtscore_1e1p = self.bdtoutput_1e1p[rsev]
            else:
                bdtscore_1e1p = dlanatree.BDTscore_1e1p

            if ( dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2 ):
                passes = True
            
            print "[first pass highE] RSE=",rse," RSEV=",rsev," Passes=",passes
            #print "  precuts: ",passprecuts==1
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  maxshrfrac: ",max(dlanatree.MaxShrFrac,-1)>0.2," (",dlanatree.MaxShrFrac,")"
            print "  electron edep: ",dlanatree.Electron_Edep>35.0," (",dlanatree.Electron_Edep,")"
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"
            print "  enu: ",dlanatree.Enu_1e1p
            print "  bdt 1e1p>0.7: ",bdtscore_1e1p>0.7," (",bdtscore_1e1p,")"

            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":bdtscore_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                # update max bdt for vertex that passes
                print "  new max bdt for passing vtx candidates"
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}

                

        print "[ dlfilter::run_1e1p_highE_filter ] cutting pass"
        # next, save only those events, whose highest bdt score pass threshold
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "[highE] RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "[highE] RSE ",rse," vtxid=",dlanatree.vtxid,": ",max_rse[rse]," -- is not max BDT vertex: ",max_rse[rse]["vtxid"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["enu"]<700.0 or max_rse[rse]["bdt"]<0.7:
                print "[highE] RSE ",rse,": enu=",max_rse[rse]["enu"]," is BELOW energy threshold, below bdt=",max_rse[rse]["bdt"]," or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "[highE] RSE ",rse,": enu=",max_rse[rse]["enu"]," passes high-E selection (",max_rse[rse]["passes"],")"
            passes = True
            
            # update event flag
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes

            # update vertex flag
            self.rsev_dict[rsev] = passes

        return

    def run_1e1p_nearE_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "[ dlfilter::run_1e1p_nearE_filter ]"

        # we first collect the highest energy vertex
        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            passes = False
            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if ( passprecuts==1
                 and dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2 ):
                passes = True
                
            print "RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  precuts: ",passprecuts==1
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  maxshrfrac: ",max(dlanatree.MaxShrFrac,-1)>0.2," (",dlanatree.MaxShrFrac,")"
            print "  electron edep: ",dlanatree.Electron_Edep>35.0," (",dlanatree.Electron_Edep,")"
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"

            print "[first pass] RSE=",rse," RSEV=",rsev
            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":0.0,"passes":False}

            if passes and max_rse[rse]["enu"]<dlanatree.Enu_1e1p:
                # update if energy is larger
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":dlanatree.BDTscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                print "RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "RSE ",rse,": ",max_rse[rse]," -- is not max energy vertex: this=",dlanatree.Enu_1e1p," max=",max_rse[rse]["enu"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["enu"]<500.0:
                print "RSE ",rse,": energy max=",max_rse[rse]["enu"]," is below energy threshold or did not pass (",max_rse[rse]["passes"],")"
                continue
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            # PASSES, MARK EVENT AND VERTEX TO BE SAVED
            self.rse_dict[rse]   = True
            self.rsev_dict[rsev] = True


        return

    def run_1e1p_lowBDT_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "[ dlfilter::run_1e1p_lowBDT_filter ]"
        #we first collect the highest bdtscore per RSE

        max_rse = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            passes = False
            # no longer using precuts
            #passprecuts = int(dlanatree.PassPMTPrecut)
            #if self.rerun_pmtprecuts:
            #    passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
            #    print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
            #    print self.PMTPrecut_Dict[rse]['_passpmtprecut']
            #    passprecuts = passrerun

            if self.rerun_1e1p_bdt:
                print "replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",dlanatree.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                bdtscore_1e1p = self.bdtoutput_1e1p[rsev]
            else:
                bdtscore_1e1p = dlanatree.BDTscore_1e1p

            if ( dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2 ):
                passes = True
                
            print "RSE=",rse," RSEV=",rsev," Passes=",passes
            #print "  precuts: ",passprecuts==1
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  maxshrfrac: ",max(dlanatree.MaxShrFrac,-1)>0.2," (",dlanatree.MaxShrFrac,")"
            print "  electron edep: ",dlanatree.Electron_Edep>35.0," (",dlanatree.Electron_Edep,")"
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"
            print "  bdt 1e1p: ",bdtscore_1e1p<=0.7," (",bdtscore_1e1p,")"

            print "[first pass] RSE=",rse," RSEV=",rsev
            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":bdtscore_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                # update max bdt for vertex that passes
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "RSE ",rse," vtxid=",dlanatree.vtxid,": ",max_rse[rse]," -- is not max BDT vertex: ",max_rse[rse]["vtxid"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]>0.7:
                print "RSE ",rse,": score max=",max_rse[rse]["bdt"]," is above threshold or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "RSE ",rse,": score max=",max_rse[rse]["bdt"]," passes low-BDT selection (",max_rse[rse]["passes"],")"
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            self.rse_dict[rse]   = True
            self.rsev_dict[rsev] = True


        # for debug only
        #rsekeys = self.rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",self.rsev_dict[k]

        return

    def run_1e1p_midBDT_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "[ dlfilter::run_1e1p_lowBDT_filter ]"
        #we first collect the highest bdtscore per RSE

        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            passes = False
            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if ( passprecuts==1
                 and dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2 ):
                passes = True
                
            print "RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  precuts: ",passprecuts==1
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  maxshrfrac: ",max(dlanatree.MaxShrFrac,-1)>0.2," (",dlanatree.MaxShrFrac,")"
            print "  electron edep: ",dlanatree.Electron_Edep>35.0," (",dlanatree.Electron_Edep,")"
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"
            print "  bdt 1e1p: ",dlanatree.BDTscore_1e1p<=0.7," (",dlanatree.BDTscore_1e1p,")"

            print "[first pass] RSE=",rse," RSEV=",rsev
            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":dlanatree.Enu_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<dlanatree.BDTscore_1e1p:
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":dlanatree.BDTscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "RSE ",rse,": ",max_rse[rse]," -- is not max BDT vertex: this=",dlanatree.BDTscore_1e1p," max=",max_rse[rse]["bdt"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]>0.7 or max_rse[rse]["bdt"]<0.4:
                print "RSE ",rse,": score max=",max_rse[rse]["bdt"]," is above threshold or did not pass (",max_rse[rse]["passes"],")"
                continue
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            self.rse_dict[rse]   = True
            self.rsev_dict[rsev] = True


        # for debug only
        #rsekeys = self.rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",self.rsev_dict[k]

        return

    def run_lowBDT_pi0_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        
        pi0cuts = 'Proton_Edep>60.0  and Electron_Edep>35.0 and PassPMTPrecut==1 and PassShowerReco==1  
        and shower1_E_Y>80 and ChargeNearTrunk >250 and Electron_ThetaRecoB_e1ep <1.5 and _shower_alpha <2.5 and _pi0mass>0   '

        """
        print "[ dlfilter::run_lowBDT_pi0_filter ] first pass"
        #we first collect the highest bdtscore per RSE

        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            passes = False
            passprecuts = int(dlanatree.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "[pi0 filter] replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                print self.PMTPrecut_Dict[rse]['_passpmtprecut']
                passprecuts = passrerun

            if self.rerun_1e1p_bdt:
                print "[pi0 filter] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",dlanatree.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                bdtscore_1e1p = self.bdtoutput_1e1p[rsev]
            else:
                bdtscore_1e1p = dlanatree.BDTscore_1e1p

            if ( dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60 
                 and dlanatree.Electron_Edep > 35
                 and dlanatree.shower1_E_Y>80.0
                 and dlanatree.ChargeNearTrunk>250.0
                 and dlanatree.Electron_ThetaRecoB_e1ep<1.5
                 and dlanatree._shower_alpha<2.5
                 and dlanatree._pi0mass>0 ):
                passes = True
                
            print "[pi0 filter first pass] RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"
            print "  bdt 1e1p<0.7: ",bdtscore_1e1p<=0.7," (",bdtscore_1e1p,")"
            print "  shower1_E_Y: ",dlanatree.shower1_E_Y>80.0,"(",dlanatree.shower1_E_Y,")"
            print "  Qtrunk: ",dlanatree.ChargeNearTrunk>250.0,"(",dlanatree.ChargeNearTrunk,")"
            print "  Electron ThetaRecoB_1e1p: ",dlanatree.Electron_ThetaRecoB_e1ep<1.5,"(",dlanatree.Electron_ThetaRecoB_e1ep,")"
            print "  _shower_alpha: ",dlanatree._shower_alpha<2.5,"(",dlanatree._shower_alpha,")"
            print "  pi0mass: ",dlanatree._pi0mass>0,"(",dlanatree._pi0mass,")"

            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":dlanatree.Enu_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        print "[ dlfilter::run_lowBDT_pi0_filter ] selection pass"
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "[pi0 filter] RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "[pi0 filter]] RSE ",rse,": ",max_rse[rse]," -- is not max BDT vertex: this=",dlanatree.BDTscore_1e1p," max=",max_rse[rse]["bdt"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]>0.7:
                print "[pi0 filter] RSE ",rse,": score max=",max_rse[rse]["bdt"]," is above threshold or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "[pi0 filter] RSE ",rse,": score max=",max_rse[rse]["bdt"]," PASSES (",max_rse[rse]["passes"],")"
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            self.rse_dict[rse]   = True
            self.rsev_dict[rsev] = True


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
        print "[ dlfilter::run_1e1p_signal_filter ] first pass to get max 1e1p BDT score for passing vtx per event"
        max_rse = {}
        rse_vtxid = {}

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

            if self.rerun_1e1p_bdt:
                print "[signal filter] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",dlanatree.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                bdtscore_1e1p = self.bdtoutput_1e1p[rsev]
            else:
                bdtscore_1e1p = dlanatree.BDTscore_1e1p            


            if ( dlanatree.PassSimpleCuts==1
                 and dlanatree.PassShowerReco==1
                 and dlanatree.Proton_Edep > 60
                 and dlanatree.Electron_Edep > 35
                 and max(dlanatree.MaxShrFrac,-1) > 0.2):
                passes = True
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug

            print "[first pass signal] RSE=",rse," RSEV=",rsev," Passes=",passes
            #print "  precuts: ",passprecuts==1
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  maxshrfrac: ",max(dlanatree.MaxShrFrac,-1)>0.2," (",dlanatree.MaxShrFrac,")"
            print "  electron edep: ",dlanatree.Electron_Edep>35.0," (",dlanatree.Electron_Edep,")"
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"
            print "  enu: ",dlanatree.Enu_1e1p
            print "  bdt 1e1p>0.7: ",bdtscore_1e1p>0.7," (",bdtscore_1e1p,")"

            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":bdtscore_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                # update max bdt for vertex that passes
                print "  new max bdt for passing vtx candidates"
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}

        print "[ dlfilter::run_1e1p_signal_filter ] cutting pass"
        # next, save only those events, whose highest bdt score pass threshold
        self.rse_dict = {}
        self.rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "[signal] RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "[signal] RSE ",rse," vtxid=",dlanatree.vtxid,": ",max_rse[rse]," -- is not max BDT vertex: ",max_rse[rse]["vtxid"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]<0.7:
                print "[signal] RSE ",rse,": below bdt=",max_rse[rse]["bdt"]," or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "[signal] RSE ",rse,": enu=",max_rse[rse]["enu"]," bdt=",max_rse[rse]["bdt"]," passes signal selection (",max_rse[rse]["passes"],")"
            passes = True
            
            # update event flag
            if rse not in self.rse_dict:
                self.rse_dict[rse]   = passes
            elif rse in self.rse_dict and passes:
                self.rse_dict[rse]   = passes

            # update vertex flag
            self.rsev_dict[rsev] = passes

        return

    def run_rse_filter(self, dlanatree ):
        """ use a (run,subrun,event) list to filter """
        print "[ dlfilter::RSE filter ]"
        self.rse_dict = {}
        self.rsev_dict = {}

        flist = open(self.filter_pars["rse-list"],'r')
        llist = flist.readlines()
        rse_list = []
        for l in llist:
            info = l.strip().split()
            rse = (int(info[0]),int(info[1]),int(info[2]))
            rse_list.append(rse)
        rse_list.sort()
        print "[ flfilter::RSE filter ] number of (rse) in list: ",len(rse_list)

        for ientry in xrange(dlanatree.GetEntries()):
            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            if rse in rse_list:
                self.rse_dict[rse] = True
                self.rsev_dict[rsev] = True
            elif rse in self.rse_dict:
                pass
            else:
                self.rse_dict[rse] = False


        return
