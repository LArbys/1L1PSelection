#! /usr/bin/env python
###############################################################################
#
# Name: dl_multi_filter.py
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
from ROOT import TFile, std

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
import mpidutil
import varutils
from lib_mpid_torch.rootdata_pid import ROOTData
import bdt1e1p_helper
import cutdefinitions
#from SelectionDefs import apply_precuts

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

    obj = DLMultiFilter(config)
    return obj

# Analyze hit class

class DLMultiFilter(RootAnalyze):

    def __init__(self, config):
        #----------------------------------------------------------------------
        #
        # Purpose: Constructor.
        #
        #----------------------------------------------------------------------

        self.input_file_list = [] # gets append to during open_input function
        self.filter_pars = config['modules']['dl_multi_filter']
        self.filter_types = self.filter_pars['filter_type']
        self.DATARUN = self.filter_pars['data_run']
        if "remove_duplicate_vertices" in self.filter_pars:
            self.REMOVE_DUP_VTX = self.filter_pars["remove_duplicate_vertices"]
        else:
            self.REMOVE_DUP_VTX = False

        if "maxentries" in self.filter_pars:
            self.MAXENTRIES = int(self.filter_pars["maxentries"])
        else:
            self.MAXENTRIES = None

        if type(self.filter_types) is str:
            self.filter_types = [self.filter_types]
            
        if type(self.filter_types) is not list:
            raise TypeError("'filter_type' parameter should be a str with single filter type or list of filter types to run")

        
        self.sample_type = self.filter_pars['sample_type']
        if 'event_tree' in self.filter_pars:
            self.event_tree  = self.filter_pars['event_tree']
        else:
            self.event_tree  = "dlana/ubdlana_id_tree"

        DEFINED_FILTERS = ["numu-sideband",
                           "pi0-lowBDT-sideband",
                           "1e1p-highE-sideband",
                           "1e1p-nearE-sideband",
                           "1e1p-far-sideband",
                           "1e1p-midBDT-sideband",
                           "1e1p-signal",
                           "1e1p-loose-signal",
                           "rse-list"]
        
        for ftype in self.filter_types:
            if ftype not in DEFINED_FILTERS:
                raise ValueError("Invalid filter type specified [{}]. Defined: {}".format(ftype,DEFINED_FILTERS))

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
            elif self.filter_pars["bdt_model_1mu1p"]=="ensemble":
                self.bdt_model_1mu1p = bdtutil.load_BDT_ensemble( "1m1p", self.bdt_1mu1p_weightfile, nbdts=10, runs=[self.DATARUN] )
            print "DLMultiFilter: RERUN 1MU1P BDT"
        else:
            self.rerun_1mu1p_bdt = False

        if "rerun_1e1p_bdt" in self.filter_pars and bool(self.filter_pars["rerun_1e1p_bdt"]):
            self.rerun_1e1p_bdt = True            
            self.bdt_1e1p_weightfile = self.filter_pars["bdt_weights_1e1p"]
            if self.filter_pars["bdt_model_1e1p"]=="single":
                self.bdt_model_1e1p = bdtutil.load_BDT_model( self.bdt_1e1p_weightfile )
            elif self.filter_pars["bdt_model_1e1p"]=="ensemble":                
                self.bdt_model_1e1p = bdtutil.load_BDT_ensemble( "1e1p", self.bdt_1e1p_weightfile, nbdts=20, runs=[self.DATARUN] )
            print "DLMultiFilter: RERUN 1e1P BDT"
        else:
            self.rerun_1e1p_bdt = False            

        if "rerun_mpid" in self.filter_pars and self.filter_pars["rerun_mpid"]==True:
            self.mpid_cfg_path = os.environ["UBMPIDNET_DIR"]+"/production_cfg/inference_config_tufts_WC.cfg"
            self.RUN_MPID = True
        else:
            self.mpid_cfg_path = False
            self.RUN_MPID = False

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
        print "DLMultiFilter::event_info = ",rse
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
        # WE DO A LOT IN THIS FUNCTION, BREAKING THE DESIGN
        # Design is event-centric
        # unfortunately, our trees are vertex candidate-centric
        #
        # We recalculate per vertex-quantities first
        # Then we determine which events and vertices pass
        #
        # Then we actually write the output files.
        #
        #----------------------------------------------------------------------

        print 'DLMultiFilter::open_input called.'
        print 'Input file: %s' % input_file.GetName()
        #input_file.ls()
        self.input_file_list.append(input_file.GetName())
        
        # determine the entry-indexed andd vertex-indexed trees we are using
        # as reference for both types

        # EVENT-INDEX TREE
        larlite_id_tree = input_file.Get(self.event_tree)

        # VERTEX-INDEX TREE
        finalvertextree = input_file.Get("dlana/FinalVertexVariables")
        self._use_ubdlana_idtree = False
        if self.event_tree=="dlana/ubdlana_id_tree":
            try:
                print "DLMultiFilter::Using dlana/ubdlana_id_tree"
                nevent_entries  = larlite_id_tree.GetEntries()
                self._use_ubdlana_idtree = True                
            except:
                print "DLMultiFilter::Could not load ubdlana. try larlite_id_tree"
                larlite_id_tree = None
                nevent_entries = 0
                self._use_ubdlana_idtree = False
        if self.event_tree=="larlite_id_tree" or larlite_id_tree is None:
            larlite_id_tree = input_file.Get("larlite_id_tree")
            try:
                nevent_entries = larlite_id_tree.GetEntries()
                print "DLMultiFilter: using larlite_id_tree"
                self._use_ubdlana_idtree = False
            except:
                print "DLMultiFilter::Could not load ubdlana_id_tree either"
                larlite_id_tree = None
                nevent_entries = 0
        print "DLMultiFilter::event id tree: ",self.event_tree
        print "DLMultiFilter::_use_ubdlana_idtree: ",self._use_ubdlana_idtree

        nevent_entries  = larlite_id_tree.GetEntries()
        nvertex_entries = finalvertextree.GetEntries()

        # SORT THE TYPE OF TREES IN THE FILE
        self.vertex_indexed_trees = []
        self.event_indexed_trees  = []
        self.other_trees = []
        self.pot_sum_tree = None
        self.fvv_tree = None        
        self.vertex_indexed_trees,self.event_indexed_trees,self.other_trees,self.fvv_tree,self.pot_sum_tree,dirdict = self.get_tree_lists(input_file,nvertex_entries,nevent_entries)

        # FOR DEBUG
        if False:
            print "[EVENT-INDEXED TREES]"
            for atree in self.event_indexed_trees:
                print atree.GetName()
        if True:
            print "[VERTEX-INDEXED TREES]"
            for atree in self.vertex_indexed_trees:
                print atree.GetName()


        # RERUN SOME VERTEX-VARIABLES
        
        # rerun precuts (Event-Tree)
        if self.rerun_pmtprecuts and nevent_entries>0:
            self.PMTPrecut_Dict = performPMTPrecuts( input_file.GetName(), **self.precutpars )

        # rerun the 1mu1p BDT (Vertex-Tree)
        # we get back a dictionary indexed by (run,subrun,event,vertex)
        # which stores replacement values for the BDT variables and BDT score
        if self.rerun_1mu1p_bdt:
            if self.filter_pars["bdt_model_1mu1p"]=="single":            
                self.bdtoutput_1mu1p = bdtutil.rerun_1mu1p_models( self.bdt_model_1mu1p, finalvertextree )
            else:
                self.bdtoutput_1mu1p = bdtutil.rerun_1mu1p_ensemble( self.bdt_model_1mu1p, finalvertextree, self.DATARUN, nbdts=10, maxentries=self.MAXENTRIES )                
        else:
            self.bdtoutput_1mu1p = {}

        # rerun the 1e1p BDT (Vertex-Tree)
        # we get back a dictionary indexed by (run,subrun,event,vertex)
        # which stores replacement values for the BDT variables and BDT score
        if self.rerun_1e1p_bdt:
            if self.filter_pars["bdt_model_1e1p"]=="single":                        
                self.bdtoutput_1e1p = bdtutil.rerun_1e1p_models( self.bdt_model_1e1p, finalvertextree )
            else:
                self.bdtoutput_1e1p, self.bdtrerun_electron_edep = \
                bdtutil.rerun_1e1p_ensemble( self.bdt_model_1e1p, finalvertextree, self.DATARUN, nbdts=10, maxentries=self.MAXENTRIES )
        else:
            self.bdtoutput_1e1p = {}
            self.bdtrerun_electron_edep = {}

        # setup the MPID CNN (Event-Tree)
        if self.RUN_MPID:
            # we'll need (1) to load the MPID model and (2) load the larcv data IOManager to read in wire images
            self.mpid, self.mpid_cfg = mpidutil.load_mpid_model( self.mpid_cfg_path )
            # larcv data
            self.iolcv = larcv.IOManager(larcv.IOManager.kREAD,"larcv")
            self.iolcv.specify_data_read( larcv.kProductImage2D, "wire" )
            self.iolcv.specify_data_read(larcv.kProductPGraph,"inter_par")
            self.iolcv.specify_data_read(larcv.kProductPixel2D,"inter_par_pixel")
            self.iolcv.specify_data_read(larcv.kProductPixel2D,"inter_img_pixel")
            self.iolcv.specify_data_read(larcv.kProductPixel2D,"inter_int_pixel")
            self.iolcv.add_in_file( input_file.GetName() )
            self.iolcv.initialize()

            # MAKE THE MPID DATACLASS AND TREE
            self.mpid_dir = self.output_file.mkdir("mpid_all")
            self.mpid_dir.cd()
            self.mpid_data = ROOTData()
            self.mpid_tree  = ROOT.TTree("multipid_tree","UB MPID score tree")
            self.mpid_data.init_tree(self.mpid_tree)
            self.mpid_data.reset()

            tot_mpid_vertices = 0
            self.mpid_results = {}
            for ientry in xrange( self.iolcv.get_n_entries() ):
                # Get the LARCV ENTRY
                self.iolcv.read_entry(ientry)
                
                # run mpid on entry
                nmpid_vertices, mpid_entry_results = mpidutil.run_mpid_on_larcv_entry( self.mpid_cfg, self.mpid, self.iolcv, self.mpid_data, self.mpid_tree,return_result_dict=True)
                for entry_data_key in mpid_entry_results:
                    self.mpid_results[entry_data_key] = mpid_entry_results[entry_data_key]
                tot_mpid_vertices += nmpid_vertices
            
            # Get MPID results for cuts
            print "==================================="
            print "MPID RESULT DUMP"
            for rsev in self.mpid_results:
                print rsev
                print self.mpid_results[rsev]
            print "==================================="
            
        else:
            self.mpid = None
            self.mpid_cfg = None
            self.iolcv = None
            self.mpid_results = None
            self.mpid_data = None
            self.mpid_tree = None
        

        # RUN THE FILTERS ON THE VERTICES, GET RESULTS, which is dictionary of RS and RSE that pass
        self.filters_results = {}
        self.filters_results = self.run_all_filters( finalvertextree, self.filter_types )

        # NOW WE MAKE THE OUTPUT FILES, ONE FOR EACH FILE
        # WE COULD ALSO TRY TO SAVE THE TREES IN ROOT FILE FOLDERS
        self.filter_trees = {}
        for filtertype in self.filter_types:
            # MAKE FILTER DIRECTORY
            filter_dir_name = filtertype.replace("-","_").replace("1e1p","nu1e1p")
            self.output_file.cd()
            try:
                filterdir = self.output_file.mkdir( filter_dir_name )
            except:
                filterdir = self.output_file.Get( filter_dir_name )
            filterdir.cd()

            # NOW WE MAKE OUTPUT TREES AND CLONE THEM FOR EACH
            print "====================================================="
            print "Cloning output trees for filter[",filtertype,"]"

            out_event_indexed_trees  = []
            out_vertex_indexed_trees = []

            outfvv_tree = None

            # clone event-indexed trees
            for tree in self.event_indexed_trees:
                if str(tree.GetName()) in dirdict:
                    try:
                        rootdir = self.output_file.mkdir( filter_dir_name+"/"+dirdict[str(tree.GetName())] )
                    except:
                        rootdir = self.output_file.Get( filter_dir_name+"/"+dirdict[str(tree.GetName())] )
                    rootdir.cd()
                out_event_indexed_trees.append( tree.CloneTree(0) )
                filterdir.cd()
                
            # clone vertex-indexed trees
            for tree in self.vertex_indexed_trees:
                if str(tree.GetName()) in dirdict:
                    rootdir = self.output_file.Get( filter_dir_name+"/"+dirdict[str(tree.GetName())] )
                    try:
                        rootdir.cd()
                    except:
                        rootdir = self.output_file.mkdir( filter_dir_name+"/"+dirdict[str(tree.GetName())] )
                        rootdir.cd()

                out_vertex_indexed_trees.append( tree.CloneTree(0) )
                if tree==self.fvv_tree:
                    outfvv_tree = out_vertex_indexed_trees[-1]
                filterdir.cd()

            # get filter dict results
            rse_dict  = self.filters_results[filtertype]["rse"]
            rsev_dict = self.filters_results[filtertype]["rsev"]
            
            neventsout = 0
            nverticesout = 0
            # loop over input event entries
            for ientry in xrange( nevent_entries ):
                
                if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                    break

                larlite_id_tree.GetEntry(ientry)
                if not self._use_ubdlana_idtree:
                    rse = ( larlite_id_tree._run_id, larlite_id_tree._subrun_id, larlite_id_tree._event_id )
                else:
                    rse = (larlite_id_tree.run,larlite_id_tree.subrun,larlite_id_tree.event)
                    
                if rse in rse_dict and rse_dict[rse]:
                    neventsout+=1
                    for tree in self.event_indexed_trees:
                        tree.GetEntry(ientry)
                    for tree in out_event_indexed_trees:
                        tree.Fill()
                else:
                    if rse in rse_dict:
                        print "[",filtertype,": do not save event] ",rse," result=",rse_dict[rse]
                    else:
                        print "[",filtertype,": do not save event] ",rse," not in RSE dict"

            # VERTEX-INDEXED TREES OUTPUT
            # if rerunning steps, we have to replace the branch addresses of the FVV with new ones
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
                bdtvars_1e1p = {}
                varnames = bdt1e1p_helper.get_bdt_var_names()
                for varname in varnames:
                    if varname is None:
                        continue
                    bdtvars_1e1p[varname] = array.array('f',[0.0])
                    outfvv_tree.SetBranchAddress(varname, bdtvars_1e1p[varname] )
            if self.RUN_MPID:
                """MuonPID_int_v:ProtonPID_int_v:EminusPID_int_v"""
                varnames = mpidutil.get_fvv_var_names()
                mpidvars_1e1p = {}
                for varname in varnames:
                    mpidvars_1e1p[varname] = std.vector("float")(3,0)
                    outfvv_tree.SetBranchAddress(varname, mpidvars_1e1p[varname])
                

            # FINAL VERTEX TREE LOOP
            for ientry in xrange( finalvertextree.GetEntries() ):
                if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                    break
                finalvertextree.GetEntry(ientry)
                rse  = ( finalvertextree.run, finalvertextree.subrun, finalvertextree.event)
                rsev = ( finalvertextree.run, finalvertextree.subrun, finalvertextree.event, finalvertextree.vtxid)

                if rse in rse_dict and rse_dict[rse]:
                    
                    if self.REMOVE_DUP_VTX  and (rsev not in rsev_dict or not rsev_dict[rsev]):
                        print "[",filtertype,": Vertex fails]",rse," vtxid=",finalvertextree.vtxid
                        continue

                    nverticesout+=1                    
                    print "[",filtertype,": Vertex passes]",rse," vtxid=",finalvertextree.vtxid
                    for tree in self.vertex_indexed_trees:
                        tree.GetEntry(ientry)

                    if self.rerun_pmtprecuts:
                        print '[',filtertype,': replacing precut results with those from rerun of rse[',rse,']: '
                        print ' totpe old=',self.fvv_tree.TotPE,' new=',self.PMTPrecut_Dict[rse]['_totpe']
                        print ' maxpefrac old=',self.fvv_tree.MaxPEFrac,' new=',self.PMTPrecut_Dict[rse]['_maxpefrac']
                        rerun_totpe[0]      = self.PMTPrecut_Dict[rse]['_totpe']
                        rerun_porchtotpe[0] = self.PMTPrecut_Dict[rse]['_porchtotpe']
                        rerun_maxpefrac[0]  = self.PMTPrecut_Dict[rse]['_maxpefrac']
                        rerun_pass[0]       = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                    if self.rerun_1mu1p_bdt:
                        print '[',filtertype,': replacing 1mu1p results from rerun of rsev[',rsev,']'
                        rerun_1mu1p_cosmic[0] = 0.0
                        rerun_1mu1p_nu[0]     = self.bdtoutput_1mu1p[rsev]

                    if self.rerun_1e1p_bdt:
                        print '[',filtertype,': replacing 1e1p results from rerun of rsev[',rsev,']'
                        rerun_1e1p[0] = self.bdtoutput_1e1p[rsev]
                        varnames = bdt1e1p_helper.get_bdt_var_names()
                        input_vars = bdt1e1p_helper.getNewShowerCalibTrainingVarbs( self.fvv_tree, newCalib=True )
                        for n,varname in enumerate(varnames):
                            if varname in bdtvars_1e1p:
                                bdtvars_1e1p[varname][0] = input_vars[n]

                    if self.RUN_MPID:
                        print '[',filtertype,': replacing mpid results of rsev[',rsev,']'
                        mpid_vtx_results = self.mpid_results[rsev]
                        for iplane in range(3):
                            mpidvars_1e1p["EminusPID_int_v"][iplane]   = mpid_vtx_results[(iplane,"int")][0]
                            mpidvars_1e1p["GammaPID_int_v"][iplane]    = mpid_vtx_results[(iplane,"int")][1]
                            mpidvars_1e1p["MuonPID_int_v"][iplane]     = mpid_vtx_results[(iplane,"int")][2]
                            mpidvars_1e1p["ProtonPID_int_v"][iplane]   = mpid_vtx_results[(iplane,"int")][4]
                            mpidvars_1e1p["EminusPID_pix_v"][iplane]   = mpid_vtx_results[(iplane,"pix")][0]
                            mpidvars_1e1p["GammaPID_pix_v"][iplane]    = mpid_vtx_results[(iplane,"pix")][1]
                            mpidvars_1e1p["MuonPID_pix_v"][iplane]     = mpid_vtx_results[(iplane,"pix")][2]
                            mpidvars_1e1p["ProtonPID_pix_v"][iplane]   = mpid_vtx_results[(iplane,"pix")][4]


                    for tree in out_vertex_indexed_trees:
                        tree.Fill()



            print "=====[SUMMARY: ",filtertype,"]=================================="
            print "Num of output event-indexed trees: ",len(out_event_indexed_trees)
            print "Num of output vertex-indexed trees: ",len(out_vertex_indexed_trees)
            print "Num of events saved: ",neventsout
            print "Num of vertices saved: ",nverticesout
            self.filter_trees[filtertype] = { "event-indexed":out_event_indexed_trees,"vertex-indexed":out_vertex_indexed_trees}
                        
        # THE POT SUMMARY TREE
        print "=== POT SUMMARY TREE ===="
        self.output_file.cd()
        if self.pot_sum_tree is not None:
            print "Copy POT summary tree"
            outpot = self.pot_sum_tree.CloneTree(0)
            abytes = self.pot_sum_tree.GetEntry(0)
            ii = 0
            while abytes>0:
                outpot.Fill()
                ii += 1
                abytes = self.pot_sum_tree.GetEntry(ii)
        else:
            print "No POT Summary Tree to copy"
                

        print "[ End input tree prep ]"
        print "================================"



    def end_job(self):
        """ end of job tasks """
        print "[End of Job]"
        self.output_file.ls()


    def run_numu_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "==== [DL FILTER: 1mu1p ] ==============="
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

	    # get rid of pmtprecuts to rely entirely on common optical filter (basically only getting rid of maxfrac requirement)

            if self.rerun_1mu1p_bdt:
                print "[1mu1p] replaced bdt scores with recalculated ones"
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
                 and bdtscore_1mu1p_nu>=0.5 ):
                passes = True
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            # update RSE dictionary if
            #  no previous entry
            #  previous entry and this passes
            if rse not in rse_dict:
                rse_dict[rse]   = passes
            elif rse in rse_dict and passes:
                rse_dict[rse]   = passes
            # add to RSEV dictionary
            rsev_dict[rsev] = passes

            print "[1m1p] RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  maxshrfrac: ",dlanatree.MaxShrFrac<0.2," (",dlanatree.MaxShrFrac,")"
            print "  opening angle: ",dlanatree.OpenAng>0.5," (",dlanatree.OpenAng,")"
            print "  chargeneartrunk: ",dlanatree.ChargeNearTrunk>0," (",dlanatree.ChargeNearTrunk,")"
            print "  failedboost_1m1p: ",dlanatree.FailedBoost_1m1p!=1," (",dlanatree.FailedBoost_1m1p,")"
            print "  bdt score: ",bdtscore_1mu1p_nu>=0.5," (",bdtscore_1mu1p_nu,")"

        # for debug only
        #rsekeys = rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",rsev_dict[k]

        return rse_dict, rsev_dict

    def run_1e1p_highE_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "////////////////////////////////////////////////////////////"
        print "[ dl_multi_filter::run_1e1p_highE_filter ] first pass to get max 1e1p BDT score for passing vtx per event"
        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            # Collect selection variables from FVV tree
            x = varutils.make_selection_vars_from_fvv( dlanatree )

            # UPDATE VALUES IF WE'VE RERUN SOME PORTIONS
            passprecuts = int(x.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                varutils.update_pmt_precuts(x,rse,self.PMTPrecut_Dict)
                passprecuts = int(x.PassPMTPrecut)

            if self.rerun_1e1p_bdt:
                print "[1e1p highE] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",x.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                varutils.update_bdt1e1p(x,rsev,self.bdtoutput_1e1p)
            bdtscore_1e1p = x.BDTscore_1e1p

            if self.RUN_MPID:
                print "[1e1p highE] replacing MPID scores with recalculated ones"
                if rsev in self.mpid_results:
                    varutils.update_mpid(x,rsev,self.mpid_results)

            # apply 1e1p cuts
            passes = cutdefinitions.precuts(x,self.DATARUN) and cutdefinitions.postcuts(x,self.DATARUN)

            # apply highE filter (subset of far-sideband)
            if x.Enu_1e1p<700.0 or x.BDTscore_1e1p<0.95:
                passes = False

            print "[1e1p highE first pass] RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  precuts: ",x.PassPMTPrecut==1
            print "  simplecuts: ",x.PassSimpleCuts==1
            print "  showerreco: ",x.PassShowerReco==1
            print "  maxshrfrac: ",max(x.MaxShrFrac,-1)>0.2," (",x.MaxShrFrac,")"
            print "  electron edep: ",x.Electron_Edep>35.0," (",x.Electron_Edep,")"
            print "  proton edep: ",x.Proton_Edep>60.0," (",x.Proton_Edep,")"
            print "  enu: ",x.Enu_1e1p>700.0," (",x.Enu_1e1p,")"
            print "  bdt 1e1p>0.95: ",bdtscore_1e1p>0.95," (",bdtscore_1e1p,")"

            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":bdtscore_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                # update max bdt for vertex that passes
                print "  new max bdt for passing vtx candidates"
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":x.Enu_1e1p,"passes":True}
                

        print "[ dl_multi_filter::run_1e1p_highE_filter ] cutting pass"
        # next, save only those events, whose highest bdt score pass threshold
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "[1e1p highE] RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "[1e1p highE] RSE ",rse," vtxid=",dlanatree.vtxid,": ",max_rse[rse]," -- is not max BDT vertex: ",max_rse[rse]["vtxid"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["enu"]<700.0 or max_rse[rse]["bdt"]<0.7:
                print "[1e1p highE] RSE ",rse,": enu=",max_rse[rse]["enu"]," is BELOW energy threshold, below bdt=",max_rse[rse]["bdt"]," or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "[1e1p highE] RSE ",rse,": enu=",max_rse[rse]["enu"]," passes high-E selection (",max_rse[rse]["passes"],")"
            passes = True
            
            # update event flag
            if rse not in rse_dict:
                rse_dict[rse]   = passes
            elif rse in rse_dict and passes:
                rse_dict[rse]   = passes

            # update vertex flag
            rsev_dict[rsev] = passes

        return rse_dict, rsev_dict

    def run_1e1p_nearE_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "/////////////////////////////////////////////////////"
        print "[ dl_multi_filter::run_1e1p_nearE_filter ]"

        # we first collect the highest energy vertex
        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            # Collect selection variables from FVV tree
            x = varutils.make_selection_vars_from_fvv( dlanatree )

            # UPDATE VALUES IF WE'VE RERUN SOME PORTIONS
            passprecuts = int(x.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                varutils.update_pmt_precuts(x,rse,self.PMTPrecut_Dict)
                passprecuts = int(x.PassPMTPrecut)

            if self.rerun_1e1p_bdt:
                print "[1e1p nearE] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",x.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                varutils.update_bdt1e1p(x,rsev,self.bdtoutput_1e1p)
            bdtscore_1e1p = x.BDTscore_1e1p

            if self.RUN_MPID:
                print "[1e1p nearE] replacing MPID scores with recalculated ones"
                if rsev in self.mpid_results:
                    varutils.update_mpid(x,rsev,self.mpid_results)

            # apply 1e1p cuts
            passes = cutdefinitions.precuts(x,self.DATARUN) and cutdefinitions.postcuts(x,self.DATARUN)

            # apply nearE filter
            if x.Enu_1e1p<500.0 or x.Enu_1e1p>700.0 or x.BDTscore_1e1p<0.95:
                passes = False

            print "[1e1p nearE first pass] RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  precuts: ",passprecuts==1
            print "  simplecuts: ",dlanatree.PassSimpleCuts==1
            print "  showerreco: ",dlanatree.PassShowerReco==1
            print "  maxshrfrac: ",max(dlanatree.MaxShrFrac,-1)>0.2," (",dlanatree.MaxShrFrac,")"
            print "  electron edep: ",dlanatree.Electron_Edep>35.0," (",dlanatree.Electron_Edep,")"
            print "  proton edep: ",dlanatree.Proton_Edep>60.0," (",dlanatree.Proton_Edep,")"
            print "  enu: ",(x.Enu_1e1p>500.0 and x.Enud_1e1p<700.0)," (",x.Enu_1e1p,")"
            print "  bdt 1e1p>0.95: ",bdtscore_1e1p>0.95," (",bdtscore_1e1p,")"


            print "[first pass] RSE=",rse," RSEV=",rsev
            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":0.0,"passes":False}

            if passes and max_rse[rse]["enu"]<x.Enu_1e1p:
                # update if energy is larger
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":x.BDTscore_1e1p,"enu":x.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

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
            rse_dict[rse]   = True
            rsev_dict[rsev] = True


        return rse_dict, rsev_dict

    def run_1e1p_farsideband_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "///////////////////////////////////////////////////////////"
        print "[ dl_multi_filter::run_1e1p_farsideband_filter ]"
        #we first collect the highest bdtscore per RSE

        max_rse = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            passes = False
            # Collect selection variables from FVV tree
            x = varutils.make_selection_vars_from_fvv( dlanatree )

            # UPDATE VALUES IF WE'VE RERUN SOME PORTIONS
            passprecuts = int(x.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                varutils.update_pmt_precuts(x,rse,self.PMTPrecut_Dict)
                passprecuts = int(x.PassPMTPrecut)

            if self.rerun_1e1p_bdt:
                print "[1e1p far-sideband] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",x.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                varutils.update_bdt1e1p(x,rsev,self.bdtoutput_1e1p)
            bdtscore_1e1p = x.BDTscore_1e1p

            if self.RUN_MPID:
                print "[1e1p far-sideband] replacing MPID scores with recalculated ones"
                if rsev in self.mpid_results:
                    varutils.update_mpid(x,rsev,self.mpid_results)

            # apply 1e1p cuts
            passes = cutdefinitions.precuts(x,self.DATARUN) and cutdefinitions.postcuts(x,self.DATARUN)

            # signal+nearbox blocker
            if x.BDTscore_1e1p>0.95 and x.Enu_1e1p<700.0:
                passes = False
                
            print "[1e1p far-sideband] RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  precuts: ",passprecuts==1
            print "  simplecuts: ",x.PassSimpleCuts==1
            print "  showerreco: ",x.PassShowerReco==1
            print "  maxshrfrac: ",max(x.MaxShrFrac,-1)>0.2," (",x.MaxShrFrac,")"
            print "  electron edep: ",x.Electron_Edep>35.0," (",x.Electron_Edep,")"
            print "  proton edep: ",x.Proton_Edep>60.0," (",x.Proton_Edep,")"
            print "  bdt 1e1p: ",bdtscore_1e1p<0.95," (",bdtscore_1e1p,")"
            print "  bdt 1e1p AND low-E: ",x.BDTscore_1e1p>0.95 and x.Enu_1e1p<500.0

            print "[1e1p far-sideband first pass] RSE=",rse," RSEV=",rsev
            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":bdtscore_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                # update max bdt for vertex that passes
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            if rse not in max_rse:
                # surprising
                print "[1e1p far-sideband] RSE ",rse," not in max_rse dict"
                continue
            if max_rse[rse]["vtxid"]!=dlanatree.vtxid:
                # ignore non-max  BDT vertex
                print "[1e1p far-sideband] RSE ",rse," vtxid=",dlanatree.vtxid,": ",max_rse[rse]," -- is not max BDT vertex: ",max_rse[rse]["vtxid"]
                continue

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]>0.7:
                print "[1e1p far-sideband] RSE ",rse,": score max=",max_rse[rse]["bdt"]," is above threshold or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "[1e1p far-sideband] RSE ",rse,": score max=",max_rse[rse]["bdt"]," passes low-BDT selection (",max_rse[rse]["passes"],")"
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            rse_dict[rse]   = True
            rsev_dict[rsev] = True


        # for debug only
        #rsekeys = rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",rsev_dict[k]

        return rse_dict, rsev_dict

    def run_1e1p_midBDT_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "[ dl_multi_filter::run_1e1p_midBDT_filter ]"
        #we first collect the highest bdtscore per RSE

        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            
            # Collect selection variables from FVV tree
            x = varutils.make_selection_vars_from_fvv( dlanatree )

            # UPDATE VALUES IF WE'VE RERUN SOME PORTIONS
            passprecuts = int(x.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                varutils.update_pmt_precuts(x,rse,self.PMTPrecut_Dict)
                passprecuts = int(x.PassPMTPrecut)

            if self.rerun_1e1p_bdt:
                print "[1e1p mid-BDT] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",x.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                varutils.update_bdt1e1p(x,rsev,self.bdtoutput_1e1p)
            bdtscore_1e1p = x.BDTscore_1e1p

            if self.RUN_MPID:
                print "[1e1p mid-BDT] replacing MPID scores with recalculated ones"
                if rsev in self.mpid_results:
                    varutils.update_mpid(x,rsev,self.mpid_results)

            # apply 1e1p cuts
            passes = cutdefinitions.precuts(x,self.DATARUN) and cutdefinitions.postcuts(x,self.DATARUN)

            # low BDT cut off
            if x.BDTscore_1e1p>0.95 or x.BDTscore<0.4:
                passes = False

            # signal blocker
            if x.BDTscore_1e1p>0.95 and x.Enu_1e1p<500.0:
                passes = False

            print "[first pass] RSE=",rse," RSEV=",rsev
            if rse not in max_rse:
                # provide default
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":dlanatree.Enu_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<dlanatree.BDTscore_1e1p:
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":dlanatree.BDTscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}



        # next, save only those events, whose highest bdt score pass threshold
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

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

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]>0.95 or max_rse[rse]["bdt"]<0.4:
                print "RSE ",rse,": score max=",max_rse[rse]["bdt"]," is above threshold or did not pass (",max_rse[rse]["passes"],")"
                continue
            
            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug
                
            rse_dict[rse]   = True
            rsev_dict[rsev] = True


        # for debug only
        #rsekeys = rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",rsev_dict[k]

        return rse_dict, rsev_dict

    def run_lowBDT_pi0_filter(self, dlanatree ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        
        pi0cuts = 'Proton_Edep>60.0  and Electron_Edep>35.0 and PassPMTPrecut==1 and PassShowerReco==1  
        and shower1_E_Y>80 and ChargeNearTrunk >250 and Electron_ThetaRecoB_e1ep <1.5 and _shower_alpha <2.5 and _pi0mass>0   '

        """
        print "[ dl_multi_filter::run_lowBDT_pi0_filter ] first pass"
        #we first collect the highest bdtscore per RSE

        max_rse = {}
        rse_vtxid = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

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
                if rsev in self.bdtoutput_1e1p:
                    print "  nu: old=",dlanatree.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                    bdtscore_1e1p = self.bdtoutput_1e1p[rsev]
                else:
                    bdtscore_1e1p = 0.0
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
        print "[ dl_multi_filter::run_lowBDT_pi0_filter ] selection pass"
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

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
                
            rse_dict[rse]   = True
            rsev_dict[rsev] = True


        # for debug only
        #rsekeys = rsev_dict.keys()
        #rsekeys.sort()
        #for k in rsekeys:
        #    print k,": ",rsev_dict[k]

        return rse_dict, rsev_dict

    def run_1e1p_signal_filter(self, dlanatree, BDTCUT=0.95 ):
        """ use the final vertex tree to make selection 
        we create an RSE and RSEV dict
        """
        print "////////////////////////////////////////////////////////////"
        print "[ dl_multi_filter::run_1e1p_signal_filter ] first pass to get max 1e1p BDT score for passing vtx per event"
        max_rse = {}
        rse_vtxid = {}

        # FIRST LOOP COMPILES WHAT PASSES, FINDS MAX BDT-SCORE PER EVENT
        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            # Collect selection variables from FVV tree
            x = varutils.make_selection_vars_from_fvv( dlanatree )

            passes = True
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)

            # PRECUTS
            passprecuts = int(x.PassPMTPrecut)
            if self.rerun_pmtprecuts:
                passrerun = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
                print "replaced precut evaluation with rerun result. old=",passprecuts," new=",passrerun,
                varutils.update_pmt_precuts(x,rse,self.PMTPrecut_Dict)
                passprecuts = int(x.PassPMTPrecut)

            if self.rerun_1e1p_bdt:
                print "[signal filter] replacing 1e1p bdt scores with recalculated ones"
                print "  nu: old=",x.BDTscore_1e1p," new=",self.bdtoutput_1e1p[rsev]
                varutils.update_bdt1e1p(x,rsev,self.bdtoutput_1e1p)
            bdtscore_1e1p = x.BDTscore_1e1p

            if self.RUN_MPID:
                print "[signal filter] replacing MPID scores with recalculated ones"
                if rsev in self.mpid_results:
                    varutils.update_mpid(x,rsev,self.mpid_results)

            # apply 1e1p cuts
            passes = cutdefinitions.precuts(x,self.DATARUN) and cutdefinitions.postcuts(x,self.DATARUN)

            # apply BDT score
            if x.BDTscore_1e1p<BDTCUT or x.Enu_1e1p>500.0:
                passes = False

            # for debug: make something pass in order to check
            if self._DEBUG_MODE_:
                passes = True # for debug

            print "[first pass signal] RSE=",rse," RSEV=",rsev," Passes=",passes
            print "  precuts: ",passprecuts==1
            print "  simplecuts: ",x.PassSimpleCuts==1
            print "  showerreco: ",x.PassShowerReco==1
            print "  electron edep: ",x.Electron_Edep>35.0," (",x.Electron_Edep,")"
            print "  proton edep: ",x.Proton_Edep>60.0," (",x.Proton_Edep,")"
            print "  maxshrfrac: ",max(x.MaxShrFrac,-1)>0.2," (",x.MaxShrFrac,")"
            print "  pi0mass: ",x.Pi0Mass<=50.0," (",x.Pi0Mass,")"
            print "  enu: ",x.Enu_1e1p
            print "  bdt 1e1p>",BDTCUT,": ",bdtscore_1e1p>BDTCUT," (",bdtscore_1e1p,") [rerun=",self.rerun_1e1p_bdt,"]"

            if rse not in max_rse:
                # provide default
                print "  [[ first passing vtx candidates for RSE ]]"
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":-1.0,"enu":bdtscore_1e1p,"passes":False}

            if passes and max_rse[rse]["bdt"]<bdtscore_1e1p:
                # update max bdt for vertex that passes
                print "  [[ new max bdt score for passing vtx candidates ]]"
                max_rse[rse] = {"vtxid":dlanatree.vtxid,"bdt":bdtscore_1e1p,"enu":dlanatree.Enu_1e1p,"passes":True}

        print "[ dl_multi_filter::run_1e1p_signal_filter ] cutting pass. entries in max_rse=",len(max_rse)
        # next, save only those events, whose highest bdt score pass threshold
        rse_dict = {}
        rsev_dict = {}

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

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

            if not max_rse[rse]["passes"] or max_rse[rse]["bdt"]<BDTCUT:
                print "[signal] RSE ",rse,": below bdt=",max_rse[rse]["bdt"]," or did not pass (",max_rse[rse]["passes"],")"
                continue

            print "[signal] RSE ",rse,": enu=",max_rse[rse]["enu"]," bdt=",max_rse[rse]["bdt"]," PASSES signal selection (",max_rse[rse]["passes"],")"
            passes = True

            # update event flag
            if rse not in rse_dict:
                rse_dict[rse]   = passes
            elif rse in rse_dict and passes:
                rse_dict[rse]   = passes

            # update vertex flag
            rsev_dict[rsev] = passes
            
        print "[signal] number of passing RSEV: ",len(rsev_dict)
        return rse_dict, rsev_dict


    def run_rse_filter(self, dlanatree ):
        """ use a (run,subrun,event) list to filter """
        print "[ dl_multi_filter::RSE filter ]"
        rse_dict = {}
        rsev_dict = {}

        flist = open(self.filter_pars["rse-list"],'r')
        llist = flist.readlines()
        rse_list = []
        for l in llist:
            info = l.strip().split()
            rse = (int(info[0]),int(info[1]),int(info[2]))
            rse_list.append(rse)
        rse_list.sort()
        print "[ dl_multi_filter::RSE filter ] number of (rse) in list: ",len(rse_list)

        for ientry in xrange(dlanatree.GetEntries()):

            if self.MAXENTRIES is not None and ientry+1>=self.MAXENTRIES:
                break

            dlanatree.GetEntry(ientry)

            passes = False
            rse  = (dlanatree.run,dlanatree.subrun,dlanatree.event)
            rsev = (dlanatree.run,dlanatree.subrun,dlanatree.event,dlanatree.vtxid)
            if rse in rse_list:
                rse_dict[rse] = True
                rsev_dict[rsev] = True
            elif rse in rse_dict:
                pass
            else:
                rse_dict[rse] = False


        return rse_dict, rsev_dict

    def get_tree_lists(self,input_file,nvertex_entries,nevent_entries):
        """
        Get list of trees in the input (DLANA) file. Split them up into certain categories
        * vertex-indexed
        * event-indexed
        * other trees
        * pot summary tree
        * THE final vertex variable tree
        """
        vertex_indexed_trees = []
        event_indexed_trees  = []
        other_trees = []
        pot_sum_tree = None
        fvv_tree = None

        dirlist = [None]
        dirdict = {}
        
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
                        pot_sum_tree = atree
                        continue

                    if is_tree_event_indexed( basetreename ):
                        event_indexed_trees.append(atree)
                        if nevent_entries>0 and nentries!=nevent_entries and basetreename!="ubdlana_id_tree":
                            raise RuntimeError("Event-indexed tree({}) nentries ({}) does not seem to match known event-indexed tree entries ({})".format(basetreename,nentries,nevent_entries))
                    else:
                        vertex_indexed_trees.append(atree)
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

        return vertex_indexed_trees,event_indexed_trees,other_trees,fvv_tree,pot_sum_tree,dirdict
        
    def run_all_filters( self, finalvertextree, filter_list ):
        """
        run the filters, collection dictionary of results
        return dict: keys=filter name, value=result dictionary list

        result dictionary list contains
        (run,subrun,event) dictionary indiciating which passes
        (run,subrun,event,vertexid) dictionary indicating which passed
        """
        filters_results = {}
        for filtertype in filter_list:

            if filtertype=="numu-sideband":
                rse,rsev = self.run_numu_filter(finalvertextree)
            elif filtertype=="1e1p-highE-sideband":
                rse,rsev = self.run_1e1p_highE_filter(finalvertextree)
            elif filtertype=="1e1p-nearE-sideband":
                rse,rsev = self.run_1e1p_nearE_filter(finalvertextree)
            elif filtertype=="1e1p-far-sideband":
                rse,rsev = self.run_1e1p_farsideband_filter(finalvertextree)
            elif filtertype=="1e1p-midBDT-sideband":
                rse,rsev = self.run_1e1p_midBDT_filter(finalvertextree)
            elif filtertype=="1e1p-signal":
                rse,rsev = self.run_1e1p_signal_filter(finalvertextree, bdtcut=0.95)
            elif filtertype=="1e1p-loose-signal":
                rse,rsev = self.run_1e1p_signal_filter(finalvertextree, bdtcut=0.7)
            elif filtertype=="pi0-lowBDT-sideband":
                rse,rsev = self.run_lowBDT_pi0_filter(finalvertextree)
            elif filtertype=="rse-list":
                rse,rsev = self.run_rse_filter(finalvertextree)
            else:
                raise ValueError("unrecognized filter type: ",filtertype)
            
            filters_results[filtertype] = {"rse":rse,"rsev":rsev}

        return filters_results

