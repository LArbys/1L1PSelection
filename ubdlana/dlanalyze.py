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

# PRECUT utils
from LEEPreCuts_Functions import makePMTpars,performPMTPrecuts,getPMTPrecutDict

# MPID utils
import mpidutil

# DL Final Vertex Variables
from showerrecodata import ShowerRecoData # provides wrapper for shower reco info
from dlanatree import DLanaTree
from makeselectionvars import make_selection_vars

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

        self.input_file_list = [] # gets append to during open_input function
        self.tracker_treename = config['modules']['dlanalyze']['tracker_tree']
        self.ismc             = config['modules']['dlanalyze']['ismc']
        self.sample_type      = config['modules']['dlanalyze']['sample_type']

        another_tree = config['modules']['dlanalyze']['another_tree']
        print 'DLAnalyze constructed with second tree = %s' % another_tree
        self.tree_name = "_recoTree"
        self.tree_obj = None

        # SETUP SSNET Shower reco
        self.llout_name = config['modules']['dlanalyze']['showerreco']['out_larlite_tree']
        shr_ana         = config['modules']['dlanalyze']['showerreco']['out_ana_tree']
        self.adc_tree   = config['modules']['dlanalyze']['showerreco']['adctree']
        self.second_shr = config['modules']['dlanalyze']['showerreco']['second_shower']
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

        # SETUP MPID
        mpid_cfg = os.environ["UBMPIDNET_DIR"]+"/production_cfg/inference_config_tufts_WC.cfg"
        self.mpid, self.mpid_cfg = mpidutil.load_mpid_model( mpid_cfg )

        # SETUP PMT Precuts module
        self.precutpars = makePMTpars( self.sample_type )
        self.precutpars['ophittree'] = config['modules']['dlanalyze']['precut_ophits']

        # SCE class
        self.sce = larutil.SpaceChargeMicroBooNEMCC9()

        # Calibration maps
        calibmap_path = os.environ["UBDLANA_DIR"]+"/CalibrationMaps_MCC9.root"
        self.calibfile = TFile.Open(calibmap_path,'read')
        calibMap0 = self.calibfile.Get("hImageCalibrationMap_00")
        calibMap1 = self.calibfile.Get("hImageCalibrationMap_01")
        calibMap2 = self.calibfile.Get("hImageCalibrationMap_02")
        self.calibMap_v = [calibMap0,calibMap1,calibMap2]

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
        """ this is a loop over vertex entries """
        entry = tree.GetReadEntry()
        self.tree_obj.GetEntry(entry)

        print "----- ANALYZE ENTRY [%d] ---------------------------------"%(entry)
        #for branch in self.tree_obj.GetListOfBranches():
        #    if self.tree_obj.GetBranchStatus(branch.GetName()):
        #        print '  %s' % branch.GetName()

        make_selection_vars( entry, self.ismc,
                             self.tree_obj, self.df_ShowerReco, self.PMTPrecut_Dict, self.MC_dict,
                             self.anatreeclass, self.calibMap_v,
                             sce = self.sce )
        self.anatreeclass.outTree.Fill()

    def analyze_larlite_and_larcv_entry(self, tree, entry):
        #----------------------------------------------------------------------
        #
        # Purpose: Analyze loaded tree (fill histograms).  Called by framework.
        #
        # Arguments: tree - Loaded tree.
        #
        #----------------------------------------------------------------------

        print 'analyze_larlite_and_larcv_entry called for entry %d.' % entry

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

        # run mpid reco
        nmpid_vertices = mpidutil.run_mpid_on_larcv_entry( self.mpid_cfg, self.mpid, self.in_lcv, self.mpid_data, self.mpid_anatree )


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
        self.input_file_list.append(input_file.GetName())

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
        print 'Looking for TTree named %s' % self.tracker_treename
        obj = input_file.Get(self.tracker_treename)
        if obj.InheritsFrom('TTree'):
            print 'Found %s with %d entries.' % (self.tracker_treename, obj.GetEntriesFast())

            # Activate all branches.

            obj.SetBranchStatus('*', 1)
            #print 'List of activated branches of tree %s' % self.tree_name
            #for branch in obj.GetListOfBranches():
            #    if obj.GetBranchStatus(branch.GetName()):
            #        print '  %s' % branch.GetName()

            self.tree_obj = obj
            print

        # LARLITE ID TREE: governs loop to make MPID, Showerreco and Precut variables
        self.larlite_id_tree = input_file.Get("larlite_id_tree")

        # make mpid and showerreco data
        self.make_morereco_variables()

        # make precut data
        self.PMTPrecut_Dict = performPMTPrecuts( input_file.GetName(), **self.precutpars )

        # OTHER SELECTION VAR TREES, WE FRIEND THEM TO THE TRACKER TREE
        self.VtxTree  = input_file.Get("VertexTree")
        self.ShpTree  = input_file.Get("ShapeAnalysis")
        self.ShrTree  = input_file.Get("SecondShowerAnalysis")
        if self.ismc:
            self.MCTree = input_file.Get("MCTree")
            self.MC_dict = {}
            for ev in self.MCTree:
                run            = ev.run
                subrun         = ev.subrun
                event          = ev.event
                IDev           = tuple((run,subrun,event))
                self.MC_dict[IDev] = dict(parentPDG=ev.parentPDG,energyInit=ev.energyInit,
                                          parentX=ev.parentX,parentY=ev.parentY,parentZ=ev.parentZ,
                                          nproton=ev.nproton,nlepton=ev.nlepton)
        else:
            self.MCTree = None

        self.tree_obj.AddFriend(self.VtxTree)
        self.tree_obj.AddFriend(self.ShpTree)
        self.tree_obj.AddFriend(self.ShrTree)
        if self.ismc:
            self.tree_obj.AddFriend(self.MCTree)

        # We attach the MPID ana tree as well
        self.tree_obj.AddFriend(self.mpid_anatree)


        print "[ End input tree prep ]"
        print "================================"

    def make_morereco_variables(self):
        """ before running event loop, we go through events and make showerreco, mpid, precut variables """
        print "make showerreco, mpid, precut variables"
        nentries = self.larlite_id_tree.GetEntries()

        # we create a dictionary where we will store shower reco variables
        self.dict_ShowerReco = {"entries":[]}

        for entry in xrange(nentries):
            self.larlite_id_tree.GetEntry(entry)
            self.analyze_larlite_and_larcv_entry(self.larlite_id_tree,entry)
            self.extract_showerreco_entry_vars(self.dict_ShowerReco,self.showerreco)

        self.df_ShowerReco = ShowerRecoData( self.dict_ShowerReco )

    def end_job(self):
        """ close larcv and larlite files. larlite output file will write."""
        self.in_lcv.finalize()
        self.io_ll.close()
        self.calibfile.Close()
        fout = open('dlanalyze_input_list.txt','w')
        for f in self.input_file_list:
            print>>fout,f.strip()
        fout.close()


    def extract_showerreco_entry_vars(self,data,showerreco):
        """ get variables from showerreco class """

        entrydata = { "run":self.in_lcv.event_id().run(),
                      "subrun":self.in_lcv.event_id().subrun(),
                      "event":self.in_lcv.event_id().event(),
                      "shower_energies":[],
                      "shower_sumQs":[],
                      "shower_shlengths":[],
                      "vertex_pos":[],
                      "shower_gap":[],
                      "shower_direction_3d":[],
                      "shower_direction_2d":[],
                      "shower_opening_2d":[],
                      "shower_start_2d":[],
                      }

        # Save first shower output
        for ivtx in xrange(showerreco.numVertices()):
            entrydata["shower_energies"].append( [ showerreco.getVertexShowerEnergy(ivtx,p) for p in xrange(3) ] )
            entrydata["shower_sumQs"].append( [ showerreco.getVertexShowerSumQ(ivtx,p) for p in xrange(3) ] )
            entrydata["shower_shlengths"].append( [ showerreco.getVertexShowerShlength(ivtx,p) for p in xrange(3) ] )
            entrydata["vertex_pos"].append( [ showerreco.getVertexPos(ivtx).at(p) for p in xrange(3) ] )
            entrydata["shower_gap"].append( [ showerreco.getVertexShowerGap(ivtx,p) for p in xrange(3) ] )
            entrydata["shower_direction_3d"].append( [ showerreco.getFirstDirection(ivtx,dir) for dir in xrange(3) ])
            entrydata["shower_direction_2d"].append( [ showerreco.getVertexShowerDirection2D(ivtx,dir) for dir in xrange(3) ])
            entrydata["shower_opening_2d"].append( [ showerreco.getVertexShowerOpening2D(ivtx,dir) for dir in xrange(3) ])
            # showerstart also needs a loop over x,y,z
            for p in xrange(3):
                entrydata["shower_start_2d"].append( [ showerreco.getShowerStart2D(ivtx,p,dir) for dir in xrange(3) ])


        if self.second_shr:
            # Save second shower output
            entrydata["secondshower_energies"] = []
            entrydata["secondshower_sumQs"] = []
            entrydata["secondshower_shlengths"] = []
            entrydata["secondshower_gap"] = []
            entrydata["pi0mass"] = []
            entrydata["opening_angle_3d"] = []
            entrydata["shower_impact"] = []
            entrydata["secondshower_impact"] =[]
            entrydata["secondshower_direction_3d"]=[]
            entrydata["secondshower_direction_2d"]=[]
            entrydata["secondshower_opening_2d"]=[]
            entrydata["secondshower_start_2d"] =[]

            for ivtx in xrange(showerreco.numVertices()):
                entrydata["secondshower_energies"].append( [ showerreco.getVertexSecondShowerEnergy(ivtx,p) for p in xrange(3) ] )
                entrydata["secondshower_sumQs"].append( [ showerreco.getVertexSecondShowerSumQ(ivtx,p) for p in xrange(3) ] )
                entrydata["secondshower_shlengths"].append( [ showerreco.getVertexSecondShowerShlength(ivtx,p) for p in xrange(3) ] )
                entrydata["secondshower_gap"].append( [ showerreco.getVertexSecondShowerGap(ivtx,p) for p in xrange(3) ] )
                entrydata["pi0mass"].append([showerreco.getPi0Mass(ivtx)])
                entrydata["opening_angle_3d"].append([showerreco.getAlpha(ivtx)])
                entrydata["shower_impact"].append([showerreco.getImpact1(ivtx)])
                entrydata["secondshower_impact"].append([showerreco.getImpact2(ivtx)])
                entrydata["secondshower_direction_2d"].append( [ showerreco.getVertexSecondShowerDirection2D(ivtx,dir) for dir in xrange(3) ])
                entrydata["secondshower_opening_2d"].append( [ showerreco.getVertexSecondShowerOpening2D(ivtx,dir) for dir in xrange(3) ])
                entrydata["secondshower_direction_3d"].append( [ showerreco.getSecondDirection(ivtx,dir) for dir in xrange(3) ])
                # showerstart also needs a loop over x,y,z
                for p in xrange(3):
                    entrydata["secondshower_start_2d"].append( [ showerreco.getSecondShowerStart2D(ivtx,p,dir) for dir in xrange(3) ])

        # Below are MC-based truth metrics

        # save vertex truth information
        if False:
            # build graph and get primary particles
            mcpg.buildgraph( self.in_lcv, ioll )
            mcpg.printGraph()
            node_v = mcpg.getPrimaryParticles()

            # determine topology, get electron energy if available
            nelectrons = 0
            nprotons   = 0
            nother     = 0
            pidX = []
            pidOther = []
            for inode in xrange(node_v.size()):
                node = node_v.at(inode)
                #print "node[",inode,"] pid=",node.pid," E=",node.E_MeV
                if abs(node.pid)==11:
                    entrydata["true_electron_energy"] = node.E_MeV
                    nelectrons += 1
                elif node.pid==2212:
                    #print "found proton: ",node.E_MeV
                    if node.E_MeV>60.0:
                        nprotons += 1
                else:
                    if node.pid in [211,-211]:
                        # charged pions
                        nother += 1
                        pidOther.append( node.pid )
                    else:
                        pidX.append(node.pid)
            entrydata["true_topology"] = "%de%dp%dX"%(nelectrons,nprotons,nother)
            print "topology",entrydata["true_topology"]
            #print "unknown PID: ",pidX

            # determine distance of reco vertices from true
            vtx_v = mcpg.findTrackID(-1).start
            entrydata["true_vertex"] = [ vtx_v[i] for i in xrange(3) ] # get the ROOT node
            offset_v = sce.GetPosOffsets( vtx_v[0], vtx_v[1], vtx_v[2] )
            vtx_sce_v = [ vtx_v[0]-offset_v[0]+0.7,
                          vtx_v[1]+offset_v[1],
                          vtx_v[2]+offset_v[2] ]
            entrydata["true_vertex_sce"] = vtx_sce_v

            entrydata["vertex_dist_from_truth"] = []
            for pos in entrydata["vertex_pos"]:
                d = 0.0
                for i in xrange(3):
                    d += (pos[i]-vtx_sce_v[i])*(pos[i]-vtx_sce_v[i])
                d = sqrt(d)
                entrydata["vertex_dist_from_truth"].append(d)

        if False:  #args.use_bnb:
            entrydata["pi0mass"] = []
            entrydata["haspi0"] = showerreco.getHasPi0()
            entrydata["ccnc"] = showerreco.getCCNC()
            entrydata["disttoint"] = []
            entrydata["impact1"] = []
            entrydata["impact2"] = []
            entrydata["alpha"] = []
            entrydata["firstdirection"] = []
            entrydata["seconddirection"] = []

            mcpg.buildgraph( self.in_lcv, ioll )
            vtx_v = mcpg.findTrackID(-1).start
            entrydata["true_vertex"] = [ vtx_v[i] for i in xrange(3) ] # get the ROOT node
            offset_v = sce.GetPosOffsets( vtx_v[0], vtx_v[1], vtx_v[2] )
            vtx_sce_v = [ vtx_v[0]-offset_v[0]+0.7,
                          vtx_v[1]+offset_v[1],
                          vtx_v[2]+offset_v[2] ]
            entrydata["true_vertex_sce"] = vtx_sce_v

            for ivtx in xrange(showerreco.numVertices()):
                entrydata["pi0mass"].append( showerreco.getPi0Mass(ivtx))
                entrydata["disttoint"].append( showerreco.getDistToInt(ivtx))
                entrydata["impact1"].append( showerreco.getImpact1(ivtx))
                entrydata["impact2"].append( showerreco.getImpact2(ivtx))
                entrydata["alpha"].append( showerreco.getAlpha(ivtx))
                for dir in xrange(3):
                    entrydata["firstdirection"].append( showerreco.getFirstDirection(ivtx,dir))
                    entrydata["seconddirection"].append( showerreco.getSecondDirection(ivtx,dir))


        if False:#args.use_ncpi0:
            entrydata["true_shower_energies"] = []
            entrydata["true_shower_starts"] = []
            entrydata["remaining_adc"] = []
            entrydata["overlap_fraction1"] = []
            entrydata["overlap_fraction2"] = []
            entrydata["purity"] = []
            entrydata["efficiency"] = []
            entrydata["pi0mass"] = []
            entrydata["useformass"] = []
            entrydata["disttoint"] = []
            entrydata["impact1"] = []
            entrydata["impact2"] = []

            for ivtx in xrange(showerreco.numVertices()):
                entrydata["useformass"].append( showerreco.getUseForMass(ivtx))
                entrydata["pi0mass"].append( showerreco.getPi0Mass(ivtx))
                entrydata["disttoint"].append( showerreco.getDistToInt(ivtx))
                entrydata["impact1"].append( showerreco.getImpact1(ivtx))
                entrydata["impact2"].append( showerreco.getImpact2(ivtx))
                entrydata["overlap_fraction1"].append( [showerreco.getOverlapFraction1(ivtx,plane,0) for plane in xrange(2) ] )
                entrydata["overlap_fraction1"].append( [showerreco.getOverlapFraction1(ivtx,plane,1) for plane in xrange(2) ] )
                entrydata["overlap_fraction2"].append( [showerreco.getOverlapFraction2(ivtx,plane,0) for plane in xrange(2) ] )
                entrydata["overlap_fraction2"].append( [showerreco.getOverlapFraction2(ivtx,plane,1) for plane in xrange(2) ] )
                entrydata["purity"].append( [showerreco.getShowerTruthMatchPur(ivtx,shower) for shower in xrange(6)])
                entrydata["efficiency"].append( [showerreco.getShowerTruthMatchEff(ivtx,shower) for shower in xrange(6)])


            for ivtx in xrange(showerreco.numShowers()):
                entrydata["true_shower_energies"].append( [ showerreco.getTrueShowerEnergy(ivtx) for shower in xrange(2) ] )
                entrydata["true_shower_starts"].append( [ showerreco.getTrueShowerStarts(ivtx).at(p) for p in xrange(3) ] )
                entrydata["remaining_adc"].append( [showerreco.getRemainingADC()])

        if False: #args.use_nueint:
            entrydata["uplane_profile"] = []
            entrydata["vplane_profile"] = []
            entrydata["yplane_profile"] = []
            entrydata["purity"] = []
            entrydata["efficiency"] = []

            for ivtx in xrange(showerreco.numVertices()):
                entrydata["purity"].append( [showerreco.getShowerTruthMatchPur(ivtx,shower) for shower in xrange(3)])
                entrydata["efficiency"].append( [showerreco.getShowerTruthMatchEff(ivtx,shower) for shower in xrange(3)])

            for ii in xrange(showerreco.numpointsU()):
                entrydata["uplane_profile"].append( [ showerreco.getUPlaneShowerProfile(ii,index) for index in xrange(2) ] )
            for ii in xrange(showerreco.numpointsV()):
                entrydata["vplane_profile"].append( [ showerreco.getVPlaneShowerProfile(ii,index) for index in xrange(2) ] )
            for ii in xrange(showerreco.numpointsY()):
                entrydata["yplane_profile"].append( [ showerreco.getYPlaneShowerProfile(ii,index) for index in xrange(2) ] )

        data["entries"].append( entrydata )

        return
