#v1: added pmt precut stuff
#v2: added mctruth info (where applicable, obvi)
#v3: added energy calibration to dqdx
#v4: added mpidtree
#v5: added new shower reco from nick

# TODO:
# - make sure CCQE energies match proper off-shell formulae  (talk to steven about that)

import os,sys
from sys import argv

import ROOT
from ROOT import TFile,TTree
import pickle
import numpy as np
import pandas as pd
from numpy import mean,asarray,matmul
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from sys import argv
from array import array
from larlite import larlite,larutil
import os,sys

from SelectionDefs import NewAng, VtxInSimpleFid, VtxInFid, GetPhiT, pTrans,pTransRat, alphaT, ECCQE, ECal, Q2, OpenAngle, PhiDiff, edgeCut, ECCQE_mom, Boost, BoostTracks, Getpz, GetCCQEDiff, SensibleMinimize, Getq3q0,GetTotPE, CorrectionFactor

if len(argv) != 6:
    print('Fuck off')
    print('argv[1]: dlmerged.root')
    print('argv[2]: calibration map')
    print('argv[3]: mpid')
    print('argv[4]: showerreco')
    print('argv[5]: destination (.)')
    sys.exit()

_tag = argv[1][-27:-5]
_dest = argv[5]
_dlmerged = argv[1]
_mpid = argv[3]
_calibmap = argv[2]

sce = larutil.SpaceChargeMicroBooNEMCC9()

# --------------------------------------------- #

def MakeTreeBranch(ttree,s_name,s_type):
    if s_type == 'int':
        _myvar = array('i',[-9999])
        ttree.Branch(s_name,_myvar,s_name+'/I')
        return _myvar

    if s_type == 'float':
        _myvar = array('f',[-9999])
        ttree.Branch(s_name,_myvar,s_name+'/F')
        return _myvar

    if s_type == 'tvector':
        _myvar = ROOT.vector('double')()
        ttree.Branch(s_name,_myvar)
        return _myvar



def PerformPMTPrecuts(s_dlmerged):
    """
    calculate PMT precuts
    """
    
    # ---------------------------------------------- #
    #print('<EVID: %s> -- First, we will figure out the PMT Precut info.'%_tag)  #gotta do this first for io reasons    
    PMTPrecut_Dict = {}

    PP_WINDOW_LENGTH = 130
    PP_COINC_WINDOW = 6
    PP_PE_THRESH = 20
    PP_PMT_MAX_FRAC_CUTOFF = 0.60
    PP_WIN_START = 190
    PP_WIN_END = 320
    PP_PORCH_WIN_START = PP_WIN_START - 130
    PP_PORCH_WIN_END = PP_WIN_END - 130
    PP_TICK_SIZE = 0.015625		# 15.625 ns per tick

    # Load up larlite
    ll_manager = larlite.storage_manager()
    ll_manager.set_io_mode(ll_manager.kREAD)
    ll_manager.add_in_filename(argv[1])
    ll_manager.set_in_rootdir("")
    ll_manager.open()

    while ll_manager.next_event():
    
        id_rse = tuple((ll_manager.run_id(),ll_manager.subrun_id(),ll_manager.event_id()))

        # ---------- Grab vectors containing particle info --------- #
        ophits   = ll_manager.get_data(larlite.data.kOpHit,"ophitBeam")
        # ---------------------------------------------------------- #

        # Torch cut and PE cut
        porchClusters = [0]*((PP_WIN_END-PP_WIN_START)/PP_COINC_WINDOW+1)
        Clusters = [0]*((PP_WIN_END-PP_WIN_START)/PP_COINC_WINDOW+1)
        for ophit in ophits:
            if PP_PORCH_WIN_START <= ophit.PeakTime()/PP_TICK_SIZE <= PP_PORCH_WIN_END:
                porchClusters[int((ophit.PeakTime()/PP_TICK_SIZE-PP_PORCH_WIN_START)/PP_COINC_WINDOW)]+=ophit.PE()
            if PP_WIN_START <= ophit.PeakTime()/PP_TICK_SIZE <= PP_WIN_END:
                Clusters[int((ophit.PeakTime()/PP_TICK_SIZE-PP_WIN_START)/PP_COINC_WINDOW)]+=ophit.PE()

        porchTotPE,porchFlashBins = GetTotPE(PP_PE_THRESH,porchClusters)
        TotPE,FlashBins = GetTotPE(PP_PE_THRESH,Clusters)

        # PMT Frac Cut
        PEbyPMT  = [0]*32
        FracDict = {}
        for ophit in ophits:
            if (int(ophit.PeakTime()/PP_TICK_SIZE - PP_WIN_START)/PP_COINC_WINDOW) in FlashBins:
                PEbyPMT[ophit.OpChannel()]+=ophit.PE()
                if ophit.OpChannel() in FracDict:
                    FracDict[ophit.OpChannel()]+=ophit.PE()/TotPE
                else:
                    FracDict[ophit.OpChannel()] = ophit.PE()/TotPE

        if TotPE >= PP_PE_THRESH:
            PMTfrac = [x / TotPE for x in PEbyPMT]
        else:
            PMTfrac = []

        PMTPrecutDict[id_rse] = dict(_totpe=TotPE,_porchtotpe=porchTotPE,_maxpefrac=max(PMTfrac),_passpmtprecut=0)

    return PMTPrecutDict
    
def GetPMTPrecutDict(s_dlmerged):
    """ get results from previously applied pmt precuts """
    PMTPrecutDict = {}

    # Load up larlite
    ll_manager = larlite.storage_manager()
    ll_manager.set_io_mode(ll_manager.kREAD)
    ll_manager.add_in_filename(s_dlmerged)
    ll_manager.set_data_to_read(larlite.data.kUserInfo,'precutresults')
    ll_manager.set_in_rootdir("")
    ll_manager.open()

    while ll_manager.next_event():
        id_rse = tuple((ll_manager.run_id(),ll_manager.subrun_id(),ll_manager.event_id()))

        # ---------- Grab vectors containing particle info --------- #
        precutresults   = ll_manager.get_data(larlite.data.kUserInfo,"precutresults")
        precutresult = precutresults[0]
        # ---------------------------------------------------------- #
        PMTPrecutDict[id_rse] = dict(_totpe=precutresult.get_double("beamPE"),
                                     _porchtotpe=precutresult.get_double("vetoPE"),
                                     _maxpefrac=precutresult.get_double("maxFrac"),
                                     _passpmtprecut=precutresult.get_int("pass"))

    ll_manager.close()
    return PMTPrecutDict

# -----------------------------------------------------------------------#
#    Get started!
# -----------------------------------------------------------------------#

print()
print('<EVID: %s> -- First, we will figure out the PMT Precut info.'%_tag)  #gotta do this first for io reasons
PMTPrecut_Dict = GetPMTPrecutDict(_dlmerged)
PMTPrecut_Dict = PerformPMTPrecuts(_dlmerged)

print()
print('<EVID: %s> -- Now read in the shower reco stuff.'%_tag)
df_ShowerReco = pd.read_csv(argv[4],sep=' ')
df_ShowerReco.set_index(['run','subrun','event','vtxid'],inplace=True)

print()
print('<EVID: %s> -- Now make sure we can read the root trees we want (and make sure they are present).'%_tag)
try:
    DLMergedFile = TFile(_dlmerged,'read')
    MPIDFile = TFile(_mpid,'read')

    TrkTree  = DLMergedFile.Get("_recoTree")
    VtxTree  = DLMergedFile.Get("VertexTree")
    ShpTree  = DLMergedFile.Get("ShapeAnalysis")
    ShrTree  = DLMergedFile.Get("SecondShowerAnalysis")
    MCTree   = DLMergedFile.Get("MCTree")
    MPIDTree = MPIDFile.Get("multipid_tree")

    TrkTree.AddFriend(VtxTree)
    TrkTree.AddFriend(ShpTree)
    TrkTree.AddFriend(ShrTree)
    TrkTree.AddFriend(MPIDTree)
except:
    print 'Fucked: %s %s' %(_dlmerged,_mpid)
    sys.exit()

print()
print('<EVID: %s> -- Load up that truth info.'%_tag)
MC_dict = {}
IsMC = False
for ev in MCTree:
    run            = ev.run
    subrun         = ev.subrun
    event          = ev.event
    IDev           = tuple((run,subrun,event))

    MC_dict[IDev] = dict(parentPDG=ev.parentPDG,energyInit=ev.energyInit,parentX=ev.parentX,parentY=ev.parentY,parentZ=ev.parentZ,nproton=ev.nproton,nlepton=ev.nlepton)
    if not np.isnan(ev.energyInit):
        IsMC = True

print()
print('<EVID: %s> -- Load up calibration map.'%_tag)
calibfile = TFile.Open(_calibmap,'read')
calibMap0 = calibfile.Get("hImageCalibrationMap_00")
calibMap1 = calibfile.Get("hImageCalibrationMap_01")
calibMap2 = calibfile.Get("hImageCalibrationMap_02")
calibMap_v = [calibMap0,calibMap1,calibMap2]

print()
print('<EVID: %s> -- Create target root file and get started!'%_tag)
outFileName = 'FinalVertexVariables-prime_'+_tag+'.root'
outFileName = os.path.join(_dest,outFileName)
outFile = TFile(outFileName,'RECREATE')
outTree = TTree('FinalVertexVariables','Final Vertex Variable Tree')

## lepton agnostic
_run = MakeTreeBranch(outTree,'run','int')
_subrun = MakeTreeBranch(outTree,'subrun','int')
_event = MakeTreeBranch(outTree,'event','int')
_vtxid = MakeTreeBranch(outTree,'vtxid','int')
_x = MakeTreeBranch(outTree,'Xreco','float')
_y = MakeTreeBranch(outTree,'Yreco','float')
_z = MakeTreeBranch(outTree,'Zreco','float')
_infiducial = MakeTreeBranch(outTree,'InFiducial','int')
_anyReco = MakeTreeBranch(outTree,'AnyReco','int')
_ntracks = MakeTreeBranch(outTree,'NTracks','int')
_n5tracks = MakeTreeBranch(outTree,'N5cmTracks','int')
_passCuts = MakeTreeBranch(outTree,'PassSimpleCuts','int')
_passShowerReco = MakeTreeBranch(outTree,'PassShowerReco','int')
_passSecShr = MakeTreeBranch(outTree,'PassSecondShower','int')
_failBoost = MakeTreeBranch(outTree,'FailedBoost','int')
_good3DReco = MakeTreeBranch(outTree,'Good3DReco','int')
_eta = MakeTreeBranch(outTree,'Eta','float')
_openAng = MakeTreeBranch(outTree,'OpenAng','float')
_thetas = MakeTreeBranch(outTree,'Thetas','float')
_phis = MakeTreeBranch(outTree,'Phis','float')
_charge_near_trunk = MakeTreeBranch(outTree,'ChargeNearTrunk','float')
_longtracklen = MakeTreeBranch(outTree,'LongTrackLen','float')
_shorttracklen = MakeTreeBranch(outTree,'ShortTrackLen','float')
_maxshrFrac = MakeTreeBranch(outTree,'MaxShrFrac','float')
_minshrFrac = MakeTreeBranch(outTree,'MinShrFrac','float')

# 1mu1p stuff
_CCQE_energy_shift_1m1p = MakeTreeBranch(outTree,'CCQEEnergyShift_1m1p','float')
_enu_1m1p = MakeTreeBranch(outTree,'Enu_1m1p','float')
_phiT_1m1p = MakeTreeBranch(outTree,'PhiT_1m1p','float')
_alphaT_1m1p = MakeTreeBranch(outTree,'AlphaT_1m1p','float')
_pT_1m1p = MakeTreeBranch(outTree,'PT_1m1p','float')
_pTRat_1m1p = MakeTreeBranch(outTree,'PTRat_1m1p','float')
_bjX_1m1p = MakeTreeBranch(outTree,'BjX_1m1p','float')
_bjY_1m1p = MakeTreeBranch(outTree,'BjY_1m1p','float')
_q2_1m1p = MakeTreeBranch(outTree,'Q2_1m1p','float')
_sph_1m1p = MakeTreeBranch(outTree,'Sph_1m1p','float')
_pzEnu_1m1p = MakeTreeBranch(outTree,'PzEnu_1m1p','float')
_q0_1m1p = MakeTreeBranch(outTree,'Q0_1m1p','float')
_q3_1m1p = MakeTreeBranch(outTree,'Q3_1m1p','float')
_openAngB_1m1p = MakeTreeBranch(outTree,'OpenAngB_1m1p','float')
_thetasB_1m1p = MakeTreeBranch(outTree,'ThetasB_1m1p','float')
_phisB_1m1p = MakeTreeBranch(outTree,'PhisB_1m1p','float')
_phiTB_1m1p = MakeTreeBranch(outTree,'PhiTB_1m1p','float')
_alphaTB_1m1p = MakeTreeBranch(outTree,'AlphaTB_1m1p','float')
_pTB_1m1p = MakeTreeBranch(outTree,'PTB_1m1p','float')
_bjXB_1m1p = MakeTreeBranch(outTree,'BjXB_1m1p','float')
_bjYB_1m1p = MakeTreeBranch(outTree,'BjYB_1m1p','float')
_q2B_1m1p = MakeTreeBranch(outTree,'Q2B_1m1p','float')
_sphB_1m1p = MakeTreeBranch(outTree,'SphB_1m1p','float')

# same stuff, but for 1e1p
_CCQE_energy_shift_1e1p = MakeTreeBranch(outTree,'CCQEEnergyShift_1e1p','float')
_enu_1e1p = MakeTreeBranch(outTree,'Enu_1e1p','float')
_phiT_1e1p = MakeTreeBranch(outTree,'PhiT_1e1p','float')
_alphaT_1e1p = MakeTreeBranch(outTree,'AlphaT_1e1p','float')
_pT_1e1p = MakeTreeBranch(outTree,'PT_1e1p','float')
_pTRat_1e1p = MakeTreeBranch(outTree,'PTRat_1e1p','float')
_bjX_1e1p = MakeTreeBranch(outTree,'BjX_1e1p','float')
_bjY_1e1p = MakeTreeBranch(outTree,'BjY_1e1p','float')
_q2_1e1p = MakeTreeBranch(outTree,'Q2_1e1p','float')
_sph_1e1p = MakeTreeBranch(outTree,'Sph_1e1p','float')
_pzEnu_1e1p = MakeTreeBranch(outTree,'PzEnu_1e1p','float')
_q0_1e1p = MakeTreeBranch(outTree,'Q0_1e1p','float')
_q3_1e1p = MakeTreeBranch(outTree,'Q3_1e1p','float')
_openAngB_1e1p = MakeTreeBranch(outTree,'OpenAngB_1e1p','float')
_thetasB_1e1p = MakeTreeBranch(outTree,'ThetasB_1e1p','float')
_phisB_1e1p = MakeTreeBranch(outTree,'PhisB_1e1p','float')
_phiTB_1e1p = MakeTreeBranch(outTree,'PhiTB_1e1p','float')
_alphaTB_1e1p = MakeTreeBranch(outTree,'AlphaTB_1e1p','float')
_pTB_1e1p = MakeTreeBranch(outTree,'PTB_1e1p','float')
_bjXB_1e1p = MakeTreeBranch(outTree,'BjXB_1e1p','float')
_bjYB_1e1p = MakeTreeBranch(outTree,'BjYB_1e1p','float')
_q2B_1e1p = MakeTreeBranch(outTree,'Q2B_1e1p','float')
_sphB_1e1p = MakeTreeBranch(outTree,'SphB_1e1p','float')

_lepton_id = MakeTreeBranch(outTree,'Lepton_ID','int')
_lepton_phi = MakeTreeBranch(outTree,'Lepton_PhiReco','float')
_lepton_theta = MakeTreeBranch(outTree,'Lepton_ThetaReco','float')
_lepton_length = MakeTreeBranch(outTree,'Lepton_TrackLength','float')
_lepton_dqdx_uncalibrated = MakeTreeBranch(outTree,'Lepton_dQdx_uncalibrated','float')
_lepton_dqdx_calibrated = MakeTreeBranch(outTree,'Lepton_dQdx','float')
_muon_E = MakeTreeBranch(outTree,'Muon_Edep','float')
_electron_E = MakeTreeBranch(outTree,'Electron_Edep','float')
_lepton_edge_dist = MakeTreeBranch(outTree,'Lepton_EdgeDist','float')
_muon_phiB_1m1p = MakeTreeBranch(outTree,'Muon_PhiRecoB_1m1p','float')
_muon_thetaB_1m1p = MakeTreeBranch(outTree,'Muon_ThetaRecoB_1m1p','float')
_muon_EB_1m1p = MakeTreeBranch(outTree,'Muon_EdepB_1m1p','float')
_electron_phiB_1e1p = MakeTreeBranch(outTree,'Electron_PhiRecoB_1e1p','float')
_electron_thetaB_1e1p = MakeTreeBranch(outTree,'Electron_ThetaRecoB_e1ep','float')
_electron_EB_1e1p = MakeTreeBranch(outTree,'Electron_EdepB_1e1p','float')
_proton_id = MakeTreeBranch(outTree,'Proton_ID','float')
_proton_phi = MakeTreeBranch(outTree,'Proton_PhiReco','float')
_proton_theta = MakeTreeBranch(outTree,'Proton_ThetaReco','float')
_proton_length = MakeTreeBranch(outTree,'Proton_TrackLength','float')
_proton_dqdx_uncalibrated = MakeTreeBranch(outTree,'Proton_dQdx_uncalibrated','float')
_proton_dqdx_calibrated = MakeTreeBranch(outTree,'Proton_dQdx','float')
_proton_E = MakeTreeBranch(outTree,'Proton_Edep','float')
_proton_edge_dist = MakeTreeBranch(outTree,'Proton_EdgeDist','float')
_proton_phiB_1m1p = MakeTreeBranch(outTree,'Proton_PhiRecoB_1m1p','float')
_proton_thetaB_1m1p = MakeTreeBranch(outTree,'Proton_ThetaRecoB_1m1p','float')
_proton_EB_1m1p = MakeTreeBranch(outTree,'Proton_EdepB_1m1p','float')
_proton_phiB_1e1p = MakeTreeBranch(outTree,'Proton_PhiRecoB_1e1p','float')
_proton_thetaB_1e1p = MakeTreeBranch(outTree,'Proton_ThetaRecoB_1e1p','float')
_proton_EB_1e1p = MakeTreeBranch(outTree,'Proton_EdepB_1e1p','float')

# Precut stuff
_totPE = MakeTreeBranch(outTree,'TotPE','float')
_porchTotPE = MakeTreeBranch(outTree,'PorchTotPE','float')
_maxPEFrac = MakeTreeBranch(outTree,'MaxPEFrac','float')
_passPMTPrecut = MakeTreeBranch(outTree,'PassPMTPrecut','int')

# MC stuff
_parentPDG = MakeTreeBranch(outTree,'MC_parentPDG','int')
_energyInit = MakeTreeBranch(outTree,'MC_energyInit','float')
_parentX = MakeTreeBranch(outTree,'MC_parentX','float')
_parentY = MakeTreeBranch(outTree,'MC_parentY','float')
_parentZ = MakeTreeBranch(outTree,'MC_parentZ','float')
_nproton = MakeTreeBranch(outTree,'MC_nproton','int')
_nlepton = MakeTreeBranch(outTree,'MC_nlepton','int')
_parentSCEX = MakeTreeBranch(outTree,'MC_parentSCEX','float')
_parentSCEY = MakeTreeBranch(outTree,'MC_parentSCEX','float')
_parentSCEZ = MakeTreeBranch(outTree,'MC_parentSCEX','float')
_scedr = MakeTreeBranch(outTree,'MC_scedr','float')

# MPID stuff
_eminusPID_int_v = MakeTreeBranch(outTree,'EminusPID_int_v','tvector')
_muonPID_int_v = MakeTreeBranch(outTree,'MuonPID_int_v','tvector')
_protonPID_int_v = MakeTreeBranch(outTree,'ProtonPID_int_v','tvector')
_gammaPID_int_v = MakeTreeBranch(outTree,'GammaPID_int_v','tvector')
_pionPID_int_v = MakeTreeBranch(outTree,'PionPID_int_v','tvector')
_eminusPID_pix_v = MakeTreeBranch(outTree,'EminusPID_pix_v','tvector')
_muonPID_pix_v = MakeTreeBranch(outTree,'MuonPID_pix_v','tvector')
_protonPID_pix_v = MakeTreeBranch(outTree,'ProtonPID_pix_v','tvector')
_gammaPID_pix_v = MakeTreeBranch(outTree,'GammaPID_pix_v','tvector')
_pionPID_pix_v = MakeTreeBranch(outTree,'PionPID_pix_v','tvector')

def clear_vertex():
    _lepton_id        = int(-9999)
    _lepton_phi       = float(-9999)
    _lepton_theta     = float(-9999)
    _lepton_length    = float(-9999)
    _lepton_dqdx_uncalibrated      = float(-9999)
    _lepton_dqdx_calibrated      = float(-9999)
    _lepton_edge_dist = float(-9999)
    _muon_E         = float(-9999)
    _muon_phiB_1m1p      = float(-9999)
    _muon_thetaB_1m1p    = float(-9999)
    _muon_EB_1m1p        = float(-9999)
    _electron_E         = float(-9999)
    _electron_phiB_1e1p      = float(-9999)
    _electron_thetaB_1e1p    = float(-9999)
    _electron_EB_1e1p        = float(-9999)

    _proton_id        = int(-9999)
    _proton_phi       = float(-9999)
    _proton_theta     = float(-9999)
    _proton_length    = float(-9999)
    _proton_dqdx_uncalibrated      = float(-9999)
    _proton_dqdx_calibrated      = float(-9999)
    _proton_E         = float(-9999)
    _proton_edge_dist = float(-9999)
    _proton_phiB_1m1p      = float(-9999)
    _proton_thetaB_1m1p    = float(-9999)
    _proton_EB_1m1p        = float(-9999)
    _proton_phiB_1e1p      = float(-9999)
    _proton_thetaB_1e1p    = float(-9999)
    _proton_EB_1e1p        = float(-9999)


print()
print('<EVID: %s> -- Great. Now loop through track events and go wild'%_tag)

for indo,ev in enumerate(TrkTree):
    BE = 29.5

    run            = ev.run
    subrun         = ev.subrun
    event          = ev.event
    IDev           = tuple((run,subrun,event))
    vtxid          = ev.vtx_id
    IDvtx          = tuple((run,subrun,event,vtxid))
    vtxX           = ev.vtx_x
    vtxY           = ev.vtx_y
    vtxZ           = ev.vtx_z

    # These are for passing simple precuts
    length_v       = ev.Length_v
    InFiducial     = VtxInSimpleFid(vtxX,vtxY,vtxZ)
    NumTracks      = len(length_v)
    Num5cmTracks   = sum(1 for x in length_v if x > 5)
    NothingRecod   = ev.nothingReconstructed
    EdgeDistance   = ev.closestWall

    PassSimpleCuts = (InFiducial and NumTracks == 2 and edgeCut(EdgeDistance))
    PassShowerReco = True
    PassSecondShower = True
    FailBoost = False

    vtxPhi_v       = ev.vertexPhi
    vtxTheta_v     = ev.vertexTheta
    dqdx_v_uncalibrated = ev.Avg_Ion_v
    dqdx_v_calibrated = ev.Avg_Ion_v

    Good3DReco     = ev.GoodVertex

    # - Get ssnet and shape analysiss stuff
    sh_foundClusY = [i for i in ev.shower_frac_Y_v]
    sh_foundClusV = [i for i in ev.shower_frac_V_v]

    if 0.0 < max(sh_foundClusY) < 1.0:
        shrFrac    = max(sh_foundClusY)
        shrFracPart = sh_foundClusY
    else:
        shrFrac     = max(sh_foundClusV)
        shrFracPart = sh_foundClusV

    try:
        eE = df_ShowerReco['e_reco'][IDvtx]
        if eE < 0:
            eE = -99999
            PassShowerReco = False
    except:
        eE = -9999
        PassShowerReco = False

    # - get MC Truth info (if applicable)
    if(IsMC):
        parentX = MC_dict[IDev]['parentX']
        parentY = MC_dict[IDev]['parentY']
        parentZ = MC_dict[IDev]['parentZ']

    # - Now, the big stuff
    if PassSimpleCuts:
        lid            = int(np.argmin(ev.Avg_Ion_v))
        pid            = int(np.argmax(ev.Avg_Ion_v))
        lTh            = ev.vertexTheta[lid]
        pTh            = ev.vertexTheta[pid]
        lPh            = ev.vertexPhi[lid]
        pPh            = ev.vertexPhi[pid]
        # Add detector rotation fix
        lTh,lPh        = NewAng(lTh,lPh)
        pTh,pPh        = NewAng(pTh,pPh)
        # ----
        mE             = ev.E_muon_v[lid]
        pE             = ev.E_proton_v[pid]
        ldq = ev.IonY_5cm_v[lid]
        pdq = ev.IonY_5cm_v[pid]

        thetas         = lTh+pTh
        phis           = PhiDiff(lPh,pPh)
        EpCCQE         = ECCQE(pE,pTh,pid="proton",B=BE)
        EmCCQE         = ECCQE(mE,lTh,pid="muon",B=BE)
        if PassShowerReco: EeCCQE = ECCQE(eE,lTh,pid="electron",B=BE)

        wallProton     = EdgeDistance[pid]
        wallLepton     = EdgeDistance[lid]

        openAng        = OpenAngle(pTh,lTh,pPh,lPh)

        if IsMC:
            for i in xrange(0,len(dqdx_v_uncalibrated)):
                dqdx_v_calibrated[i] = dqdx_v_uncalibrated[i] * CorrectionFactor(vtxX,vtxY,vtxZ,vtxTheta_v[i],vtxPhi_v[i],length_v[i],calibMap_v)

        eta = abs(dqdx_v_calibrated[pid]-dqdx_v_calibrated[lid])/(dqdx_v_calibrated[pid]+dqdx_v_calibrated[lid])

        longtracklen   = max(ev.Length_v)
        shorttracklen  = min(ev.Length_v)
        maxshrFrac     = max(shrFracPart)
        minshrFrac     = min(shrFracPart)

        # - Correct SCE of coords
        if(IsMC):
            if(VtxInSimpleFid(parentX,parentY,parentZ,5.0)):
                sce_offsets_v = sce.GetPosOffsets( parentX, parentY , parentZ)
                parentSCEX = parentX - sce_offsets_v[0] + 0.6
                parentSCEY = parentY + sce_offsets_v[1]
                parentSCEZ = parentZ + sce_offsets_v[2]
                scedr = np.sqrt((parentSCEX-vtxX)**2 + (parentSCEY-vtxY)**2 + (parentSCEZ-vtxZ)**2)
            else:
                parentSCEX = parentX
                parentSCEY = parentY
                parentSCEZ = parentZ
                scedr = 99997

        # for 1mu1p (only difference is energy used)
        Ecal_1m1p              = ECal(mE,pE,'muon',B=BE)
        dEp_1m1p               = EpCCQE - Ecal_1m1p
        dEm_1m1p               = EmCCQE - Ecal_1m1p
        dEmp_1m1p              = EpCCQE-EmCCQE
        pTRat_1m1p             = pTransRat(mE,pE,lTh,pTh,lPh,pPh,'muon')
        Q2cal_1m1p             = Q2(Ecal_1m1p,mE,lTh,'muon')
        Q3_1m1p,Q0_1m1p        = Getq3q0(pE,mE,pTh,lTh,pPh,lPh,'muon',B=BE)
        EHad_1m1p              = (Ecal_1m1p - mE - 105.66)
        y_1m1p                 = EHad_1m1p/Ecal_1m1p
        x_1m1p                 = Q2cal_1m1p/(2*939.5654*EHad_1m1p)
        sph_1m1p               = sqrt(dEp_1m1p**2+dEm_1m1p**2+dEmp_1m1p**2)
        phiT_1m1p              = GetPhiT(mE,pE,lTh,pTh,lPh,pPh,'muon')
        pzenu_1m1p             = Getpz(mE,pE,lTh,pTh,'muon') - Ecal_1m1p
        pT_1m1p                = pTrans(mE,pE,lTh,pTh,lPh,pPh,'muon')
        alphT_1m1p             = alphaT(mE,pE,lTh,pTh,lPh,pPh,'muon')
        CCQE_energy_shift_1m1p = SensibleMinimize(mE,pE,lTh,pTh,lPh,pPh,'muon',B=BE)

        # Now boost these badbois
        try:
            pEB_1m1p,mEB_1m1p,pThB_1m1p,mThB_1m1p,pPhB_1m1p,mPhB_1m1p,EcalB_1m1p,EpCCQEB_1m1p,EmCCQEB_1m1p,sphB_1m1p = BoostTracks(mE,pE,lTh,pTh,lPh,pPh,'muon',B=BE)
            Q2calB_1m1p          = Q2(EcalB_1m1p,mEB_1m1p,mThB_1m1p)
            openAngB_1m1p        = OpenAngle(pThB_1m1p,mThB_1m1p,pPhB_1m1p,mPhB_1m1p)
            thetasB_1m1p         = mThB_1m1p+pThB_1m1p
            phisB_1m1p           = PhiDiff(mPhB_1m1p,pPhB_1m1p)
            EHadB_1m1p           = (EcalB_1m1p - mEB_1m1p - 105.66)
            yB_1m1p              = EHadB_1m1p/EcalB_1m1p
            xB_1m1p              = Q2calB_1m1p/(2*939.5654*EHadB_1m1p)
            phiTB_1m1p           = GetPhiT(mEB_1m1p,pEB_1m1p,mThB_1m1p,pThB_1m1p,mPhB_1m1p,pPhB_1m1p,'muon')
            pTB_1m1p             = pTrans(mEB_1m1p,pEB_1m1p,mThB_1m1p,pThB_1m1p,mPhB_1m1p,pPhB_1m1p,'muon')
            alphTB_1m1p          = alphaT(mEB_1m1p,pEB_1m1p,mThB_1m1p,pThB_1m1p,mPhB_1m1p,pPhB_1m1p,'muon')
            FailBoost = False
        except:
            FailBoost = True
            print 'BADBOOST: ',mE,pE,lTh,pTh,lPh,pPh

        #Now, let's hack in the electron stuff for now... this'll need to be changed later, but it'll suffice for the time being
        if PassShowerReco:
            # Second Shower Cut
            if ev.secondshower == 1:
                secsh_OpenAng = OpenAngle(ev.shr_theta,lTh,ev.shr_phi,lPh)
                if ev.shr_rad_pts > 15 and secsh_OpenAng < 0.866:
                    PassSecondShower = False

            Ecal_1e1p              = ECal(eE,pE,'electron',B=BE)
            dEp_1e1p               = EpCCQE - Ecal_1e1p
            dEe_1e1p               = EeCCQE - Ecal_1e1p
            dEep_1e1p              = EpCCQE-EeCCQE
            pTRat_1e1p             = pTransRat(eE,pE,lTh,pTh,lPh,pPh,'electron')
            Q2cal_1e1p             = Q2(Ecal_1e1p,eE,lTh,'electron')
            Q3_1e1p,Q0_1e1p      = Getq3q0(pE,eE,pTh,lTh,pPh,lPh,'electron',B=BE)
            EHad_1e1p              = (Ecal_1e1p - eE - .511)
            y_1e1p                 = EHad_1e1p/Ecal_1e1p
            x_1e1p                 = Q2cal_1e1p/(2*939.5654*EHad_1e1p)
            sph_1e1p               = sqrt(dEp_1e1p**2+dEe_1e1p**2+dEep_1e1p**2)
            phiT_1e1p              = GetPhiT(eE,pE,lTh,pTh,lPh,pPh,'electron')
            pT_1e1p                = pTrans(eE,pE,lTh,pTh,lPh,pPh,'electron')
            pzenu_1e1p             = Getpz(eE,pE,lTh,pTh,'electron') - Ecal_1e1p
            alphT_1e1p             = alphaT(eE,pE,lTh,pTh,lPh,pPh,'electron')
            CCQE_energy_shift_1e1p = SensibleMinimize(eE,pE,lTh,pTh,lPh,pPh,'electron',B=BE)

            # now booooost
            try:
                pEB_1e1p,eEB_1e1p,pThB_1e1p,eThB_1e1p,pPhB_1e1p,ePhB_1e1p,EcalB_1e1p,EpCCQEB_1e1p,EeCCQEB_1e1p,sphB_1e1p = BoostTracks(eE,pE,lTh,pTh,lPh,pPh,'electron',B=BE)
                Q2calB_1e1p          = Q2(EcalB_1e1p,eEB_1e1p,eThB_1e1p)
                openAngB_1e1p        = OpenAngle(pThB_1e1p,eThB_1e1p,pPhB_1e1p,ePhB_1e1p)
                thetasB_1e1p         = eThB_1e1p+pThB_1e1p
                phisB_1e1p           = PhiDiff(ePhB_1e1p,pPhB_1e1p)
                EHadB_1e1p           = (EcalB_1e1p - eEB_1e1p - .511)
                yB_1e1p              = EHadB_1e1p/EcalB_1e1p
                xB_1e1p              = Q2calB_1e1p/(2*939.5654*EHadB_1e1p)
                phiTB_1e1p           = GetPhiT(eEB_1e1p,pEB_1e1p,eThB_1e1p,pThB_1e1p,ePhB_1e1p,pPhB_1e1p,'electron')
                pTB_1e1p             = pTrans(eEB_1e1p,pEB_1e1p,eThB_1e1p,pThB_1e1p,ePhB_1e1p,pPhB_1e1p,'electron')
                alphTB_1e1p          = alphaT(eEB_1e1p,pEB_1e1p,eThB_1e1p,pThB_1e1p,ePhB_1e1p,pPhB_1e1p,'electron')
                FailBoost = False
            except:
                FailBoost = True
                print 'BADBOOST: ',eE,pE,lTh,pTh,lPh,pPh

    _run[0]          = run
    _subrun[0]       = subrun
    _event[0]        = event
    _vtxid[0]        = vtxid
    _x[0]            = vtxX
    _y[0]            = vtxY
    _z[0]            = vtxZ
    _anyReco[0]      = not(NothingRecod)
    _infiducial[0]   = InFiducial
    _ntracks[0]      = NumTracks
    _n5tracks[0]     = Num5cmTracks
    _passCuts[0]     = PassSimpleCuts
    _passShowerReco[0] = PassShowerReco
    _passSecShr[0]   = PassSecondShower
    _failBoost[0]    = FailBoost
    _good3DReco[0]   = Good3DReco
    _eta[0]          = eta                                      if PassSimpleCuts else -99999
    _openAng[0]      = openAng                                  if PassSimpleCuts else -99999
    _thetas[0]       = thetas                                   if PassSimpleCuts else -99999
    _phis[0]         = phis                                     if PassSimpleCuts else -99999
    _charge_near_trunk[0] = ldq + pdq                           if PassSimpleCuts else -99999
    _longtracklen[0]   = longtracklen                           if PassSimpleCuts else -99999
    _shorttracklen[0]  = shorttracklen                          if PassSimpleCuts else -99999
    _maxshrFrac[0]     = maxshrFrac                             if PassSimpleCuts else -99999
    _minshrFrac[0]     = minshrFrac                             if PassSimpleCuts else -99999

    _lepton_id[0] = int(lid)                                    if PassSimpleCuts else -99999
    _lepton_phi[0] = float(lPh)                                 if PassSimpleCuts else -99999
    _lepton_theta[0] = float(lTh)                               if PassSimpleCuts else -99999
    _lepton_length[0] = float(length_v[lid])                    if PassSimpleCuts else -99999
    _lepton_dqdx_calibrated[0] = float(dqdx_v_calibrated[lid])  if PassSimpleCuts else -99999
    _lepton_dqdx_uncalibrated[0] = float(dqdx_v_uncalibrated[lid])  if PassSimpleCuts else -99999
    _lepton_edge_dist[0] = float(wallLepton)                    if PassSimpleCuts else -99999
    _muon_E[0] = float(mE)                                      if PassSimpleCuts else -99999
    _electron_E[0] = float(eE)                                  if PassSimpleCuts else -99999
    _muon_phiB_1m1p[0] = float(mPhB_1m1p)                       if PassSimpleCuts and not FailBoost else -99999
    _muon_thetaB_1m1p[0] = float(mThB_1m1p)                     if PassSimpleCuts and not FailBoost else -99999
    _muon_EB_1m1p[0] = float(mEB_1m1p)                          if PassSimpleCuts and not FailBoost else -99999
    _electron_phiB_1e1p[0] = float(ePhB_1e1p)                   if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    _electron_thetaB_1e1p[0] = float(eThB_1e1p)                 if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    _electron_EB_1e1p[0] = float(eEB_1e1p)                      if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    _proton_id[0] = int(pid)                                    if PassSimpleCuts else -99999
    _proton_phi[0] = float(pPh)                                 if PassSimpleCuts else -99999
    _proton_theta[0] = float(pTh)                               if PassSimpleCuts else -99999
    _proton_length[0] = float(length_v[pid])                    if PassSimpleCuts else -99999
    _proton_dqdx_calibrated[0] = float(dqdx_v_calibrated[pid])  if PassSimpleCuts else -99999
    _proton_dqdx_uncalibrated[0] = float(dqdx_v_uncalibrated[pid])  if PassSimpleCuts else -99999
    _proton_edge_dist[0] = float(wallProton)                    if PassSimpleCuts else -99999
    _proton_E[0] = float(pE)                                    if PassSimpleCuts else -99999
    _proton_phiB_1m1p[0] = float(pPhB_1m1p)                     if PassSimpleCuts and not FailBoost else -99999
    _proton_thetaB_1m1p[0] = float(pThB_1m1p)                   if PassSimpleCuts and not FailBoost else -99999
    _proton_EB_1m1p[0] = float(pEB_1m1p)                        if PassSimpleCuts and not FailBoost else -99999
    _proton_phiB_1e1p[0] = float(pPhB_1e1p)                     if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    _proton_thetaB_1e1p[0] = float(pThB_1e1p)                   if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    _proton_EB_1e1p[0] = float(pEB_1e1p)                        if PassSimpleCuts and PassShowerReco and not FailBoost else -99999

    _CCQE_energy_shift_1m1p[0] = CCQE_energy_shift_1m1p         if PassSimpleCuts else -99999
    _enu_1m1p[0]          = Ecal_1m1p                           if PassSimpleCuts else -99999
    _phiT_1m1p[0]         = phiT_1m1p                           if PassSimpleCuts else -99999
    _alphaT_1m1p[0]       = alphT_1m1p                          if PassSimpleCuts else -99999
    _pT_1m1p[0]           = pT_1m1p                             if PassSimpleCuts else -99999
    _pTRat_1m1p[0]        = pTRat_1m1p                          if PassSimpleCuts else -99999
    _bjX_1m1p[0]          = x_1m1p                              if PassSimpleCuts else -99999
    _sph_1m1p[0]          = sph_1m1p                            if PassSimpleCuts else -99999
    _pzEnu_1m1p[0]        = pzenu_1m1p                          if PassSimpleCuts else -99999
    _q2_1m1p[0]           = Q2cal_1m1p                          if PassSimpleCuts else -99999
    _q0_1m1p[0]           = Q0_1m1p                             if PassSimpleCuts else -99999
    _q3_1m1p[0]           = Q3_1m1p                             if PassSimpleCuts else -99999
    _bjY_1m1p[0]          = y_1m1p                              if PassSimpleCuts else -99999
    _openAngB_1m1p[0]     = openAngB_1m1p                       if PassSimpleCuts and not FailBoost else -99995
    _thetasB_1m1p[0]      = thetasB_1m1p                        if PassSimpleCuts and not FailBoost else -99995
    _phisB_1m1p[0]        = phisB_1m1p                          if PassSimpleCuts and not FailBoost else -99995
    _phiTB_1m1p[0]        = phiTB_1m1p                          if PassSimpleCuts and not FailBoost else -99995
    _alphaTB_1m1p[0]      = alphTB_1m1p                         if PassSimpleCuts and not FailBoost else -99995
    _pTB_1m1p[0]          = pTB_1m1p                            if PassSimpleCuts and not FailBoost else -99995
    _bjXB_1m1p[0]         = xB_1m1p                             if PassSimpleCuts and not FailBoost else -99995
    _sphB_1m1p[0]         = sphB_1m1p                           if PassSimpleCuts and not FailBoost else -99995
    _q2B_1m1p[0]          = Q2calB_1m1p                         if PassSimpleCuts and not FailBoost else -99995
    _bjYB_1m1p[0]         = yB_1m1p                             if PassSimpleCuts and not FailBoost else -99995

    _CCQE_energy_shift_1e1p[0] = CCQE_energy_shift_1e1p         if PassSimpleCuts and PassShowerReco else -99999
    _enu_1e1p[0]          = Ecal_1e1p                           if PassSimpleCuts and PassShowerReco else -99999
    _phiT_1e1p[0]         = phiT_1e1p                           if PassSimpleCuts and PassShowerReco else -99999
    _alphaT_1e1p[0]       = alphT_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    _pT_1e1p[0]           = pT_1e1p                             if PassSimpleCuts and PassShowerReco else -99999
    _pTRat_1e1p[0]        = pTRat_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    _bjX_1e1p[0]          = x_1e1p                              if PassSimpleCuts and PassShowerReco else -99999
    _sph_1e1p[0]          = sph_1e1p                            if PassSimpleCuts and PassShowerReco else -99999
    _pzEnu_1e1p[0]        = pzenu_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    _q2_1e1p[0]           = Q2cal_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    _q0_1e1p[0]           = Q0_1e1p                             if PassSimpleCuts and PassShowerReco else -99999
    _q3_1e1p[0]           = Q3_1e1p                             if PassSimpleCuts and PassShowerReco else -99999
    _bjY_1e1p[0]          = y_1e1p                              if PassSimpleCuts and PassShowerReco else -99999
    _openAngB_1e1p[0]     = openAngB_1e1p                       if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _thetasB_1e1p[0]      = thetasB_1e1p                        if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _phisB_1e1p[0]        = phisB_1e1p                          if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _phiTB_1e1p[0]        = phiTB_1e1p                          if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _alphaTB_1e1p[0]      = alphTB_1e1p                         if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _pTB_1e1p[0]          = pTB_1e1p                            if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _bjXB_1e1p[0]         = xB_1e1p                             if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _sphB_1e1p[0]         = sphB_1e1p                           if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _q2B_1e1p[0]          = Q2calB_1e1p                         if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    _bjYB_1e1p[0]         = yB_1e1p                             if PassSimpleCuts and PassShowerReco and not FailBoost else -99995

    _totPE[0] = PMTPrecut_Dict[IDev]['_totpe']
    _porchTotPE[0] = PMTPrecut_Dict[IDev]['_porchtotpe']
    _maxPEFrac[0] = PMTPrecut_Dict[IDev]['_maxpefrac']
    _passPMTPrecut[0] = PMTPrecut_Dict[IDev]['_passpmtprecut']
    _parentPDG[0] = MC_dict[IDev]['parentPDG']                  if IsMC else -99998
    _energyInit[0] = MC_dict[IDev]['energyInit']                if IsMC else -99998
    _nproton[0] = MC_dict[IDev]['nproton']                      if IsMC else -99998
    _nlepton[0] = MC_dict[IDev]['nlepton']                      if IsMC else -99998
    _parentX[0] = parentX                                       if IsMC else -99998
    _parentY[0] = parentY                                       if IsMC else -99998
    _parentZ[0] = parentZ                                       if IsMC else -99998
    _parentSCEX[0] = parentSCEX                                 if IsMC and PassSimpleCuts else -99998
    _parentSCEY[0] = parentSCEY                                 if IsMC and PassSimpleCuts else -99998
    _parentSCEZ[0] = parentSCEZ                                 if IsMC and PassSimpleCuts else -99998
    _scedr[0] = scedr                                           if IsMC and PassSimpleCuts else -99998

    _eminusPID_int_v.clear()
    _muonPID_int_v.clear()
    _protonPID_int_v.clear()
    _gammaPID_int_v.clear()
    _pionPID_int_v.clear()
    _eminusPID_pix_v.clear()
    _muonPID_pix_v.clear()
    _protonPID_pix_v.clear()
    _gammaPID_pix_v.clear()
    _pionPID_pix_v.clear()
    for i in ev.eminus_int_score_torch:  _eminusPID_int_v.push_back(i)
    for i in ev.muon_int_score_torch:    _muonPID_int_v.push_back(i)
    for i in ev.proton_int_score_torch:  _protonPID_int_v.push_back(i)
    for i in ev.gamma_int_score_torch:   _gammaPID_int_v.push_back(i)
    for i in ev.pion_int_score_torch:    _pionPID_int_v.push_back(i)
    for i in ev.eminus_pix_score_torch:  _eminusPID_pix_v.push_back(i)
    for i in ev.muon_pix_score_torch:    _muonPID_pix_v.push_back(i)
    for i in ev.proton_pix_score_torch:  _protonPID_pix_v.push_back(i)
    for i in ev.gamma_pix_score_torch:   _gammaPID_pix_v.push_back(i)
    for i in ev.pion_pix_score_torch:    _pionPID_pix_v.push_back(i)

    outTree.Fill()
    clear_vertex()

DLMergedFile.Close()
MPIDFile.Close()

print()
print('<EVID: %s> -- Excellent! Now just write it up and we are done. Thanks for being patient'%_tag)
outTree.Write()
outFile.Close()
sys.exit(0)
