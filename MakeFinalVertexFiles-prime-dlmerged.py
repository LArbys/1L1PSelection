#v0: numu only so we can run on the dlmerged files without having to deal with shower energy for the moment
#v1: added pmt precut stuff
#v2: added 1e1p stuff back in with hack for shower energies


## # TODO:
# 1) various fixes for 1e1p (why use 2 tracks when we're looking for 1trk1sh -- could increase efficiency)

import ROOT
from ROOT import TFile,TTree
import matplotlib.pyplot as plt
import pickle
import numpy as np
from numpy import mean,asarray,matmul
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from sys import argv
from array import array
import os,sys

from SelectionDefs import NewAng, VtxInSimpleFid, VtxInFid, GetPhiT, pTrans,pTransRat, alphaT, ECCQE, ECal, Q2, OpenAngle, PhiDiff, edgeCut, ECCQE_mom, Boost, BoostTracks, Getpz, GetCCQEDiff, SensibleMinimize, Getq3q0

if len(argv) != 4:
    print('Fuck off')
    print('argv[1]: dlmerged.root')
    print('argv[2]: showerreco pickle')
    print('argv[3]: destination (.)')

_tag = argv[1][-27:-5]
_dest = argv[3]

# ---------------------------------------------- #

print()
print('<EVID: %s> -- First, we will figure out the PMT Precut info.'%_tag)  #gotta do this first for io reasons

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

    PMTfrac = [x / TotPE for x in PEbyPMT]
    if len(PMTfrac) == 0:
        MaxPEFrac = -1
    else:
        MaxPEFrac = max(PMTfrac)

    PassPMTPrecut = porchTotPE < PP_PE_THRESH and TotPE > PP_PE_THRESH and 0 < MaxPEFrac < PP_PMT_MAX_FRAC_CUTOFF

    PMTPrecut_Dict[id_rse] = {TotPE,porchTotPE,MaxPEFrac,PassPMTPrecut}

ll_manager.close()

# --- Open pickle file and take out the good stuff
with open(argv[2],"rb") as handle:
    [sh_theta_par1, sh_theta_par2, sh_phi_par1, sh_phi_par2, sh_energy_par1, sh_energy_par2] = pickle.load(handle,encoding="latin1")

# --- Open Ana File (hadd vertexana.root + tracker_anaout.root)
DLMergedFile = TFile(argv[1])

# --- Load relevant analysis trees from track and vertex ana files ----- #
try:
    TrkTree  = DLMergedFile.Get("_recoTree")
    vtxTree  = DLMergedFile.Get("VertexTree")
    shpTree  = DLMergedFile.Get("ShapeAnalysis")
    shrTree  = DLMergedFile.Get("SecondShowerAnalysis")

    TrkTree.AddFriend(vtxTree)
    TrkTree.AddFriend(shpTree)
    TrkTree.AddFriend(shrTree)
except:
    print "FUCKED: %s"%argv[1]
    sys.exit()

# ---------------------------------------------------------------------- #

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

# --- Create output ROOT file and initialize variables ----------------- #

outFileName = 'FinalVertexVariables-prime_'+_tag+'_NUMU.root'
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
_passCuts = MakeTreeBranch(outTree,'PassCuts','int')
_passSecShr = MakeTreeBranch(outTree,'PassSecShr','int')
_good3DReco = MakeTreeBranch(outTree,'Good3DReco','float')
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
_lepton_dqdx = MakeTreeBranch(outTree,'Lepton_dQdx','float')
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
_proton_dqdx = MakeTreeBranch(outTree,'Proton_dQdx','float')
_proton_E = MakeTreeBranch(outTree,'Proton_Edep','float')
_proton_edge_dist = MakeTreeBranch(outTree,'Proton_EdgeDist','float')
_proton_phiB_1m1p = MakeTreeBranch(outTree,'Proton_PhiRecoB_1m1p','float')
_proton_thetaB_1m1p = MakeTreeBranch(outTree,'Proton_ThetaRecoB_1m1p','float')
_proton_EB_1m1p = MakeTreeBranch(outTree,'Proton_EdepB_1m1p','float')
_proton_phiB_1e1p = MakeTreeBranch(outTree,'Proton_PhiRecoB_1e1p','float')
_proton_thetaB_1e1p = MakeTreeBranch(outTree,'Proton_ThetaRecoB_1e1p','float')
_proton_EB_1e1p = MakeTreeBranch(outTree,'Proton_EdepB_1e1p','float')

_phi_v = MakeTreeBranch(outTree,'phi_v','tvector')
_theta_v = MakeTreeBranch(outTree,'theta_v','tvector')
_length_v = MakeTreeBranch(outTree,'length_v','tvector')
_dqdx_v = MakeTreeBranch(outTree,'dqdx_v','tvector')
_EifP_v = MakeTreeBranch(outTree,'EifP_v','tvector')
_EifMu_v = MakeTreeBranch(outTree,'EifMu_v','tvector')

# Precut stuff
_totPE = MakeTreeBranch(outTree,'TotPE','float')
_porchTotPE = MakeTreeBranch(outTree,'PorchTotPE','float')
_maxPEFrac = MakeTreeBranch(outTree,'MaxPEFrac','float')
_passPMTPrecut = MakeTreeBranch(outTree,'PassPMTPrecut','int')

def clear_vertex():
    _lepton_id        = int(-9999)
    _lepton_phi       = float(-9999)
    _lepton_theta     = float(-9999)
    _lepton_length    = float(-9999)
    _lepton_dqdx      = float(-9999)
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
    _proton_dqdx      = float(-9999)
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

    passCuts = True
    passShowerReco = True
    passSecShr = True

    if not InFiducial: passCuts = False
    elif (NumTracks != 2 and Num5cmTracks != 2): passCuts = False
    elif edgeCut(EdgeDistance): passCuts = False

    vtxPhi_v       = ev.vertexPhi
    vtxTheta_v     = ev.vertexTheta
    dqdx_v         = ev.Avg_Ion_v
    iondlen_v      = ev.IondivLength_v
    Good3DReco     = ev.GoodVertex

    EifP_v         = ev.E_proton_v
    EifMu_v        = ev.E_muon_v

    # - Get ssnet and shape analysiss stuff
    sh_foundClusY = [i for i in ev.shower_frac_Y_v]
    sh_foundClusV = [i for i in ev.shower_frac_V_v]

    if 0.0 < max(sh_foundClusY) < 1.0:
        shrFrac    = max(foundClusY)
        shrFracPart = foundClus1
    else:
        shrFrac     = max(foundClusV)
        shrFracPart = foundClus2

    if NumTracks == 2:
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

        # Gotta find which particle is the proton in shower reco
        test_proton_openang_par1 = OpenAngle(pTh,sh_theta_par1[IDvtx],pPh,sh_phi_par1[IDvtx])
        test_proton_openang_par2 = OpenAngle(pTh,sh_theta_par2[IDvtx],pPh,sh_phi_par2[IDvtx])
        if test_proton_openang_par1 < test_proton_openang_par2:         # then par2 corresponds to the lepton
            try:
                eE = sh_energy_par2[IDvtx]
                if eE < 0: passShowerReco = False
            except:
                eE = -9999
                passShowerReco = False
        else:
            try:
                eE = sh_energy_par1[IDvtx]
                if eE < 0: passShowerReco = False
            except:
                eE = -9999
                passShowerReco = False

        thetas         = lTh+pTh
        phis           = PhiDiff(lPh,pPh)
        EpCCQE         = ECCQE(pE,pTh,pid="proton",B=BE)
        EmCCQE         = ECCQE(mE,lTh,pid="muon",B=BE)
        if(passShowerReco): EeCCQE = ECCQE(eE,lTh,pid="electron",B=BE)

        wallProton     = EdgeDistance[pid]
        wallLepton     = EdgeDistance[lid]

        openAng        = OpenAngle(pTh,lTh,pPh,lPh)
        eta            = abs(ev.Avg_Ion_v[pid]-ev.Avg_Ion_v[lid])/(ev.Avg_Ion_v[pid]+ev.Avg_Ion_v[lid])

        longtracklen   = max(ev.Length_v)
        shorttracklen  = min(ev.Length_v)
        maxshrFrac     = max(shrFracPart)
        minshrFrac     = min(shrFracPart)

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

        #Now, let's hack in the electron stuff for now... this'll need to be changed later, but it'll suffice for the time being
        if passShowerReco == 1:

            # Second Shower Cut
            if ev.secondshower == 1:
                secsh_OpenAng = OpenAngle(ev.shr_theta,lTh,ev.shr_phi,lPhi)
                if ev.shr_rad_pts > 15 and secsh_OpenAng < 0.866:
                    passSecShr = False

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
    _passCuts[0]     = passCuts
    _passShowerReco[0] = passShowerReco
    _passSecShr[0]   = passSecShr
    _good3DReco[0]   = Good3DReco
    _eta[0]          = eta       if NumTracks == 2 else -99999
    _openAng[0]      = openAng   if NumTracks == 2 else -99999
    _thetas[0]       = thetas    if NumTracks == 2 else -99999
    _phis[0]         = phis      if NumTracks == 2 else -99999
    _charge_near_trunk[0] = ldq + pdq if NumTracks == 2 else -99999
    _longtracklen[0]   = longtracklen if NumTracks == 2 else -99999
    _shorttracklen[0]  = shorttracklen if NumTracks == 2 else -99999
    _maxshrFrac[0]     = maxshrFrac if NumTracks==2 else -99999
    _minshrFrac[0]     = minshrFrac if NumTracks==2 else -99999

    _lepton_id[0] = int(lid) if NumTracks == 2 else -99999
    _lepton_phi[0] = float(vtxPhi_v[lid]) if NumTracks == 2 else -99999
    _lepton_theta[0] = float(vtxTheta_v[lid]) if NumTracks == 2 else -99999
    _lepton_length[0] = float(length_v[lid]) if NumTracks == 2 else -99999
    _lepton_dqdx[0] = float(dqdx_v[lid]) if NumTracks == 2 else -99999
    _lepton_edge_dist[0] = float(wallLepton) if NumTracks == 2 else -99999
    _muon_E[0] = float(EifMu_v.at(lid)) if NumTracks == 2 else -99999
    _electron_E[0] = float(eE) if NumTracks==2 else -99999
    _muon_phiB_1m1p[0] = float(mPhB_1m1p) if NumTracks == 2 else -99999
    _muon_thetaB_1m1p[0] = float(mThB_1m1p) if NumTracks == 2 else -99999
    _muon_EB_1m1p[0] = float(mEB_1m1p) if NumTracks == 2 else -99999
    _proton_id[0] = int(pid) if NumTracks == 2 else -99999
    _proton_phi[0] = float(vtxPhi_v[pid]) if NumTracks == 2 else -99999
    _proton_theta[0] = float(vtxTheta_v[pid]) if NumTracks == 2 else -99999
    _proton_length[0] = float(length_v[pid]) if NumTracks == 2 else -99999
    _proton_dqdx[0] = float(dqdx_v[pid]) if NumTracks == 2 else -99999
    _proton_edge_dist[0] = float(wallProton) if NumTracks == 2 else -99999
    _proton_E[0] = float(EifP_v.at(pid)) if NumTracks == 2 else -99999
    _proton_phiB_1m1p[0] = float(pPhB_1m1p) if NumTracks == 2 else -99999
    _proton_thetaB_1m1p[0] = float(pThB_1m1p) if NumTracks == 2 else -99999
    _proton_EB_1m1p[0] = float(pEB_1m1p) if NumTracks == 2 else -99999
    _proton_phiB_1e1p[0] = float(pPhB_1e1p) if (NumTracks==2 and passShowerReco) else -99999
    _proton_thetaB_1e1p[0] = float(pThB_1e1p) if (NumTracks==2 and passShowerReco) else -99999
    _proton_EB_1e1p[0] = float(pEB_1e1p) if (NumTracks==2 and passShowerReco) else -99999

    _CCQE_energy_shift_1m1p[0] = CCQE_energy_shift_1m1p if NumTracks == 2 else -99999
    _enu_1m1p[0]          = Ecal_1m1p      if NumTracks == 2 else -99999
    _phiT_1m1p[0]         = phiT_1m1p      if NumTracks == 2 else -99999
    _alphaT_1m1p[0]       = alphT_1m1p     if NumTracks == 2 else -99999
    _pT_1m1p[0]           = pT_1m1p        if NumTracks == 2 else -99999
    _pTRat_1m1p[0]        = pTRat_1m1p     if NumTracks == 2 else -99999
    _bjX_1m1p[0]          = x_1m1p         if NumTracks == 2 else -99999
    _sph_1m1p[0]          = sph_1m1p       if NumTracks == 2 else -99999
    _pzEnu_1m1p[0]        = pzenu_1m1p     if NumTracks == 2 else -99999
    _q2_1m1p[0]           = Q2cal_1m1p     if NumTracks == 2 else -99999
    _q0_1m1p[0]           = Q0_1m1p     if NumTracks == 2 else -99999
    _q3_1m1p[0]           = Q3_1m1p     if NumTracks == 2 else -99999
    _bjY_1m1p[0]          = y_1m1p         if NumTracks == 2 else -99999
    _openAngB_1m1p[0]     = openAngB_1m1p  if NumTracks == 2 else -99999
    _thetasB_1m1p[0]      = thetasB_1m1p   if NumTracks == 2 else -99999
    _phisB_1m1p[0]        = phisB_1m1p     if NumTracks == 2 else -99999
    _phiTB_1m1p[0]        = phiTB_1m1p     if NumTracks == 2 else -99999
    _alphaTB_1m1p[0]      = alphTB_1m1p    if NumTracks == 2 else -99999
    _pTB_1m1p[0]          = pTB_1m1p       if NumTracks == 2 else -99999
    _bjXB_1m1p[0]         = xB_1m1p        if NumTracks == 2 else -99999
    _sphB_1m1p[0]         = sphB_1m1p      if NumTracks == 2 else -99999
    _q2B_1m1p[0]          = Q2calB_1m1p    if NumTracks == 2 else -99999
    _bjYB_1m1p[0]         = yB_1m1p        if NumTracks == 2 else -99999

    _CCQE_energy_shift_1e1p[0] = CCQE_energy_shift_1e1p if (NumTracks==2 and passShowerReco) else -99999
    _enu_1e1p[0]          = Ecal_1e1p      if (NumTracks==2 and passShowerReco) else -99999
    _phiT_1e1p[0]         = phiT_1e1p      if (NumTracks==2 and passShowerReco) else -99999
    _alphaT_1e1p[0]       = alphT_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _pT_1e1p[0]           = pT_1e1p        if (NumTracks==2 and passShowerReco) else -99999
    _pTRat_1e1p[0]        = pTRat_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _bjX_1e1p[0]          = x_1e1p         if (NumTracks==2 and passShowerReco) else -99999
    _sph_1e1p[0]          = sph_1e1p       if (NumTracks==2 and passShowerReco) else -99999
    _pzEnu_1e1p[0]        = pzenu_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _q2_1e1p[0]           = Q2cal_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _q0_1e1p[0]           = Q0_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _q3_1e1p[0]           = Q3_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _bjY_1e1p[0]          = y_1e1p         if (NumTracks==2 and passShowerReco) else -99999
    _openAngB_1e1p[0]     = openAngB_1e1p  if (NumTracks==2 and passShowerReco) else -99999
    _thetasB_1e1p[0]      = thetasB_1e1p   if (NumTracks==2 and passShowerReco) else -99999
    _phisB_1e1p[0]        = phisB_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _phiTB_1e1p[0]        = phiTB_1e1p     if (NumTracks==2 and passShowerReco) else -99999
    _alphaTB_1e1p[0]      = alphTB_1e1p    if (NumTracks==2 and passShowerReco) else -99999
    _pTB_1e1p[0]          = pTB_1e1p       if (NumTracks==2 and passShowerReco) else -99999
    _bjXB_1e1p[0]         = xB_1e1p        if (NumTracks==2 and passShowerReco) else -99999
    _sphB_1e1p[0]         = sphB_1e1p      if (NumTracks==2 and passShowerReco) else -99999
    _q2B_1e1p[0]          = Q2calB_1e1p    if (NumTracks==2 and passShowerReco) else -99999
    _bjYB_1e1p[0]         = yB_1e1p        if (NumTracks==2 and passShowerReco) else -99999

    _totPE[0] = TotPE if NumTracks == 2 else -99999
    _porchTotPE[0] = porchTotPE if NumTracks == 2 else -99999
    _maxPEFrac[0] = MaxPEFrac if NumTracks == 2 else -99999
    _passPMTPrecut[0] = PassPMTPrecut if NumTracks == 2 else -99999

    _phi_v.clear()
    _theta_v.clear()
    _length_v.clear()
    _dqdx_v.clear()
    _EifP_v.clear()
    _EifMu_v.clear()
    for i in vtxPhi_v:      _phi_v.push_back(i)
    for i in vtxTheta_v:    _theta_v.push_back(i)
    for i in length_v:      _length_v.push_back(i)
    for i in dqdx_v:        _dqdx_v.push_back(i)
    for i in EifP_v:        _EifP_v.push_back(i)
    for i in EifMu_v:       _EifMu_v.push_back(i)

    outTree.Fill()
    clear_vertex()

print()
print('<EVID: %s> -- Excellent! Now just write it up and we are done. Thanks for being patient'%_tag)
outTree.Write()
outFile.Close()
sys.exit(0)
































































#print('Grab Info From Vtx Tree')
#for j,ev in enumerate(VtxTree):
#    print(str(j)+"/"+str(VtxTree.GetEntries())+"\r",end='')
#    idx = tuple((ev.run,ev.subrun,ev.event,ev.vtxid))
#    idxv = tuple((ev.run,ev.subrun,ev.event,ev.vtxid))

 #    if abs(ev.muon_int_score[2] - 0.24659) > 0.001:
#        muonPID[idxv] = ev.muon_int_score[2]
#    else:
#        muonPID[idxv] = ev.muon_int_score[1]
#
#    if abs(ev.proton_int_score[2] - 0.9999992) > 0.0001:
#        protPID[idxv] = ev.proton_int_score[2]
#    else:
#        protPID[idxv] = ev.proton_int_score[1]

#    if abs(ev.gamma_int_score[2] - 0.9506025) > 0.0001:
#        shwrPID[idxv] = max([ev.eminus_int_score[2],ev.gamma_int_score[2]])
#        gammaPID[idxv] = ev.gamma_int_score[2]
#        elecPID[idxv] = ev.eminus_int_score[2]
#    else:
#        shwrPID[idxv] = max([ev.eminus_int_score[1],ev.gamma_int_score[1]])
#        gammaPID[idxv] = ev.gamma_int_score[1]
#        elecPID[idxv] = ev.eminus_int_score[1]

#    secShr[idxv] = [ev.secondshower,ev.shr_min_dist,ev.shr_theta,ev.shr_phi,ev.shr_rad_pts]
#    VtxType[idxv] = ev.vertex_type

#    foundClus1 = [i for i in ev.shower_frac_Y_v]
#    foundClus2 = [i for i in ev.shower_frac_V_v]

#    if max(foundClus1) < 1.0 and max(foundClus1) > 0:
#        shrFrac[idxv]     = max(foundClus1)
#        shrFracPart[idxv] = foundClus1
#    else:
#        shrFrac[idxv]     = max(foundClus2)
#        shrFracPart[idxv] = foundClus2

#    foundClus = [i for i in ev.shower_frac_Y_v if i >=0]

#print()
#print('Okay. Now open up the trimmed pickle files')


#with open(argv[3],"rb") as handle:
#    [showerReco,showerRecoDir,vicsBDTScore] = pickle.load(handle,encoding="latin1")

print()
print('Great. Now loop through track events and go wild')

for indo,vtx in enumerate(TrkTree):
    print(str(indo)+'/'+str(TrkTree.GetEntries())+'\r',end='')

    BE = 29.5

    run            = vtx.run
    subrun         = vtx.subrun
    event          = vtx.event
    vtxid          = vtx.vtx_id
    IDvtx          = tuple((run,subrun,event,vtxid))
    vtxX           = vtx.vtx_x
    vtxY           = vtx.vtx_y
    vtxZ           = vtx.vtx_z

    # These are for passing simple precuts
    length_v            = vtx.Length_v
    InFiducialSimple    = VtxInSimpleFid(vtxX,vtxY,vtxZ)
    InFiducialTight     = VtxInFid(vtxX,vtxY,vtxZ)
    NumTracks           = len(length_v)
    Num5cmTracks        = sum(1 for x in length_v if x > 5)
    NothingRecod        = vtx.nothingReconstructed
    EdgeDistance        = vtx.closestWall

    PassCuts1m1p    = True
    PassCuts1e1p    = True
    PassShowerReco  = True
    PassSecShr      = True

    PassCuts1m1p = InFiducialSimple and (NumTracks == 2) and (Num5cmTracks == 2) and edgecut(EdgeDistance)

    NumShowers = vtx.nshowers

    PassCuts1e1p = InFiducialSimple and (NumShowers > 0 and NumShowers < 3) and (Num5cmTracks > 0 and Num5cmTracks <= 3) and edgecut(EdgeDistance)

    vtxPhi_v       = vtx.vertexPhi
    vtxTheta_v     = vtx.vertexTheta
    dqdx_v         = vtx.Avg_Ion_v
    iondlen_v      = vtx.IondivLength_v
    Good3DReco     = vtx.GoodVertex

    EifP_v         = vtx.E_proton_v
    EifMu_v        = vtx.E_muon_v

    if NumTracks == 2:
        mid            = int(np.argmin(vtx.Avg_Ion_v))
        pid            = int(np.argmax(vtx.Avg_Ion_v))
        mTh            = vtx.vertexTheta[mid]
        pTh            = vtx.vertexTheta[pid]
        mPh            = vtx.vertexPhi[mid]
        pPh            = vtx.vertexPhi[pid]
        # Add detector rotation fix
        mTh,mPh        = NewAng(mTh,mPh)
        pTh,pPh        = NewAng(pTh,pPh)
        # ----
        mE             = vtx.E_muon_v[mid]
        pE             = vtx.E_proton_v[pid]
        mdq = vtx.IonY_5cm_v[mid]
        pdq = vtx.IonY_5cm_v[pid]

        thetas         = mTh+pTh
        phis           = PhiDiff(mPh,pPh)
        EpCCQE         = ECCQE(pE,pTh,pid="proton",B=BE)
        EmCCQE         = ECCQE(mE,mTh,pid="muon",B=BE)

        wallProton     = EdgeDistance[pid]
        wallMuon     = EdgeDistance[mid]

        openAng        = OpenAngle(pTh,mTh,pPh,mPh)
        shwFrac        = shrFrac[IDvtx]
        eta            = abs(vtx.Avg_Ion_v[pid]-vtx.Avg_Ion_v[mid])/(vtx.Avg_Ion_v[pid]+vtx.Avg_Ion_v[mid])

        longtracklen   = max(vtx.Length_v)
        shorttracklen  = min(vtx.Length_v)

        # for 1mu1p (only difference is energy used)
        Ecal_1m1p              = ECal(mE,pE,'muon',B=BE)
        dEp_1m1p               = EpCCQE - Ecal_1m1p
        dEm_1m1p               = EmCCQE - Ecal_1m1p
        dEmp_1m1p              = EpCCQE-EmCCQE
        pTRat_1m1p             = pTransRat(mE,pE,mTh,pTh,mPh,pPh,'muon')
        Q2cal_1m1p             = Q2(Ecal_1m1p,mE,mTh,'muon')
        Q3_1m1p,Q0_1m1p      = Getq3q0(pE,mE,pTh,mTh,pPh,mPh,'muon',B=BE)
        EHad_1m1p              = (Ecal_1m1p - mE - 105.66)
        y_1m1p                 = EHad_1m1p/Ecal_1m1p
        x_1m1p                 = Q2cal_1m1p/(2*939.5654*EHad_1m1p)
        sph_1m1p               = sqrt(dEp_1m1p**2+dEm_1m1p**2+dEmp_1m1p**2)
        phiT_1m1p              = GetPhiT(mE,pE,mTh,pTh,mPh,pPh,'muon')
        pzenu_1m1p             = Getpz(mE,pE,mTh,pTh,'muon') - Ecal_1m1p
        pT_1m1p                = pTrans(mE,pE,mTh,pTh,mPh,pPh,'muon')
        alphT_1m1p             = alphaT(mE,pE,mTh,pTh,mPh,pPh,'muon')
        CCQE_energy_shift_1m1p = SensibleMinimize(mE,pE,mTh,pTh,mPh,pPh,'muon',B=BE)

        # Now boost these badbois
        pEB_1m1p,mEB_1m1p,pThB_1m1p,mThB_1m1p,pPhB_1m1p,mPhB_1m1p,EcalB_1m1p,EpCCQEB_1m1p,EmCCQEB_1m1p,sphB_1m1p = BoostTracks(mE,pE,mTh,pTh,mPh,pPh,'muon',B=BE)
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






if(passShowerReco): EeCCQE = ECCQE(eE,lTh,pid="electron",B=BE)
    maxshrFrac     = max(shrFracPart[IDvtx])
    minshrFrac     = min(shrFracPart[IDvtx])


#    if passShowerReco == 1:
#
#        if secShr[idxv][0] ==1:
#            shTheta = secShr[idxv][2]
#            shPhi   = secShr[idxv][3]
#            shOpen  = OpenAngle(shTheta,lTh,shPhi,lPh)
#            shDist  = secShr[idxv][1]
#            shSize  = secShr[idxv][4]
#
#            if shSize > 15 and shOpen < 0.866: passSecShr = False
#
#        try:
#            AngleDiffMu = acos(OpenAngle(lTh,acos(showerRecoDir[idxv][2]),lPh,atan2(showerRecoDir[idxv][1],showerRecoDir[idxv][0])))
#            AngleDiffP  = acos(OpenAngle(pTh,acos(showerRecoDir[idxv][2]),pPh,atan2(showerRecoDir[idxv][1],showerRecoDir[idxv][0])))
#        except:
#            AngleDiffMu = -5
#            AngleDiffP  = -5
#        if min(AngleDiffMu,AngleDiffP) > 0.5: passShowerReco = False
#
#        Ecal_1e1p              = ECal(eE,pE,'electron',B=BE)
#        dEp_1e1p               = EpCCQE - Ecal_1e1p
#        dEe_1e1p               = EeCCQE - Ecal_1e1p
#        dEep_1e1p              = EpCCQE-EeCCQE
#        pTRat_1e1p             = pTransRat(eE,pE,lTh,pTh,lPh,pPh,'electron')
#        Q2cal_1e1p             = Q2(Ecal_1e1p,eE,lTh,'electron')
#        Q3_1e1p,Q0_1e1p      = Getq3q0(pE,eE,pTh,lTh,pPh,lPh,'electron',B=BE)
#        EHad_1e1p              = (Ecal_1e1p - eE - .511)
#        y_1e1p                 = EHad_1e1p/Ecal_1e1p
#        x_1e1p                 = Q2cal_1e1p/(2*939.5654*EHad_1e1p)
#        sph_1e1p               = sqrt(dEp_1e1p**2+dEe_1e1p**2+dEep_1e1p**2)
#        phiT_1e1p              = GetPhiT(eE,pE,lTh,pTh,lPh,pPh,'electron')
#        pT_1e1p                = pTrans(eE,pE,lTh,pTh,lPh,pPh,'electron')
#        pzenu_1e1p             = Getpz(eE,pE,lTh,pTh,'electron') - Ecal_1e1p
#        alphT_1e1p             = alphaT(eE,pE,lTh,pTh,lPh,pPh,'electron')
#        CCQE_energy_shift_1e1p = SensibleMinimize(eE,pE,lTh,pTh,lPh,pPh,'electron',B=BE)
#
#        # now booost!
#        pEB_1e1p,eEB_1e1p,pThB_1e1p,eThB_1e1p,pPhB_1e1p,ePhB_1e1p,EcalB_1e1p,EpCCQEB_1e1p,EeCCQEB_1e1p,sphB_1e1p = BoostTracks(eE,pE,lTh,pTh,lPh,pPh,'electron',B=BE)
#        Q2calB_1e1p          = Q2(EcalB_1e1p,eEB_1e1p,eThB_1e1p)
#        openAngB_1e1p        = OpenAngle(pThB_1e1p,eThB_1e1p,pPhB_1e1p,ePhB_1e1p)
#        thetasB_1e1p         = eThB_1e1p+pThB_1e1p
#        phisB_1e1p           = PhiDiff(ePhB_1e1p,pPhB_1e1p)
#        EHadB_1e1p           = (EcalB_1e1p - eEB_1e1p - .511)
#        yB_1e1p              = EHadB_1e1p/EcalB_1e1p
#        xB_1e1p              = Q2calB_1e1p/(2*939.5654*EHadB_1e1p)
#        phiTB_1e1p           = GetPhiT(eEB_1e1p,pEB_1e1p,eThB_1e1p,pThB_1e1p,ePhB_1e1p,pPhB_1e1p,'electron')
#        pTB_1e1p             = pTrans(eEB_1e1p,pEB_1e1p,eThB_1e1p,pThB_1e1p,ePhB_1e1p,pPhB_1e1p,'electron')
#        alphTB_1e1p          = alphaT(eEB_1e1p,pEB_1e1p,eThB_1e1p,pThB_1e1p,ePhB_1e1p,pPhB_1e1p,'electron')

    _run[0]          = run
    _subrun[0]       = subrun
    _event[0]        = event
    _vtxid[0]        = vtxid
    _x[0]            = vtxX
    _y[0]            = vtxY
    _z[0]            = vtxZ
    _anyReco[0]      = not(NothingRecod)
    _infiducialsimple[0]   = InFiducialSimple
    _infiducialtight[0] = InFiducialTight
    _ntracks[0]      = NumTracks
    _n5tracks[0]     = Num5cmTracks

    _passCuts1e1p[0]     = PassCuts1e1p
    _passCuts1m1p[0]     = PassCuts1m1p
