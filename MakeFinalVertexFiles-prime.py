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

from SelectionDefs import NewAng, VtxInSimpleFid, VtxInFid, GetPhiT, pTrans,pTransRat, alphaT, ECCQE, ECal, Q2, OpenAngle, PhiDiff, edgeCut, ECCQE_mom, Boost, BoostTracks, Getpz, GetCCQEDiff, SensibleMinimize

if len(argv) != 6:
    print 'Fuck off'
    print  'argv[1]: dllee_vertex.root'
    print  'argv[2]: multipid_tree.root'
    print  'argv[3]: trimmed vertex pickle'
    print  'argv[4]: tag'
    print  'argv[5]: destination (.)'
    sys.exit()

# --- Open Ana File (hadd vertexana.root + tracker_anaout.root)
DLLEEVtxFile = TFile(argv[1])
MPIDRecoFile = TFile(argv[2])

# --- Load relevant analysis trees from track and vertex ana files ----- #
TrkTree  = DLLEEVtxFile.Get("_recoTree")
VtxTree  = DLLEEVtxFile.Get("VertexTree")
shpTree  = DLLEEVtxFile.Get("ShapeAnalysis")
shrTree  = DLLEEVtxFile.Get("SecondShowerAnalysis")
mpidTree = MPIDRecoFile.Get("multipid_tree")

VtxTree.AddFriend(shpTree)
VtxTree.AddFriend(shrTree)
VtxTree.AddFriend(mpidTree)
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
outFileName = 'FinalVertexVariables-prime'+argv[4]+'.root'
outFileName = os.path.join(sys.argv[5],outFileName)
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
_vtxalgo = MakeTreeBranch(outTree,'VtxAlgo','int')
_infiducial = MakeTreeBranch(outTree,'InFiducial','int')
_anyReco = MakeTreeBranch(outTree,'AnyReco','int')
_ntracks = MakeTreeBranch(outTree,'NTracks','int')
_n5tracks = MakeTreeBranch(outTree,'N5cmTracks','int')
_passCuts = MakeTreeBranch(outTree,'PassCuts','int')
_passShowerReco = MakeTreeBranch(outTree,'PassShowerReco','int')
_passSecShr = MakeTreeBranch(outTree,'PassSecShr','int')
_good3DReco = MakeTreeBranch(outTree,'Good3DReco','float')
_eta = MakeTreeBranch(outTree,'Eta','float')
_openAng = MakeTreeBranch(outTree,'OpenAng','float')
_thetas = MakeTreeBranch(outTree,'Thetas','float')
_phis = MakeTreeBranch(outTree,'Phis','float')
_charge_near_trunk = MakeTreeBranch(outTree,'ChargeNearTrunk','float')
_longtracklen = MakeTreeBranch(outTree,'LongTrackLen','float')
_shorttracklen = MakeTreeBranch(outTree,'ShortTrackLen','float')
_PIDmu = MakeTreeBranch(outTree,'PIDmu','float')
_PIDsh = MakeTreeBranch(outTree,'PIDsh','float')
_PIDp = MakeTreeBranch(outTree,'PIDp','float')
_PIDg = MakeTreeBranch(outTree,'PIDg','float')
_PIDe = MakeTreeBranch(outTree,'PIDe','float')
_maxshrFrac = MakeTreeBranch(outTree,'MaxShrFrac','float')
_minshrFrac = MakeTreeBranch(outTree,'MinShrFrac','float')
_BDTprob = MakeTreeBranch(outTree,'BDTprob','float')

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

shrFrac = {}
shrFracPart = {}
muonPID = {}
protPID = {}
shwrPID = {}
gammaPID = {}
elecPID = {}
secShr = {}
goodVtx = {}
VtxType = {}

print 'Grab Info From Vtx Tree'
for j,ev in enumerate(VtxTree):
    print str(j)+"/"+str(VtxTree.GetEntries())+"\r",
    idx = tuple((ev.run,ev.subrun,ev.event,ev.vtxid))
    idxv = tuple((ev.run,ev.subrun,ev.event,ev.vtxid))

    if abs(ev.muon_int_score[2] - 0.24659) > 0.001:
        muonPID[idxv] = ev.muon_int_score[2]
    else:
        muonPID[idxv] = ev.muon_int_score[1]

    if abs(ev.proton_int_score[2] - 0.9999992) > 0.0001:
        protPID[idxv] = ev.proton_int_score[2]
    else:
        protPID[idxv] = ev.proton_int_score[1]

    if abs(ev.gamma_int_score[2] - 0.9506025) > 0.0001:
        shwrPID[idxv] = max([ev.eminus_int_score[2],ev.gamma_int_score[2]])
        gammaPID[idxv] = ev.gamma_int_score[2]
        elecPID[idxv] = ev.eminus_int_score[2]
    else:
        shwrPID[idxv] = max([ev.eminus_int_score[1],ev.gamma_int_score[1]])
        gammaPID[idxv] = ev.gamma_int_score[1]
        elecPID[idxv] = ev.eminus_int_score[1]

    secShr[idxv] = [ev.secondshower,ev.shr_min_dist,ev.shr_theta,ev.shr_phi,ev.shr_rad_pts]
    VtxType[idxv] = ev.vertex_type

    foundClus1 = [i for i in ev.shower_frac_Y_v]
    foundClus2 = [i for i in ev.shower_frac_V_v]

    if max(foundClus1) < 1.0 and max(foundClus1) > 0:
        shrFrac[idxv]     = max(foundClus1)
        shrFracPart[idxv] = foundClus1
    else:
        shrFrac[idxv]     = max(foundClus2)
        shrFracPart[idxv] = foundClus2

    foundClus = [i for i in ev.shower_frac_Y_v if i >=0]

print
print 'Okay. Now open up the trimmed pickle files'
with open(argv[3],"rb") as handle: [showerReco,showerRecoDir,vicsBDTScore] = pickle.load(handle)

print
print 'Great. Now loop through track events and go wild'

for indo,ev in enumerate(TrkTree):
    print str(indo)+'/'+str(TrkTree.GetEntries())+'\r',

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
    InFiducial     = VtxInFid(vtxX,vtxY,vtxZ)
    NumTracks      = len(length_v)
    Num5cmTracks   = sum(1 for x in length_v if x > 5)
    NothingRecod   = ev.nothingReconstructed
    EdgeDistance   = ev.closestWall

    passCuts = True
    passShowerReco = True
    passSecShr = True

    if not InFiducial: passCuts = False
    elif (NumTracks != 2): passCuts = False
    elif edgeCut(EdgeDistance): passCuts = False
    try:
        eE = showerReco[IDvtx]
        if eE < 0: passShowerReco = False
    except:
        eE = -9999
        passShowerReco = False

    vtxPhi_v       = ev.vertexPhi
    vtxTheta_v     = ev.vertexTheta
    dqdx_v         = ev.Avg_Ion_v
    iondlen_v      = ev.IondivLength_v
    Good3DReco     = ev.GoodVertex

    VertexType     = VtxType[IDvtx]
    EifP_v         = ev.E_proton_v
    EifMu_v        = ev.E_muon_v

    if passCuts:
        lid            = np.argmin(ev.Avg_Ion_v)
        pid            = np.argmax(ev.Avg_Ion_v)
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
        EpCCQE         = ECCQE(pE,pTh,pid="proton")
        EmCCQE         = ECCQE(mE,lTh,pid="muon")
        if(passShowerReco): EeCCQE = ECCQE(eE,lTh,pid="electron")

        try:
            BDTprob    = -8888 if isnan(vicsBDTScore[IDvtx][0]) else vicsBDTScore[IDvtx][0]
        except:
            BDTprob    = -9999

        wallProton     = EdgeDistance[pid]
        wallLepton     = EdgeDistance[lid]

        openAng        = OpenAngle(pTh,lTh,pPh,lPh)
        shwFrac        = shrFrac[IDvtx]
        eta            = abs(ev.Avg_Ion_v[pid]-ev.Avg_Ion_v[lid])/(ev.Avg_Ion_v[pid]+ev.Avg_Ion_v[lid])

        longtracklen   = max(ev.Length_v)
        shorttracklen  = min(ev.Length_v)
        maxshrFrac     = max(shrFracPart[IDvtx])
        minshrFrac     = min(shrFracPart[IDvtx])

        # for 1mu1p (only difference is energy used)
        Ecal_1m1p              = ECal(mE,pE,'muon')
        dEp_1m1p               = EpCCQE - Ecal_1m1p
        dEm_1m1p               = EmCCQE - Ecal_1m1p
        dEmp_1m1p              = EpCCQE-EmCCQE
        pTRat_1m1p             = pTransRat(mE,pE,lTh,pTh,lPh,pPh,'muon')
        Q2cal_1m1p             = Q2(Ecal_1m1p,mE,lTh,'muon')
        EHad_1m1p              = (Ecal_1m1p - mE - 105.66)
        y_1m1p                 = EHad_1m1p/Ecal_1m1p
        x_1m1p                 = Q2cal_1m1p/(2*939.5654*EHad_1m1p)
        sph_1m1p               = sqrt(dEp_1m1p**2+dEm_1m1p**2+dEmp_1m1p**2)
        phiT_1m1p              = GetPhiT(mE,pE,lTh,pTh,lPh,pPh,'muon')
        pT_1m1p                = pTrans(mE,pE,lTh,pTh,lPh,pPh,'muon')
        pzenu_1m1p             = Getpz(mE,pE,lTh,pTh,'muon') - Ecal_1m1p
        alphT_1m1p             = alphaT(mE,pE,lTh,pTh,lPh,pPh,'muon')
        CCQE_energy_shift_1m1p = SensibleMinimize(mE,pE,lTh,pTh,lPh,pPh,'muon')

        # Now boost these badbois
        pEB_1m1p,mEB_1m1p,pThB_1m1p,mThB_1m1p,pPhB_1m1p,mPhB_1m1p,EcalB_1m1p,EpCCQEB_1m1p,EmCCQEB_1m1p,sphB_1m1p = BoostTracks(mE,pE,lTh,pTh,lPh,pPh,'muon')
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

        if passShowerReco == 1:

            if secShr[idxv][0] ==1:
                shTheta = secShr[idxv][2]
                shPhi   = secShr[idxv][3]
                shOpen  = OpenAngle(shTheta,lTh,shPhi,lPh)
                shDist  = secShr[idxv][1]
                shSize  = secShr[idxv][4]

                if shSize > 15 and shOpen < 0.866: passSecShr = False

            try:
                AngleDiffMu = acos(OpenAngle(lTh,acos(showerRecoDir[idxv][2]),lPh,atan2(showerRecoDir[idxv][1],showerRecoDir[idxv][0])))
                AngleDiffP  = acos(OpenAngle(pTh,acos(showerRecoDir[idxv][2]),pPh,atan2(showerRecoDir[idxv][1],showerRecoDir[idxv][0])))
            except:
                AngleDiffMu = -5
                AngleDiffP  = -5
            if min(AngleDiffMu,AngleDiffP) > 0.5: passShowerReco = False

            Ecal_1e1p              = ECal(eE,pE,'electron')
            dEp_1e1p               = EpCCQE - Ecal_1e1p
            dEe_1e1p               = EeCCQE - Ecal_1e1p
            dEep_1e1p              = EpCCQE-EeCCQE
            pTRat_1e1p             = pTransRat(eE,pE,lTh,pTh,lPh,pPh,'electron')
            Q2cal_1e1p             = Q2(Ecal_1e1p,eE,lTh,'electron')
            EHad_1e1p              = (Ecal_1e1p - eE - .511)
            y_1e1p                 = EHad_1e1p/Ecal_1e1p
            x_1e1p                 = Q2cal_1e1p/(2*939.5654*EHad_1e1p)
            sph_1e1p               = sqrt(dEp_1e1p**2+dEe_1e1p**2+dEep_1e1p**2)
            phiT_1e1p              = GetPhiT(eE,pE,lTh,pTh,lPh,pPh,'electron')
            pT_1e1p                = pTrans(eE,pE,lTh,pTh,lPh,pPh,'electron')
            pzenu_1e1p             = Getpz(eE,pE,lTh,pTh,'electron') - Ecal_1e1p
            alphT_1e1p             = alphaT(eE,pE,lTh,pTh,lPh,pPh,'electron')
            CCQE_energy_shift_1e1p = SensibleMinimize(eE,pE,lTh,pTh,lPh,pPh,'electron')

            # now booost!
            pEB_1e1p,eEB_1e1p,pThB_1e1p,eThB_1e1p,pPhB_1e1p,ePhB_1e1p,EcalB_1e1p,EpCCQEB_1e1p,EeCCQEB_1e1p,sphB_1e1p = BoostTracks(eE,pE,lTh,pTh,lPh,pPh,'electron')
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
    _vtxalgo[0]      = VertexType
    _anyReco[0]      = not(NothingRecod)
    _infiducial[0]   = InFiducial
    _ntracks[0]      = NumTracks
    _n5tracks[0]     = Num5cmTracks
    _passCuts[0]     = passCuts
    _passSecShr[0]   = passSecShr
    _passShowerReco[0] = passShowerReco
    _good3DReco[0]   = Good3DReco
    _eta[0]          = eta       if passCuts else -99999
    _openAng[0]      = openAng   if passCuts else -99999
    _thetas[0]       = thetas    if passCuts else -99999
    _phis[0]         = phis      if passCuts else -99999
    _charge_near_trunk[0] = ldq + pdq if passCuts else -99999
    _BDTprob[0]     = BDTprob if passCuts else -99999
    _longtracklen[0]   = longtracklen if passCuts else -99999
    _shorttracklen[0]  = shorttracklen if passCuts else -99999
    _PIDmu[0]          = muonPID[IDvtx] if passCuts else -99999
    _PIDsh[0]          = shwrPID[IDvtx] if passCuts else -99999
    _PIDp[0]           = protPID[IDvtx] if passCuts else -99999
    _PIDg[0]           = gammaPID[IDvtx] if passCuts else -99999
    _PIDe[0]           = elecPID[IDvtx] if passCuts else -99999
    _maxshrFrac[0]     = maxshrFrac if passCuts else -99999
    _minshrFrac[0]     = minshrFrac if passCuts else -99999

    _lepton_id[0] = int(lid) if passCuts else -99999
    _lepton_phi[0] = float(vtxPhi_v[lid]) if passCuts else -99999
    _lepton_theta[0] = float(vtxTheta_v[lid]) if passCuts else -99999
    _lepton_length[0] = float(length_v[lid]) if passCuts else -99999
    _lepton_dqdx[0] = float(dqdx_v[lid]) if passCuts else -99999
    _lepton_edge_dist[0] = float(wallLepton) if passCuts else -99999
    _muon_E[0] = float(EifMu_v.at(lid)) if passCuts else -99999
    _electron_E[0] = float(eE) if passCuts else -99999
    _muon_phiB_1m1p[0] = float(mPhB_1m1p) if passCuts else -99999
    _muon_thetaB_1m1p[0] = float(mThB_1m1p) if passCuts else -99999
    _muon_EB_1m1p[0] = float(mEB_1m1p) if passCuts else -99999
    _electron_phiB_1e1p[0] = float(ePhB_1e1p) if (passCuts and passShowerReco) else -99999
    _electron_thetaB_1e1p[0] = float(eThB_1e1p) if (passCuts and passShowerReco) else -99999
    _electron_EB_1e1p[0] = float(eEB_1e1p) if (passCuts and passShowerReco) else -99999
    _proton_id[0] = int(pid) if passCuts else -99999
    _proton_phi[0] = float(vtxPhi_v[pid]) if passCuts else -99999
    _proton_theta[0] = float(vtxTheta_v[pid]) if passCuts else -99999
    _proton_length[0] = float(length_v[pid]) if passCuts else -99999
    _proton_dqdx[0] = float(dqdx_v[pid]) if passCuts else -99999
    _proton_edge_dist[0] = float(wallProton) if passCuts else -99999
    _proton_E[0] = float(EifP_v.at(pid)) if passCuts else -99999
    _proton_phiB_1m1p[0] = float(pPhB_1m1p) if passCuts else -99999
    _proton_thetaB_1m1p[0] = float(pThB_1m1p) if passCuts else -99999
    _proton_EB_1m1p[0] = float(pEB_1m1p) if passCuts else -99999
    _proton_phiB_1e1p[0] = float(pPhB_1e1p) if (passCuts and passShowerReco) else -99999
    _proton_thetaB_1e1p[0] = float(pThB_1e1p) if (passCuts and passShowerReco) else -99999
    _proton_EB_1e1p[0] = float(pEB_1e1p) if (passCuts and passShowerReco) else -99999

    _CCQE_energy_shift_1m1p[0] = CCQE_energy_shift_1m1p if passCuts else -99999
    _enu_1m1p[0]          = Ecal_1m1p      if passCuts else -99999
    _phiT_1m1p[0]         = phiT_1m1p      if passCuts else -99999
    _alphaT_1m1p[0]       = alphT_1m1p     if passCuts else -99999
    _pT_1m1p[0]           = pT_1m1p        if passCuts else -99999
    _pTRat_1m1p[0]        = pTRat_1m1p     if passCuts else -99999
    _bjX_1m1p[0]          = x_1m1p         if passCuts else -99999
    _sph_1m1p[0]          = sph_1m1p       if passCuts else -99999
    _pzEnu_1m1p[0]        = pzenu_1m1p     if passCuts else -99999
    _q2_1m1p[0]           = Q2cal_1m1p     if passCuts else -99999
    _bjY_1m1p[0]          = y_1m1p         if passCuts else -99999
    _openAngB_1m1p[0]     = openAngB_1m1p  if passCuts else -99999
    _thetasB_1m1p[0]      = thetasB_1m1p   if passCuts else -99999
    _phisB_1m1p[0]        = phisB_1m1p     if passCuts else -99999
    _phiTB_1m1p[0]        = phiTB_1m1p     if passCuts else -99999
    _alphaTB_1m1p[0]      = alphTB_1m1p    if passCuts else -99999
    _pTB_1m1p[0]          = pTB_1m1p       if passCuts else -99999
    _bjXB_1m1p[0]         = xB_1m1p        if passCuts else -99999
    _sphB_1m1p[0]         = sphB_1m1p      if passCuts else -99999
    _q2B_1m1p[0]          = Q2calB_1m1p    if passCuts else -99999
    _bjYB_1m1p[0]         = yB_1m1p        if passCuts else -99999

    _CCQE_energy_shift_1e1p[0] = CCQE_energy_shift_1e1p if (passCuts and passShowerReco) else -99999
    _enu_1e1p[0]          = Ecal_1e1p      if (passCuts and passShowerReco) else -99999
    _phiT_1e1p[0]         = phiT_1e1p      if (passCuts and passShowerReco) else -99999
    _alphaT_1e1p[0]       = alphT_1e1p     if (passCuts and passShowerReco) else -99999
    _pT_1e1p[0]           = pT_1e1p        if (passCuts and passShowerReco) else -99999
    _pTRat_1e1p[0]        = pTRat_1e1p     if (passCuts and passShowerReco) else -99999
    _bjX_1e1p[0]          = x_1e1p         if (passCuts and passShowerReco) else -99999
    _sph_1e1p[0]          = sph_1e1p       if (passCuts and passShowerReco) else -99999
    _pzEnu_1e1p[0]        = pzenu_1e1p     if (passCuts and passShowerReco) else -99999
    _q2_1e1p[0]           = Q2cal_1e1p     if (passCuts and passShowerReco) else -99999
    _bjY_1e1p[0]          = y_1e1p         if (passCuts and passShowerReco) else -99999
    _openAngB_1e1p[0]     = openAngB_1e1p  if (passCuts and passShowerReco) else -99999
    _thetasB_1e1p[0]      = thetasB_1e1p   if (passCuts and passShowerReco) else -99999
    _phisB_1e1p[0]        = phisB_1e1p     if (passCuts and passShowerReco) else -99999
    _phiTB_1e1p[0]        = phiTB_1e1p     if (passCuts and passShowerReco) else -99999
    _alphaTB_1e1p[0]      = alphTB_1e1p    if (passCuts and passShowerReco) else -99999
    _pTB_1e1p[0]          = pTB_1e1p       if (passCuts and passShowerReco) else -99999
    _bjXB_1e1p[0]         = xB_1e1p        if (passCuts and passShowerReco) else -99999
    _sphB_1e1p[0]         = sphB_1e1p      if (passCuts and passShowerReco) else -99999
    _q2B_1e1p[0]          = Q2calB_1e1p    if (passCuts and passShowerReco) else -99999
    _bjYB_1e1p[0]         = yB_1e1p        if (passCuts and passShowerReco) else -99999

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

print
print 'Excellent! Now just write it up and we are done. Thanks for being patient'
outTree.Write()
outFile.Close()
sys.exit(0)
