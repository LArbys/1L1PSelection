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

if len(argv) != 3:
    print('Fuck off')
    print('argv[1]: dlmerged.root')
    print('argv[2]: destination (.)')

# --- Open Ana File (hadd vertexana.root + tracker_anaout.root)
DLMergedFile = TFile(argv[1])

# --- Load relevant analysis trees from track and vertex ana files ----- #
try:
    TrkTree  = DLMergedFile.Get("_recoTree")
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
_tag = argv[1][-27:-5]
_dest = argv[2]

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

_lepton_id = MakeTreeBranch(outTree,'Lepton_ID','int')
_lepton_phi = MakeTreeBranch(outTree,'Lepton_PhiReco','float')
_lepton_theta = MakeTreeBranch(outTree,'Lepton_ThetaReco','float')
_lepton_length = MakeTreeBranch(outTree,'Lepton_TrackLength','float')
_lepton_dqdx = MakeTreeBranch(outTree,'Lepton_dQdx','float')
_muon_E = MakeTreeBranch(outTree,'Muon_Edep','float')
_lepton_edge_dist = MakeTreeBranch(outTree,'Lepton_EdgeDist','float')
_muon_phiB_1m1p = MakeTreeBranch(outTree,'Muon_PhiRecoB_1m1p','float')
_muon_thetaB_1m1p = MakeTreeBranch(outTree,'Muon_ThetaRecoB_1m1p','float')
_muon_EB_1m1p = MakeTreeBranch(outTree,'Muon_EdepB_1m1p','float')
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

        thetas         = lTh+pTh
        phis           = PhiDiff(lPh,pPh)
        EpCCQE         = ECCQE(pE,pTh,pid="proton",B=BE)
        EmCCQE         = ECCQE(mE,lTh,pid="muon",B=BE)

        wallProton     = EdgeDistance[pid]
        wallLepton     = EdgeDistance[lid]

        openAng        = OpenAngle(pTh,lTh,pPh,lPh)
        eta            = abs(ev.Avg_Ion_v[pid]-ev.Avg_Ion_v[lid])/(ev.Avg_Ion_v[pid]+ev.Avg_Ion_v[lid])

        longtracklen   = max(ev.Length_v)
        shorttracklen  = min(ev.Length_v)

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
    _passSecShr[0]   = passSecShr
    _good3DReco[0]   = Good3DReco
    _eta[0]          = eta       if NumTracks == 2 else -99999
    _openAng[0]      = openAng   if NumTracks == 2 else -99999
    _thetas[0]       = thetas    if NumTracks == 2 else -99999
    _phis[0]         = phis      if NumTracks == 2 else -99999
    _charge_near_trunk[0] = ldq + pdq if NumTracks == 2 else -99999
    _longtracklen[0]   = longtracklen if NumTracks == 2 else -99999
    _shorttracklen[0]  = shorttracklen if NumTracks == 2 else -99999

    _lepton_id[0] = int(lid) if NumTracks == 2 else -99999
    _lepton_phi[0] = float(vtxPhi_v[lid]) if NumTracks == 2 else -99999
    _lepton_theta[0] = float(vtxTheta_v[lid]) if NumTracks == 2 else -99999
    _lepton_length[0] = float(length_v[lid]) if NumTracks == 2 else -99999
    _lepton_dqdx[0] = float(dqdx_v[lid]) if NumTracks == 2 else -99999
    _lepton_edge_dist[0] = float(wallLepton) if NumTracks == 2 else -99999
    _muon_E[0] = float(EifMu_v.at(lid)) if NumTracks == 2 else -99999
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
