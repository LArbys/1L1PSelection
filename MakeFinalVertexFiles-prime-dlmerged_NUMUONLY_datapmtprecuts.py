#v1: added pmt precut stuff
#v2: added mctruth info (where applicable, obvi)
#v3: added energy calibration to dqdx
#v4: added mpidtree
#v5: added new shower reco from nick

# TODO:
# - make sure CCQE energies match proper off-shell formulae  (talk to steven about that)

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
from LEEPreCuts_Functions import GetPMTPrecutDict, PerformPMTPrecuts

from SelectionDefs import NewAng, VtxInSimpleFid, VtxInFid, GetPhiT, pTrans,pTransRat, alphaT, ECCQE, ECal, Q2, OpenAngle, PhiDiff, edgeCut, ECCQE_mom, Boost, BoostTracks, Getpz, GetCCQEDiff, SensibleMinimize, Getq3q0,GetTotPE, CorrectionFactor

if len(argv) != 4:
    print('Fuck off')
    print('argv[1]: dlmerged.root')
    print('argv[2]: calibration map')
    print('argv[3]: destination (.)')
    sys.exit()

_tag = argv[1][-27:-5]
_dest = argv[3]
_dlmerged = argv[1]
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



# -----------------------------------------------------------------------#
#    Get started!
# -----------------------------------------------------------------------#

print()
print('<EVID: %s> -- First, we will figure out the PMT Precut info.'%_tag)  #gotta do this first for io reasons
# PMTPrecut_Dict = GetPMTPrecutDict(_dlmerged)
PMTPrecut_Dict = PerformPMTPrecuts(_dlmerged)

print()
print('<EVID: %s> -- Now make sure we can read the root trees we want (and make sure they are present).'%_tag)
try:
    DLMergedFile = TFile(_dlmerged,'read')

    TrkTree  = DLMergedFile.Get("_recoTree")
    VtxTree  = DLMergedFile.Get("VertexTree")
    ShpTree  = DLMergedFile.Get("ShapeAnalysis")
    ShrTree  = DLMergedFile.Get("SecondShowerAnalysis")
    MCTree   = DLMergedFile.Get("MCTree")

    TrkTree.AddFriend(VtxTree)
    TrkTree.AddFriend(ShpTree)
    TrkTree.AddFriend(ShrTree)
except:
    print 'Fucked: %s' %(_dlmerged)
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
outFileName = 'FinalVertexVariables-prime_'+_tag+'_1M1P.root'
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

_lepton_id = MakeTreeBranch(outTree,'Lepton_ID','int')
_lepton_phi = MakeTreeBranch(outTree,'Lepton_PhiReco','float')
_lepton_theta = MakeTreeBranch(outTree,'Lepton_ThetaReco','float')
_lepton_length = MakeTreeBranch(outTree,'Lepton_TrackLength','float')
_lepton_dqdx_uncalibrated = MakeTreeBranch(outTree,'Lepton_dQdx_uncalibrated','float')
_lepton_dqdx_calibrated = MakeTreeBranch(outTree,'Lepton_dQdx','float')
_muon_E = MakeTreeBranch(outTree,'Muon_Edep','float')
_lepton_edge_dist = MakeTreeBranch(outTree,'Lepton_EdgeDist','float')
_muon_phiB_1m1p = MakeTreeBranch(outTree,'Muon_PhiRecoB_1m1p','float')
_muon_thetaB_1m1p = MakeTreeBranch(outTree,'Muon_ThetaRecoB_1m1p','float')
_muon_EB_1m1p = MakeTreeBranch(outTree,'Muon_EdepB_1m1p','float')
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
    _muon_phiB_1m1p[0] = float(mPhB_1m1p)                       if PassSimpleCuts and not FailBoost else -99999
    _muon_thetaB_1m1p[0] = float(mThB_1m1p)                     if PassSimpleCuts and not FailBoost else -99999
    _muon_EB_1m1p[0] = float(mEB_1m1p)                          if PassSimpleCuts and not FailBoost else -99999
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

    outTree.Fill()
    clear_vertex()

DLMergedFile.Close()

print()
print('<EVID: %s> -- Excellent! Now just write it up and we are done. Thanks for being patient'%_tag)
outTree.Write()
outFile.Close()
sys.exit(0)
