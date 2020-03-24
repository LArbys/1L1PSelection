import numpy as np
from math import sqrt
from SelectionDefs import NewAng, VtxInSimpleFid, VtxInFid, GetPhiT, pTrans,pTransRat
from SelectionDefs import alphaT, ECCQE, ECal, Q2, OpenAngle, PhiDiff, edgeCut
from SelectionDefs import ECCQE_mom, Boost, BoostTracks, Getpz, GetCCQEDiff, SensibleMinimize, Getq3q0,GetTotPE, CorrectionFactor

def make_selection_vars( indo, ismc,
                         ev, df_ShowerReco, PMTPrecut_Dict, MC_dict,
                         dlanavars, calibMap_v,
                         sce=None ):
    """
    inputs
    -------
    indo [int] entry number
    ev [TTree] merged event vertex tree

    """

    if sce is None:
        sce = larutil.SpaceChargeMicroBooNEMCC9()

    # Hard-coded binding energy!!
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

    vtxShowerData  = df_ShowerReco.get_vertex( run, subrun, event, vtxid )

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
        eE  = vtxShowerData['shower_energies'][2]
        eE2 = vtxShowerData['secondshower_energies'][2]
        if eE < 0:
            eE = -99999
            PassShowerReco = False
    except:
        eE  = -9999
        eE2 = -9999
        PassShowerReco = False

    # - get MC Truth info (if applicable)
    if ismc:
        parentX = MC_dict[IDev]['parentX']
        parentY = MC_dict[IDev]['parentY']
        parentZ = MC_dict[IDev]['parentZ']

    # - Now, the big stuff
    lid = -1
    pid = -1
    if PassSimpleCuts:
        lid            = int(np.argmin(ev.Avg_Ion_v))
        pid            = int(np.argmax(ev.Avg_Ion_v))
        print "[passed simple cuts] lid=",lid," pid=",pid
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

        if ismc:
            for i in xrange(0,len(dqdx_v_uncalibrated)):
                dqdx_v_calibrated[i] = dqdx_v_uncalibrated[i] * CorrectionFactor(vtxX,vtxY,vtxZ,vtxTheta_v[i],vtxPhi_v[i],length_v[i],calibMap_v)

        eta = abs(dqdx_v_calibrated[pid]-dqdx_v_calibrated[lid])/(dqdx_v_calibrated[pid]+dqdx_v_calibrated[lid])

        longtracklen   = max(ev.Length_v)
        shorttracklen  = min(ev.Length_v)
        maxshrFrac     = max(shrFracPart)
        minshrFrac     = min(shrFracPart)

        # - Correct SCE of coords
        if(ismc):
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

    dlanavars._run[0]          = run
    dlanavars._subrun[0]       = subrun
    dlanavars._event[0]        = event
    dlanavars._vtxid[0]        = vtxid
    dlanavars._x[0]            = vtxX
    dlanavars._y[0]            = vtxY
    dlanavars._z[0]            = vtxZ
    dlanavars._anyReco[0]      = not(NothingRecod)
    dlanavars._infiducial[0]   = InFiducial
    dlanavars._ntracks[0]      = NumTracks
    dlanavars._n5tracks[0]     = Num5cmTracks
    dlanavars._passCuts[0]     = PassSimpleCuts
    dlanavars._passShowerReco[0] = PassShowerReco
    dlanavars._passSecShr[0]   = PassSecondShower
    dlanavars._failBoost[0]    = FailBoost
    dlanavars._good3DReco[0]   = Good3DReco
    dlanavars._eta[0]          = eta                                      if PassSimpleCuts else -99999
    dlanavars._openAng[0]      = openAng                                  if PassSimpleCuts else -99999
    dlanavars._thetas[0]       = thetas                                   if PassSimpleCuts else -99999
    dlanavars._phis[0]         = phis                                     if PassSimpleCuts else -99999
    dlanavars._charge_near_trunk[0] = ldq + pdq                           if PassSimpleCuts else -99999
    dlanavars._longtracklen[0]   = longtracklen                           if PassSimpleCuts else -99999
    dlanavars._shorttracklen[0]  = shorttracklen                          if PassSimpleCuts else -99999
    dlanavars._maxshrFrac[0]     = maxshrFrac                             if PassSimpleCuts else -99999
    dlanavars._minshrFrac[0]     = minshrFrac                             if PassSimpleCuts else -99999

    dlanavars._lepton_id[0] = int(lid)                                    if PassSimpleCuts else -99999
    dlanavars._lepton_phi[0] = float(lPh)                                 if PassSimpleCuts else -99999
    dlanavars._lepton_theta[0] = float(lTh)                               if PassSimpleCuts else -99999
    dlanavars._lepton_length[0] = float(length_v[lid])                    if PassSimpleCuts else -99999
    dlanavars._lepton_dqdx_calibrated[0] = float(dqdx_v_calibrated[lid])  if PassSimpleCuts else -99999
    dlanavars._lepton_dqdx_uncalibrated[0] = float(dqdx_v_uncalibrated[lid])  if PassSimpleCuts else -99999
    dlanavars._lepton_edge_dist[0] = float(wallLepton)                    if PassSimpleCuts else -99999
    dlanavars._muon_E[0] = float(mE)                                      if PassSimpleCuts else -99999
    dlanavars._electron_E[0] = float(eE)                                  if PassSimpleCuts else -99999
    dlanavars._muon_phiB_1m1p[0] = float(mPhB_1m1p)                       if PassSimpleCuts and not FailBoost else -99999
    dlanavars._muon_thetaB_1m1p[0] = float(mThB_1m1p)                     if PassSimpleCuts and not FailBoost else -99999
    dlanavars._muon_EB_1m1p[0] = float(mEB_1m1p)                          if PassSimpleCuts and not FailBoost else -99999
    dlanavars._electron_phiB_1e1p[0] = float(ePhB_1e1p)                   if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    dlanavars._electron_thetaB_1e1p[0] = float(eThB_1e1p)                 if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    dlanavars._electron_EB_1e1p[0] = float(eEB_1e1p)                      if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    dlanavars._proton_id[0] = int(pid)                                    if PassSimpleCuts else -99999
    dlanavars._proton_phi[0] = float(pPh)                                 if PassSimpleCuts else -99999
    dlanavars._proton_theta[0] = float(pTh)                               if PassSimpleCuts else -99999
    dlanavars._proton_length[0] = float(length_v[pid])                    if PassSimpleCuts else -99999
    dlanavars._proton_dqdx_calibrated[0] = float(dqdx_v_calibrated[pid])  if PassSimpleCuts else -99999
    dlanavars._proton_dqdx_uncalibrated[0] = float(dqdx_v_uncalibrated[pid])  if PassSimpleCuts else -99999
    dlanavars._proton_edge_dist[0] = float(wallProton)                    if PassSimpleCuts else -99999
    dlanavars._proton_E[0] = float(pE)                                    if PassSimpleCuts else -99999
    dlanavars._proton_phiB_1m1p[0] = float(pPhB_1m1p)                     if PassSimpleCuts and not FailBoost else -99999
    dlanavars._proton_thetaB_1m1p[0] = float(pThB_1m1p)                   if PassSimpleCuts and not FailBoost else -99999
    dlanavars._proton_EB_1m1p[0] = float(pEB_1m1p)                        if PassSimpleCuts and not FailBoost else -99999
    dlanavars._proton_phiB_1e1p[0] = float(pPhB_1e1p)                     if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    dlanavars._proton_thetaB_1e1p[0] = float(pThB_1e1p)                   if PassSimpleCuts and PassShowerReco and not FailBoost else -99999
    dlanavars._proton_EB_1e1p[0] = float(pEB_1e1p)                        if PassSimpleCuts and PassShowerReco and not FailBoost else -99999

    dlanavars._CCQE_energy_shift_1m1p[0] = CCQE_energy_shift_1m1p         if PassSimpleCuts else -99999
    dlanavars._enu_1m1p[0]          = Ecal_1m1p                           if PassSimpleCuts else -99999
    dlanavars._phiT_1m1p[0]         = phiT_1m1p                           if PassSimpleCuts else -99999
    dlanavars._alphaT_1m1p[0]       = alphT_1m1p                          if PassSimpleCuts else -99999
    dlanavars._pT_1m1p[0]           = pT_1m1p                             if PassSimpleCuts else -99999
    dlanavars._pTRat_1m1p[0]        = pTRat_1m1p                          if PassSimpleCuts else -99999
    dlanavars._bjX_1m1p[0]          = x_1m1p                              if PassSimpleCuts else -99999
    dlanavars._sph_1m1p[0]          = sph_1m1p                            if PassSimpleCuts else -99999
    dlanavars._pzEnu_1m1p[0]        = pzenu_1m1p                          if PassSimpleCuts else -99999
    dlanavars._q2_1m1p[0]           = Q2cal_1m1p                          if PassSimpleCuts else -99999
    dlanavars._q0_1m1p[0]           = Q0_1m1p                             if PassSimpleCuts else -99999
    dlanavars._q3_1m1p[0]           = Q3_1m1p                             if PassSimpleCuts else -99999
    dlanavars._bjY_1m1p[0]          = y_1m1p                              if PassSimpleCuts else -99999
    dlanavars._openAngB_1m1p[0]     = openAngB_1m1p                       if PassSimpleCuts and not FailBoost else -99995
    dlanavars._thetasB_1m1p[0]      = thetasB_1m1p                        if PassSimpleCuts and not FailBoost else -99995
    dlanavars._phisB_1m1p[0]        = phisB_1m1p                          if PassSimpleCuts and not FailBoost else -99995
    dlanavars._phiTB_1m1p[0]        = phiTB_1m1p                          if PassSimpleCuts and not FailBoost else -99995
    dlanavars._alphaTB_1m1p[0]      = alphTB_1m1p                         if PassSimpleCuts and not FailBoost else -99995
    dlanavars._pTB_1m1p[0]          = pTB_1m1p                            if PassSimpleCuts and not FailBoost else -99995
    dlanavars._bjXB_1m1p[0]         = xB_1m1p                             if PassSimpleCuts and not FailBoost else -99995
    dlanavars._sphB_1m1p[0]         = sphB_1m1p                           if PassSimpleCuts and not FailBoost else -99995
    dlanavars._q2B_1m1p[0]          = Q2calB_1m1p                         if PassSimpleCuts and not FailBoost else -99995
    dlanavars._bjYB_1m1p[0]         = yB_1m1p                             if PassSimpleCuts and not FailBoost else -99995

    dlanavars._CCQE_energy_shift_1e1p[0] = CCQE_energy_shift_1e1p         if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._enu_1e1p[0]          = Ecal_1e1p                           if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._phiT_1e1p[0]         = phiT_1e1p                           if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._alphaT_1e1p[0]       = alphT_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._pT_1e1p[0]           = pT_1e1p                             if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._pTRat_1e1p[0]        = pTRat_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._bjX_1e1p[0]          = x_1e1p                              if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._sph_1e1p[0]          = sph_1e1p                            if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._pzEnu_1e1p[0]        = pzenu_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._q2_1e1p[0]           = Q2cal_1e1p                          if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._q0_1e1p[0]           = Q0_1e1p                             if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._q3_1e1p[0]           = Q3_1e1p                             if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._bjY_1e1p[0]          = y_1e1p                              if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._openAngB_1e1p[0]     = openAngB_1e1p                       if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._thetasB_1e1p[0]      = thetasB_1e1p                        if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._phisB_1e1p[0]        = phisB_1e1p                          if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._phiTB_1e1p[0]        = phiTB_1e1p                          if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._alphaTB_1e1p[0]      = alphTB_1e1p                         if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._pTB_1e1p[0]          = pTB_1e1p                            if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._bjXB_1e1p[0]         = xB_1e1p                             if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._sphB_1e1p[0]         = sphB_1e1p                           if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._q2B_1e1p[0]          = Q2calB_1e1p                         if PassSimpleCuts and PassShowerReco and not FailBoost else -99995
    dlanavars._bjYB_1e1p[0]         = yB_1e1p                             if PassSimpleCuts and PassShowerReco and not FailBoost else -99995

    # shower reco variables
    dlanavars._shower1_E_U[0]       = float(vtxShowerData["shower_energies"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_E_V[0]       = float(vtxShowerData["shower_energies"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_E_Y[0]       = float(vtxShowerData["shower_energies"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_E_U[0]       = float(vtxShowerData["secondshower_energies"][0]) if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_E_V[0]       = float(vtxShowerData["secondshower_energies"][1]) if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_E_Y[0]       = float(vtxShowerData["secondshower_energies"][2]) if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_gap_U[0]     = int(vtxShowerData["shower_gap"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_gap_V[0]     = int(vtxShowerData["shower_gap"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_gap_Y[0]     = int(vtxShowerData["shower_gap"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_dir_3d_X[0]  = float(vtxShowerData["shower_direction_3d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_dir_3d_Y[0]  = float(vtxShowerData["shower_direction_3d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_dir_3d_Z[0]  = float(vtxShowerData["shower_direction_3d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_dir_2d_U[0]  = float(vtxShowerData["shower_direction_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_dir_2d_V[0]  = float(vtxShowerData["shower_direction_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_dir_2d_Y[0]  = float(vtxShowerData["shower_direction_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_op_2d_U[0]   = float(vtxShowerData["shower_opening_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_op_2d_V[0]   = float(vtxShowerData["shower_opening_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_op_2d_Y[0]   = float(vtxShowerData["shower_opening_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_U_X[0] = int(vtxShowerData["shower_start_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_U_Y[0] = int(vtxShowerData["shower_start_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_U_Z[0] = int(vtxShowerData["shower_start_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_V_X[0] = int(vtxShowerData["shower_start_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_V_Y[0] = int(vtxShowerData["shower_start_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_V_Z[0] = int(vtxShowerData["shower_start_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_Y_X[0] = int(vtxShowerData["shower_start_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_Y_Y[0] = int(vtxShowerData["shower_start_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_start_2d_Y_Z[0] = int(vtxShowerData["shower_start_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower1_impact[0]   = float(vtxShowerData["shower_impact"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_gap_U[0]    = int(vtxShowerData["secondshower_gap"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_gap_V[0]    = int(vtxShowerData["secondshower_gap"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_gap_Y[0]    = int(vtxShowerData["secondshower_gap"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_dir_3d_X[0] = float(vtxShowerData["secondshower_direction_3d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_dir_3d_Y[0] = float(vtxShowerData["secondshower_direction_3d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_dir_3d_Z[0] = float(vtxShowerData["secondshower_direction_3d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_dir_2d_U[0] = float(vtxShowerData["secondshower_direction_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_dir_2d_V[0] = float(vtxShowerData["secondshower_direction_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_dir_2d_Y[0] = float(vtxShowerData["secondshower_direction_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_op_2d_U[0]  = float(vtxShowerData["secondshower_opening_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_op_2d_V[0]  = float(vtxShowerData["secondshower_opening_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_op_2d_Y[0]  = float(vtxShowerData["secondshower_opening_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_U_X[0] = int(vtxShowerData["secondshower_start_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_U_Z[0] = int(vtxShowerData["secondshower_start_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_U_Y[0] = int(vtxShowerData["secondshower_start_2d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_V_X[0] = int(vtxShowerData["secondshower_start_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_V_Y[0] = int(vtxShowerData["secondshower_start_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_V_Z[0] = int(vtxShowerData["secondshower_start_2d"][1])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_Y_X[0] = int(vtxShowerData["secondshower_start_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_Y_Y[0] = int(vtxShowerData["secondshower_start_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_start_2d_Y_Z[0] = int(vtxShowerData["secondshower_start_2d"][2])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower2_impact[0]  = float(vtxShowerData["secondshower_impact"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._shower_alpha[0]    = float(vtxShowerData["opening_angle_3d"][0])       if PassSimpleCuts and PassShowerReco else -99999
    dlanavars._pi0mass[0]         = float(vtxShowerData["pi0mass"][0])       if PassSimpleCuts and PassShowerReco else -99999

    dlanavars._totPE[0] = PMTPrecut_Dict[IDev]['_totpe']
    dlanavars._porchTotPE[0] = PMTPrecut_Dict[IDev]['_porchtotpe']
    dlanavars._maxPEFrac[0] = PMTPrecut_Dict[IDev]['_maxpefrac']
    dlanavars._passPMTPrecut[0] = PMTPrecut_Dict[IDev]['_passpmtprecut']
    dlanavars._parentPDG[0] = MC_dict[IDev]['parentPDG']                  if ismc else -99998
    dlanavars._energyInit[0] = MC_dict[IDev]['energyInit']                if ismc else -99998
    dlanavars._nproton[0] = MC_dict[IDev]['nproton']                      if ismc else -99998
    dlanavars._nlepton[0] = MC_dict[IDev]['nlepton']                      if ismc else -99998
    dlanavars._parentX[0] = parentX                                       if ismc else -99998
    dlanavars._parentY[0] = parentY                                       if ismc else -99998
    dlanavars._parentZ[0] = parentZ                                       if ismc else -99998
    dlanavars._parentSCEX[0] = parentSCEX                                 if ismc and PassSimpleCuts else -99998
    dlanavars._parentSCEY[0] = parentSCEY                                 if ismc and PassSimpleCuts else -99998
    dlanavars._parentSCEZ[0] = parentSCEZ                                 if ismc and PassSimpleCuts else -99998
    dlanavars._scedr[0] = scedr                                           if ismc and PassSimpleCuts else -99998

    dlanavars._eminusPID_int_v.clear()
    dlanavars._muonPID_int_v.clear()
    dlanavars._protonPID_int_v.clear()
    dlanavars._gammaPID_int_v.clear()
    dlanavars._pionPID_int_v.clear()
    dlanavars._eminusPID_pix_v.clear()
    dlanavars._muonPID_pix_v.clear()
    dlanavars._protonPID_pix_v.clear()
    dlanavars._gammaPID_pix_v.clear()
    dlanavars._pionPID_pix_v.clear()
    for i in ev.eminus_int_score_torch:  dlanavars._eminusPID_int_v.push_back(i)
    for i in ev.muon_int_score_torch:    dlanavars._muonPID_int_v.push_back(i)
    for i in ev.proton_int_score_torch:  dlanavars._protonPID_int_v.push_back(i)
    for i in ev.gamma_int_score_torch:   dlanavars._gammaPID_int_v.push_back(i)
    for i in ev.pion_int_score_torch:    dlanavars._pionPID_int_v.push_back(i)
    for i in ev.eminus_pix_score_torch:  dlanavars._eminusPID_pix_v.push_back(i)
    for i in ev.muon_pix_score_torch:    dlanavars._muonPID_pix_v.push_back(i)
    for i in ev.proton_pix_score_torch:  dlanavars._protonPID_pix_v.push_back(i)
    for i in ev.gamma_pix_score_torch:   dlanavars._gammaPID_pix_v.push_back(i)
    for i in ev.pion_pix_score_torch:    dlanavars._pionPID_pix_v.push_back(i)

    print "[MakeSelectionVars] Done"
    return
