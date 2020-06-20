import os,sys,pickle
import scipy
import xgboost
import numpy as np

"""
utility functions to deploy BDT models on the selection variables
"""


def load_BDT_model( weightfile_name ):
    """ model is loaded from a pickle """
    weight_path = os.environ["UBDLANA_DIR"]+"/bdt_models/"+weightfile_name
    with open(weight_path,"rb") as handle: 
        MyBDT = pickle.load(handle)
    print "[bdtutil::load_BDT_model] loaded ",weight_path
    return MyBDT

def apply_1e1p_model( model, dlvars ):
    """ apply the 1e1p BDT 
    
    inputs
    ------
    model: xgboost model loaded from pickle file
    dlvars: instance of DLanaTree (defined dlanatree.py module)

    outputs
    -------
    score [float] 1e1p score
    """

    shower_charge_ratio = dlvars._shower1_smallq_Y[0]/dlvars._shower1_sumq_Y[0] if dlvars._shower1_sumq_Y[0]>0 else 0.0

    # make input vars
    input_vars = [[
        dlvars._enu_1e1p[0],    # Enu_1e1p
        dlvars._electron_E[0],  # Electron_Edep,
        dlvars._eta[0],         # Eta,
        dlvars._pT_1e1p[0],     # PT_1e1p,
        dlvars._alphaT_1e1p[0], # AlphaT_1e1p,
        dlvars._sphB_1e1p[0],   # SphB_1e1p,
        dlvars._pzEnu_1e1p[0],  # PzEnu_1e1p,
        dlvars._charge_near_trunk[0]*dlvars._qcorrection_factor_vtx[0], #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
        dlvars._q0_1e1p[0],     # Q0_1e1p,
        dlvars._q3_1e1p[0],     # Q3_1e1p,
        dlvars._thetas[0],      # Thetas,
        dlvars._phis[0],        # Phis,
        dlvars._pTRat_1e1p[0],  # PTRat_1e1p,
        dlvars._proton_length[0], # Proton_TrackLength,
        dlvars._lepton_length[0], # Lepton_TrackLength,
        dlvars._proton_theta[0],  # Proton_ThetaReco,
        dlvars._proton_phi[0],    # Proton_PhiReco,
        dlvars._lepton_theta[0],  # Lepton_ThetaReco,
        dlvars._lepton_phi[0],    # Lepton_PhiReco,
        max(dlvars._minshrFrac[0],-1), # MinShrFrac
        max(dlvars._maxshrFrac[0],-1), # MaxShrFrac
        max(shower_charge_ratio,-1.0), # SHOWER_CHARGE_IN_SHOWER_CLUSTER/SHOWER_CHARGE_IN_IMAGE
        dlvars._openAng[0] ]] #OpenAng ]

    vars_np = np.asarray( input_vars )
    #print vars_np

    probs   = model.predict_proba(vars_np)[0]
    print "BDT[1e1p] output: ",probs
    return probs

def apply_1mu1p_model( model, dlvars ):
    """ apply the 1mu1p cosmic and nu BDTs
    
    inputs
    ------
    model: xgboost model loaded from pickle file
    dlvars: instance of DLanaTree (defined dlanatree.py module)

    outputs
    -------
    score [float] 1m1p bkg score
    """

    """
vars_june1 = ['Phis',
'ChargeNearTrunk',
'Enu_1m1p',
'PhiT_1m1p',
'AlphaT_1m1p',
'PT_1m1p',
'PTRat_1m1p',
'BjXB_1m1p',
'BjYB_1m1p',
'SphB_1m1p',
'Q0_1m1p',
'Q3_1m1p',
'Lepton_PhiReco',
'Lepton_TrackLength',
'Proton_PhiReco',
'Proton_ThetaReco']
    """

    # make input vars
    input_vars = [[
        dlvars._phis[0],          # Phis,
        dlvars._charge_near_trunk[0], #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
        dlvars._enu_1m1p[0],				# reco nu energy
        dlvars._phiT_1m1p[0],     # PhiT_1m1p
        dlvars._alphaT_1m1p[0],   # AlphaT_1m1p,
        dlvars._pT_1m1p[0],       # PT_1m1p,
        dlvars._pTRat_1m1p[0],    # PTRat_1m1p,
        dlvars._bjXB_1m1p[0],
        dlvars._bjYB_1m1p[0],
        dlvars._sphB_1m1p[0],     # Sph_1m1p,
        dlvars._q0_1m1p[0],       # Q0_1m1p,
        dlvars._q3_1m1p[0],       # Q3_1m1p,        
        dlvars._lepton_phi[0],    # Lepton_PhiReco,
        dlvars._lepton_length[0], # Lepton_TrackLength,
        dlvars._proton_phi[0],    # Proton_PhiReco,
        dlvars._proton_theta[0], # Proton_ThetaReco,
        ]]

    vars_np = np.asarray( input_vars )
    #print vars_np

    probs = model.predict_proba(vars_np)[0]    
    print "BDT[1mu1p] output: ",probs
    return probs

def rerun_1mu1p_models( model, fvv ):
    """ apply the 1mu1p cosmic and nu BDTs
    
    inputs
    ------
    cosmicBDT: xgboost model loaded from pickle file
    nuBDT:     xgboost model loaded from pickle file
    dlvars:    finalvertexvariable tree (the DL ana tree)

    outputs
    -------
    score [float] 1m1p bdt score
    score [float] 1m1p nu score
    """

    """
vars_june1 = ['Phis',
'ChargeNearTrunk',
'Enu_1m1p',
'PhiT_1m1p',
'AlphaT_1m1p',
'PT_1m1p',
'PTRat_1m1p',
'BjXB_1m1p',
'BjYB_1m1p',
'SphB_1m1p',
'Q0_1m1p',
'Q3_1m1p',
'Lepton_PhiReco',
'Lepton_TrackLength',
'Proton_PhiReco',
'Proton_ThetaReco']
    """

    bdtout_dict = {}
    nentries = fvv.GetEntries()
    print "[bdtutil::rerun_1mu1p_models] rerun on ",nentries
    for ientry in range(nentries):
        fvv.GetEntry(ientry)

        rsev = (fvv.run,fvv.subrun,fvv.event,fvv.vtxid)
        
        # make input vars
        input_vars = [[
            fvv.Phis,
            fvv.ChargeNearTrunk,
            fvv.Enu_1m1p,
            fvv.PhiT_1m1p,
            fvv.AlphaT_1m1p,
            fvv.PT_1m1p,      
            fvv.PTRat_1m1p,
            fvv.BjXB_1m1p,        
            fvv.BjYB_1m1p,
            fvv.SphB_1m1p,
            fvv.Q0_1m1p,
            fvv.Q3_1m1p,        
            fvv.Lepton_PhiReco,
            fvv.Lepton_TrackLength,
            fvv.Proton_PhiReco,
            fvv.Proton_ThetaReco,
        ]]

        vars_np = np.asarray( input_vars )
        #print vars_np

        probs   = model.predict_proba(vars_np)[0]
        print "[bdtutil::rerun_1m1p_bdt] rsev=(",rsev,")"
        print "  BDT[1mu1p] output: ",probs
        bdtout_dict[rsev] = probs[0]
        
    return bdtout_dict

def rerun_1e1p_models( model, fvv ):
    """ apply the 1e1p BDT
    
    inputs
    ------
    1e1p:      xgboost model loaded from pickle file
    dlvars:    finalvertexvariable tree (the DL ana tree)

    outputs
    -------
    1e1p bdt score dictionary
    """

    """
[Enu_1e1p,
Electron_Edep,
PT_1e1p,
AlphaT_1e1p,
SphB_1e1p,
PzEnu_1e1p,
ChargeNearTrunk,
Q0_1e1p,
Q3_1e1p,
Thetas,
Phis,
PTRat_1e1p,
Proton_TrackLength,
Lepton_TrackLength,
Proton_ThetaReco,
Proton_PhiReco,
Lepton_ThetaReco,
Lepton_PhiReco,
MinShrFrac,
MaxShrFrac,
shower1_smallQ_Y/shower1_sumQ_Y ,
BjX_1e1p,
BjY_1e1p]
    """

    bdtout_dict = {}
    nentries = fvv.GetEntries()
    print "[bdtutil::rerun_1e1p_models] rerun on ",nentries
    for ientry in range(nentries):
        fvv.GetEntry(ientry)

        rsev = (fvv.run,fvv.subrun,fvv.event,fvv.vtxid)
        
        shower_charge_ratio = fvv.shower1_smallQ_Y/fvv.shower1_sumQ_Y if fvv.shower1_sumQ_Y>0 else 0.0

        # make input vars
        input_vars = [[
                fvv.Enu_1e1p,
                fvv.Electron_Edep,
                fvv.PT_1e1p,
                fvv.AlphaT_1e1p,
                fvv.SphB_1e1p,
                fvv.PzEnu_1e1p,
                #fvv.ChargeNearTrunk*fvv.QCorrectionFactorVertex,
                fvv.ChargeNearTrunk,
                fvv.Q0_1e1p,
                fvv.Q3_1e1p,
                fvv.Thetas,
                fvv.Phis,
                fvv.PTRat_1e1p,
                fvv.Proton_TrackLength,
                fvv.Lepton_TrackLength,
                fvv.Proton_ThetaReco,
                fvv.Proton_PhiReco,
                fvv.Lepton_ThetaReco,
                fvv.Lepton_PhiReco,
                fvv.MinShrFrac,
                fvv.MaxShrFrac,
                shower_charge_ratio,
                fvv.BjX_1e1p,
                fvv.BjY_1e1p]]

        vars_np = np.asarray( input_vars )
        #print vars_np

        probs   = model.predict_proba(vars_np)[0]
        print "[bdtutil::rerun_1e1p_bdt] rsev=(",rsev,")"
        print "  BDT[1e1p] output: ",probs
        bdtout_dict[rsev] = probs[0]
        
    return bdtout_dict
