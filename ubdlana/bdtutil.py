import os,sys,pickle
import scipy
import xgboost
import numpy as np

"""
utility functions to deploy BDT models on the selection variables
"""


def load_1e1p_model( weightfile_name ):
    """ model is loaded from a pickle """
    weight_path = os.environ["UBDLANA_DIR"]+"/bdt_models/"+weightfile_name
    with open(weight_path,"rb") as handle: 
        MyBDT = pickle.load(handle)
    return MyBDT

def load_1mu1p_model():
    pass

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

    # make input vars
    input_vars = [[
        dlvars._enu_1e1p[0],    # Enu_1e1p
        dlvars._electron_E[0],  # Electron_Edep,
        dlvars._eta[0],         # Eta,
        dlvars._pT_1e1p[0],     # PT_1e1p,
        dlvars._alphaT_1e1p[0], # AlphaT_1e1p,
        dlvars._sphB_1e1p[0],   # SphB_1e1p,
        dlvars._pzEnu_1e1p[0],  # PzEnu_1e1p,
        dlvars._charge_near_trunk[0], #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
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
        -1.0, #SHOWER_CHARGE_IN_IMAGE / SHOWER_CHARGE_IN_SHOWER_CLUSTER,
        dlvars._openAng[0] ]] #OpenAng ]

    vars_np = np.asarray( input_vars )
    #print vars_np

    probs   = model.predict_proba(vars_np)[0]
    print "BDT[1e1p] output: ",probs
    return probs

def apply_1m1p_model( model, dlvars ):
    pass
