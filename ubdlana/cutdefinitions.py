from __future__ import print_function
import os,sys
import numpy as np


def precuts(x,run,cutMode=0):
    """ 
    copied from Nick's SelectionHelper.py 
    inputs
    ------
    x: Instance of SelectionVars class defined in varutils.py
    run: Data run
    cutMode: not used for now
    """

    if x.PassSimpleCuts == 0: return False
    if x.PassShowerReco ==0: return False
    if (x.TotPE < 20 or x.PorchTotPE > 20): return False
    if x.Proton_Edep < 50 or x.Electron_Edep < 35: return False
    if max(x.MaxShrFrac,-1) < 0.2: return False
    if x.ShowerConsistency > 2: return False
    if x.OpenAng < 0.5: return False
    if x.FailedBoost_1e1p: return False
    if x.Proton_ThetaReco > np.pi/2: return False
    #if not x.run in goodruns: return False
        
    if cutMode==1:
        if x.Enu_1e1p < 200 or x.Enu_1e1p > 1200: return False
        if x.Lepton_PhiReco < np.pi/2.0 + 0.25 and x.Lepton_PhiReco > np.pi/2.0 - 0.25: return False
        if x.Thetas < 0.75 or x.Thetas > 2.5: return False
        if x.Q0_1e1p > 350: return False
        if x.Q3_1e1p < 250 or x.Q3_1e1p > 850: return False
        if x.ProtonPID_int_plane2<0.1: return False  
    if cutMode==2:
        if x.Enu_1e1p < 700: return False
        if x.Enu_1e1p > 1200: return False
    if cutMode==3:
        if x.BDTscore_1e1p > 0.7 or x.BDTscore_1e1p < 0.01: return False
        if x.ProtonPID_int_plane2 < 0.0: return False
        if abs(x.Lepton_PhiReco-np.pi/2.0) < 0.25: return False
        
    return True

def postcuts(x,cutMode):
    if x.GammaPID_pix_plane2/(x.EminusPID_pix_plane2+0.0001) > 2: return False
    if x.Electron_Edep > 100 and x.MuonPID_int_plane2 > 0.2: return False
    if x.Pi0Mass > 50: return False
    #if cutMode in [0,2] and x.ProtonPID_pix_v[2] < 0.2: return False
    return True
