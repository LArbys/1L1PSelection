import os,sys
from SelectionDefs import GetShCons

"""
Provide tools to extract selection variables from DLanaTree class (dlanatree.py) and 
the output FinalVertexVariables tree
"""


class SelectionVars:
    def __init__(self):
        # common precuts
        self.PassPMTPrecut = 0
        self.PassSimpleCuts = 0

        # 1e1p precuts
        self.PassShowerReco = 0
        self.TotPE = 0.0
        self.PorchTotPE = 0.0
        self.Proton_Edep = 0.0
        self.Electron_Edep = 0.0
        self.MaxShrFrac = 0.0
        self.ShowerConsistency = 0.0
        self.OpenAng = 0.0
        self.FailedBoost_1e1p = 0.0
        self.Proton_ThetaReco = 0.0

        # 1e1p signal variables
        self.Proton_Edep = 0.0
        self.Electron_Edep = 0.0
        self.Pi0Mass = 0.0
        self.GammaPID_pix_plane2 = 0.0
        self.EminusPID_pix_plane2 = 0.0
        self.ProtonPID_pix_plane2 = 0.0
        self.MuonPID_int_plane2 = 0.0
        self.Enu_1e1p = 0.0

        
def make_selection_vars_from_dlanatree( dlvars ):
    """ using variables in DLanaTree """
    self.PassSimpleCuts = dlvars._passCuts[0]
    self.TotPE = dlvars._totPE[0]
    self.PorchTotPE = dlvars._porchTotPE[0]

def make_selection_vars_from_fvv( fvv ):

    selvars = SelectionVars()

    # common precuts
    selvars.PassPMTPrecut  = fvv.PassPMTPrecut
    selvars.PassSimpleCuts = fvv.PassSimpleCuts

    # 1e1p precuts
    selvars.PassShowerReco = fvv.PassShowerReco
    selvars.TotPE = fvv.TotPE
    selvars.PorchTotPE = fvv.PorchTotPE
    selvars.Proton_Edep = fvv.Proton_Edep
    selvars.Electron_Edep = fvv.Electron_Edep
    selvars.MaxShrFrac = fvv.MaxShrFrac
    #selvars.ShowerConsistency = GetShCons( fvv.shower1_sumQ_V, fvv.shower1_sumQ_U, fvv.shower1_sumQ_Y )
    selvars.ShowerConsistency = GetShCons( fvv.shower1_sumQ_U, fvv.shower1_sumQ_V, fvv.shower1_sumQ_Y )
    selvars.OpenAng = fvv.OpenAng
    selvars.FailedBoost_1e1p = fvv.FailedBoost_1e1p
    selvars.Proton_ThetaReco = fvv.Proton_ThetaReco
    
    # 1e1p signal variables
    selvars.Proton_Edep = fvv.Proton_Edep
    selvars.Electron_Edep = fvv.Electron_Edep
    selvars.Pi0Mass = fvv._pi0mass
    selvars.GammaPID_pix_plane2 = fvv.GammaPID_pix_v[2]
    selvars.EminusPID_pix_plane2 = fvv.EminusPID_pix_v[2]
    selvars.ProtonPID_pix_plane2 = fvv.ProtonPID_pix_v[2]
    selvars.MuonPID_int_plane2 = fvv.MuonPID_int_v[2]
    selvars.BDTscore_1e1p = fvv.BDTscore_1e1p
    selvars.Enu_1e1p = fvv.Enu_1e1p

    return selvars

def update_bdt1e1p(x,rsev,bdt_results):
    """
    Update the values in an instance of SelectionVars for BDT1e1p values.
    Intended for when we rerun the 1e1p bdt on dlana output files.
    """
    x.BDTscore_1e1p = bdt_results[rsev]

def update_mpid(x,rsev,mpid_results):
    """
    Update the values in an instance of SelectionVars for MPID values.
    Intended for when we rerun the MPID on dlana output files
    """
    mpid_data = mpid_results[rsev] # get mpid data for RSEV entry
    x.GammaPID_pix_plane2  = mpid_data[(2,"pix")][1]
    x.EminusPID_pix_plane2 = mpid_data[(2,"pix")][0]
    x.ProtonPID_pix_plane2 = mpid_data[(2,"pix")][4]
    x.MuonPID_int_plane2   = mpid_data[(2,"int")][2]

def update_pmt_precuts(x,rse,pmtprecut_results):
    """
    Update the values in an instance of SelectionVars for PMT precut variables.
    Intended for when we rerun the PMT Precuts on dlana output files.
    """
    x.PassPMTPrecut = 1 if self.PMTPrecut_Dict[rse]['_passpmtprecut'] else 0
    
    

    
    
