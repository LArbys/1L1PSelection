import os,sys,pickle
import scipy
import xgboost
import numpy as np
import bdt1e1p_helper

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

def load_BDT_ensemble( model_type, bdt_dir, nbdts=10, runs=[1,2,3] ):
    """ used to load bdt ensembles used in v1_1_4 ubdlana code"""
    MODELS=  ["1e1p","1m1p"]

    if model_type not in ["1e1p", "1m1p"]:
        raise ValueError("[ubdlana::bdtutil::load_BDT_ensemble] Did not recognize model, \"{}\". Choices: {}".format(model_type,MODELS))

    bdt_dict = {}
    for run in runs:
        bdt_dict[run] = {}
        for b in range(nbdts):
            if model_type == "1m1p":
                print "Load into ensemble: ",'BDTweights_R%i_%i_py2.pickle'%(run,b)
                bdt_dict[run][b] = pickle.load(open(bdt_dir+'/BDTweights_R%i_%i_py2.pickle'%(run,b),'rb'))
            elif model_type=="1e1p":
                print "Load into ensemble: ",'BDTweights_1e1p_R%i_%i_py2.pickle'%(run,b)
                bdt_dict[run][b] = pickle.load(open(bdt_dir+'/BDTweights_1e1p_R%i_%i_py2.pickle'%(run,b),'rb'))
            else:
                raise ValueError("no model")
    return bdt_dict

    

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
    # input_vars = [[
    #     dlvars._enu_1e1p[0],    # Enu_1e1p
    #     dlvars._electron_E[0],  # Electron_Edep,
    #     dlvars._eta[0],         # Eta,
    #     dlvars._pT_1e1p[0],     # PT_1e1p,
    #     dlvars._alphaT_1e1p[0], # AlphaT_1e1p,
    #     dlvars._sphB_1e1p[0],   # SphB_1e1p,
    #     dlvars._pzEnu_1e1p[0],  # PzEnu_1e1p,
    #     dlvars._charge_near_trunk[0]*dlvars._qcorrection_factor_vtx[0], #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
    #     dlvars._q0_1e1p[0],     # Q0_1e1p,
    #     dlvars._q3_1e1p[0],     # Q3_1e1p,
    #     dlvars._thetas[0],      # Thetas,
    #     dlvars._phis[0],        # Phis,
    #     dlvars._pTRat_1e1p[0],  # PTRat_1e1p,
    #     dlvars._proton_length[0], # Proton_TrackLength,
    #     dlvars._lepton_length[0], # Lepton_TrackLength,
    #     dlvars._proton_theta[0],  # Proton_ThetaReco,
    #     dlvars._proton_phi[0],    # Proton_PhiReco,
    #     dlvars._lepton_theta[0],  # Lepton_ThetaReco,
    #     dlvars._lepton_phi[0],    # Lepton_PhiReco,
    #     max(dlvars._minshrFrac[0],-1), # MinShrFrac
    #     max(dlvars._maxshrFrac[0],-1), # MaxShrFrac
    #     max(shower_charge_ratio,-1.0), # SHOWER_CHARGE_IN_SHOWER_CLUSTER/SHOWER_CHARGE_IN_IMAGE
    #     dlvars._openAng[0] ]] #OpenAng ]
    input_vars = [[
        dlvars._enu_1e1p[0],    # Enu_1e1p
        dlvars._electron_E[0],  # Electron_Edep,
        dlvars._pT_1e1p[0],     # PT_1e1p,
        dlvars._alphaT_1e1p[0], # AlphaT_1e1p,
        dlvars._sphB_1e1p[0],   # SphB_1e1p,
        dlvars._pzEnu_1e1p[0],  # PzEnu_1e1p,
        #dlvars._charge_near_trunk[0]*dlvars._qcorrection_factor_vtx[0], #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
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
        max(shower_charge_ratio,-1.0), # SHOWER_CHARGE_IN_SHOWER_CLUSTER/SHOWER_CHARGE_IN_IMAGE
        dlvars._bjX_1e1p[0],   #BjX_1e1p ]
        dlvars._bjY_1e1p[0] ]] #BjY_1e1p ]

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

def apply_1mu1p_ensemble_model( model, dlvars, DATARUN, nbdts=10 ):
    """ apply the 1mu1p BDT ensemble
    
    inputs
    ------
    model: dictionary of xgboost model loaded from pickle file. See load_BDT_ensemble above.
    dlvars: instance of DLanaTree (defined dlanatree.py module)

    outputs
    -------
    score [float] 1m1p bkg score
    """

    # make input vars
    input_vars = [[
        dlvars._phis[0],          # Phis,
        dlvars._charge_near_trunk[0], #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
        dlvars._enu_1m1p[0],	  # reco nu energy
        dlvars._phiT_1m1p[0],     # PhiT_1m1p
        dlvars._alphaT_1m1p[0],   # AlphaT_1m1p,
        dlvars._pT_1m1p[0],       # PT_1m1p,
        dlvars._pTRat_1m1p[0],    # PTRat_1m1p,
        dlvars._bjXB_1m1p[0],     # Bjorken X
        dlvars._bjYB_1m1p[0],     # Byorken Y
        dlvars._sphB_1m1p[0],     # Sph_1m1p,
        dlvars._q0_1m1p[0],       # Q0_1m1p,
        dlvars._q3_1m1p[0],       # Q3_1m1p,        
        dlvars._lepton_phi[0],    # Lepton_PhiReco,
        dlvars._lepton_length[0], # Lepton_TrackLength,
        dlvars._proton_phi[0],    # Proton_PhiReco,
        dlvars._proton_theta[0],  # Proton_ThetaReco,
        ]]

    
    vars_np = np.asarray( input_vars )

    scores = np.zeros((nbdts))
    for b in range(nbdts):
        sigprob = model[DATARUN][b].predict_proba(vars_np)[:,1]
        scores[b] = sigprob

    sigavg = np.average(scores)
    sigmedian = np.median(scores)
    sigmax = np.max(scores)
    
    #print vars_np
    print "BDT[1mu1p]-ENSEMBLE output: ave=%.1f median=%.1f max=%.1f"%(sigavg,sigmedian,sigmax)
    return {"ave":sigavg,"median":sigmedian,"max":sigmax}

class BDT1e1pVariables:
    def __init__(self, dlvars):

        self.Enu_1e1p = dlvars._enu_1e1p[0]         # Enu_1e1p
        self.Electron_Edep = dlvars._electron_E[0]  # Electron_Edep,
        self.Proton_Edep   = dlvars._proton_E[0]    # Proton_Edep,
        self.PT_1e1p = dlvars._pT_1e1p[0]           # PT_1e1p,
        self.AlphaT_1e1p = dlvars._alphaT_1e1p[0]   # AlphaT_1e1p,
        self.SphB_1e1p   = dlvars._sphB_1e1p[0]     # SphB_1e1p,
        self.PzEnu_1e1p  = dlvars._pzEnu_1e1p[0]    # PzEnu_1e1p,
        self.ChargeNearTrunk = dlvars._charge_near_trunk[0] #ChargeNearTrunk_UniformityCalibrated,  #!!! calibrated?
        self.Q0_1e1p = dlvars._q0_1e1p[0]           # Q0_1e1p,
        self.Q3_1e1p = dlvars._q3_1e1p[0]           # Q3_1e1p,
        self.Thetas  = dlvars._thetas[0]            # Thetas,
        self.Phis    = dlvars._phis[0]              # Phis,
        self.PTRat_1e1p = dlvars._pTRat_1e1p[0]     # PTRat_1e1p,
        self.BjX_1e1p = dlvars._bjX_1e1p[0]         #BjX_1e1p
        self.BhY_1e1p = dlvars._bjY_1e1p[0]         #BjY_1e1p ]
        self.Proton_TrackLength = dlvars._proton_length[0] # Proton_TrackLength,
        self.Lepton_TrackLength = dlvars._lepton_length[0] # Lepton_TrackLength,
        self.Proton_ThetaReco   = dlvars._proton_theta[0]  # Proton_ThetaReco,
        self.Proton_PhiReco     = dlvars._proton_phi[0]    # Proton_PhiReco,
        self.Lepton_ThetaReco   = dlvars._lepton_theta[0]  # Lepton_ThetaReco,
        self.Lepton_PhiReco     = dlvars._lepton_phi[0]    # Lepton_PhiReco,
        self.MinShrFrac         = max(dlvars._minshrFrac[0],-1) # MinShrFrac
        self.MaxShrFrac         = max(dlvars._maxshrFrac[0],-1) # MaxShrFrac
        self.shower1_sumQ_Y     = dlvars._shower1_sumq_Y[0] # shower1_sumQ_Y
        self.shower1_smallQ_Y   = dlvars._shower1_smallq_Y[0] # shower1_smallQ_Y
        

def apply_1e1p_ensemble_model( model, dlvars, DATARUN, nbdts=20, maxentries=None ):
    """ apply the 1e1p BDT
    
    inputs
    ------
    1e1p:      xgboost model loaded from pickle file
    fvv:    finalvertexvariable tree (the DL ana tree)

    outputs
    -------
    1e1p bdt score dictionary
    """
    rsev = (dlvars._run[0],dlvars._subrun[0],dlvars._event[0],dlvars._vtxid[0])
    fvv = BDT1e1pVariables(dlvars)
    
    # make input vars
    input_vars = bdt1e1p_helper.getNewShowerCalibTrainingVarbs( fvv, newCalib=True )
    vars_np = np.asarray( [input_vars] )

    scores = np.zeros((nbdts))
    for b in range(nbdts):
        sigprob = model[DATARUN][b].predict_proba(vars_np)[:,1]
        scores[b] = sigprob

    sigavg = np.average(scores)
    sigmedian = np.median(scores)
    sigmax = np.max(scores)
    print "[bdtutil::apply_1e1p_ensemble_model] rsev=(",rsev,") Electron_Edep=",input_vars[1]," SphB_1e1p=",input_vars[4]
    print "  BDT[1e1p]-ensemble ave=",sigavg," median=",sigmedian," max=",sigmax
        
    return {"ave":sigavg,"median":sigmedian,"max":sigmax}


def rerun_1mu1p_ensemble( model, fvv, DATARUN, nbdts=10, maxentries=None ):
    """ rerun the 1mu1p Ensemble BDT on a tree
    
    inputs
    ------
    model:  BDT model. See load_BDT_ensemble above.
    dlvars: finalvertexvariable tree (the DL ana tree)

    outputs
    -------
    score [float] 1m1p nu score
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

        scores = np.zeros((nbdts))
        for b in range(nbdts):
            sigprob = model[DATARUN][b].predict_proba(vars_np)[:,1]
            scores[b] = sigprob

        sigavg = np.average(scores)
        sigmedian = np.median(scores)
        sigmax = np.max(scores)
        print "[bdtutil::rerun_1m1p_bdt] rsev=(",rsev,")"
        print "  BDT[1mu1p]-ensemble ave=",sigavg," median=",sigmedian," max=",sigmax
        bdtout_dict[rsev] = sigavg

        if maxentries is not None and ientry+1==maxentries:
            break

        
    return bdtout_dict


def rerun_1e1p_ensemble( model, fvv, DATARUN, nbdts=20, maxentries=None ):
    """ apply the 1e1p BDT
    
    inputs
    ------
    1e1p:      xgboost model loaded from pickle file
    dlvars:    finalvertexvariable tree (the DL ana tree)

    outputs
    -------
    1e1p bdt score dictionary
    """

    bdtout_dict = {}
    electron_edep_dict = {}
    nentries = fvv.GetEntries()
    print "[bdtutil::rerun_1e1p_models] rerun on ",nentries
    for ientry in range(nentries):
        fvv.GetEntry(ientry)

        rsev = (fvv.run,fvv.subrun,fvv.event,fvv.vtxid)
        
        # make input vars
        input_vars = bdt1e1p_helper.getNewShowerCalibTrainingVarbs( fvv, newCalib=True )
        vars_np = np.asarray( [input_vars] )

        scores = np.zeros((nbdts))
        for b in range(nbdts):
            sigprob = model[DATARUN][b].predict_proba(vars_np)[:,1]
            scores[b] = sigprob

        sigavg = np.average(scores)
        sigmedian = np.median(scores)
        sigmax = np.max(scores)
        print "[bdtutil::rerun_1e1p_ensemble] rsev=(",rsev,") Electron_Edep=",input_vars[1]," SphB_1e1p=",input_vars[4]
        print "  BDT[1e1p]-ensemble ave=",sigavg," median=",sigmedian," max=",sigmax
        bdtout_dict[rsev] = sigavg
        electron_edep_dict[rsev] = input_vars[1]
        if maxentries is not None and ientry+1==maxentries:
            break
        
    return bdtout_dict, electron_edep_dict
