from __future__ import print_function
from math import sqrt
from SelectionDefs import NewAng, VtxInSimpleFid, VtxInFid, GetPhiT, pTrans,pTransRat
from SelectionDefs import alphaT, ECCQE, ECal, Q2, OpenAngle, PhiDiff, edgeCut
from SelectionDefs import ECCQE_mom, Boost, BoostTracks, Getpz, Getq3q0

def getBjXYboost(x,Electron_Edep):
    try:
        (pEB_1e1p,eEB_1e1p,pThB_1e1p,
         eThB_1e1p,pPhB_1e1p,ePhB_1e1p,
         EcalB_1e1p,EpCCQEB_1e1p,
         EeCCQEB_1e1p,sphB_1e1p) = BoostTracks(Electron_Edep,x.Proton_Edep,
                                               x.Lepton_ThetaReco,x.Proton_ThetaReco,
                                               x.Lepton_PhiReco,x.Proton_PhiReco,'electron')
        Q2calB_1e1p          = Q2(EcalB_1e1p,eEB_1e1p,eThB_1e1p)
        EHadB_1e1p           = (EcalB_1e1p - eEB_1e1p - .511)
        yB_1e1p              = EHadB_1e1p/EcalB_1e1p
        xB_1e1p              = Q2calB_1e1p/(2*939.5654*EHadB_1e1p)
        return xB_1e1p,yB_1e1p
    except:
        return -9999,-9999

def getSphB_boost(x,Electron_Edep):
    try:
        # we catch the calculations because sometimes we take sqrt of a negative number
        #print("sphb_boost input: ",Electron_Edep,x.Proton_Edep,
        #      x.Lepton_ThetaReco,x.Proton_ThetaReco,
        #      x.Lepton_PhiReco,x.Proton_PhiReco)

        (pEB_1e1p,eEB_1e1p,pThB_1e1p,
         eThB_1e1p,pPhB_1e1p,ePhB_1e1p,
         EcalB_1e1p,EpCCQEB_1e1p,
         EeCCQEB_1e1p,sphB_1e1p) = BoostTracks(Electron_Edep,x.Proton_Edep,
                                               x.Lepton_ThetaReco,x.Proton_ThetaReco,
                                               x.Lepton_PhiReco,x.Proton_PhiReco,'electron')
        dEp_1e1p            = EpCCQEB_1e1p - EcalB_1e1p
        dEe_1e1p            = EeCCQEB_1e1p - EcalB_1e1p
        dEep_1e1p           = EpCCQEB_1e1p - EeCCQEB_1e1p
        #print("sphB: ",sphB_1e1p," vs ",sqrt(dEp_1e1p**2+dEe_1e1p**2+dEep_1e1p**2))
        
        return sqrt(dEp_1e1p**2+dEe_1e1p**2+dEep_1e1p**2)
    except:
        return -9999


def getNewShowerCalibTrainingVarbs(x,newCalib=True,
                                   newCalib_m=0.01255796, newCalib_b=0.0,
                                   BE = 29.5):
    """
    From Nick's BDTHelper.
    newCalib_m, newCalib_b: shower pixelsum to MeV conversion factors
    BE: binding energy in MeV
    """

    if newCalib:
        Electron_Edep       = x.shower1_sumQ_Y*newCalib_m + newCalib_b
        Enu_1e1p            = ECal(Electron_Edep,x.Proton_Edep,'electron',B=BE)
        PT_1e1p             = pTrans(Electron_Edep,x.Proton_Edep,
                                     x.Lepton_ThetaReco,x.Proton_ThetaReco,
                                     x.Lepton_PhiReco,x.Proton_PhiReco,'electron')
        AlphaT_1e1p         = alphaT(Electron_Edep,x.Proton_Edep,
                                     x.Lepton_ThetaReco,x.Proton_ThetaReco,
                                     x.Lepton_PhiReco,x.Proton_PhiReco,'electron')
        PzEnu_1e1p          = Getpz(Electron_Edep,x.Proton_Edep,
                                    x.Lepton_ThetaReco,x.Proton_ThetaReco,'electron') - Enu_1e1p
        Q3_1e1p,Q0_1e1p     = Getq3q0(x.Proton_Edep,Electron_Edep,
                                      x.Proton_ThetaReco,x.Lepton_ThetaReco,
                                      x.Proton_PhiReco,x.Lepton_PhiReco,'electron',B=BE)
        pTRat_1e1p          = pTransRat(Electron_Edep,x.Proton_Edep,
                                        x.Lepton_ThetaReco,x.Proton_ThetaReco,
                                        x.Lepton_PhiReco,x.Proton_PhiReco,'electron')
        BjXB_1e1p,BjYB_1e1p    = getBjXYboost(x,Electron_Edep)
    else:
        Electron_Edep       = x.Electron_Edep
        Enu_1e1p            = x.Enu_1e1p
        PT_1e1p             = x.PT_1e1p
        AlphaT_1e1p         = x.AlphaT_1e1p
        PzEnu_1e1p          = x.PzEnu_1e1p
        Q3_1e1p,Q0_1e1p     = x.Q3_1e1p,x.Q0_1e1p
        pTRat_1e1p          = x.PTRat_1e1p
        BjXB_1e1p,BjYB_1e1p   = x.BjXB_1e1p,x.BjYB_1e1p

    EpCCQE              = ECCQE(x.Proton_Edep,x.Proton_ThetaReco,pid="proton",B=BE)
    EeCCQE              = ECCQE(Electron_Edep,x.Lepton_ThetaReco,pid="electron",B=BE)
    SphB_1e1p           = getSphB_boost(x,Electron_Edep) 

    #Standard varbs
    training_varbs = [Enu_1e1p, Electron_Edep, PT_1e1p, AlphaT_1e1p, SphB_1e1p, PzEnu_1e1p, x.ChargeNearTrunk, 
                      Q0_1e1p, Q3_1e1p, x.Thetas, x.Phis, pTRat_1e1p, x.Proton_TrackLength, x.Lepton_TrackLength, 
                      x.Proton_ThetaReco, x.Proton_PhiReco, x.Lepton_ThetaReco, x.Lepton_PhiReco, max(x.MinShrFrac,-1),
                      max(x.MaxShrFrac,-1), x.shower1_smallQ_Y/(x.shower1_sumQ_Y+1e-6), BjXB_1e1p, BjYB_1e1p]
            
    return training_varbs

def get_bdt_var_names():
    names = \
        ["Enu_1e1p", 
         "Electron_Edep", 
         "PT_1e1p", 
         "AlphaT_1e1p", 
         "SphB_1e1p", 
         "PzEnu_1e1p", 
         "ChargeNearTrunk", 
         "Q0_1e1p", 
         "Q3_1e1p", 
         "Thetas", 
         "Phis", 
         "pTRat_1e1p", 
         "Proton_TrackLength",
         "Lepton_TrackLength", 
         "Proton_ThetaReco", 
         "Proton_PhiReco", 
         "Lepton_ThetaReco", 
         "Lepton_PhiReco", 
         "MinShrFrac",
         "MaxShrFrac", 
         None,
         "BjXB_1e1p", 
         "BjYB_1e1p"]
    return names

