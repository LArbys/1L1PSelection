# BDT Models 

Model weights go here.

List:

* `BDTweights_R[1/2/3]_[0-9].pickle`:  Ensemble trained by Josh Mills for 1m1p. Used in v1_1_4.

* `BDT_1e1p_Run[1/3]_6-2-20.pickle`: from Jarrett for 1e1p. Used in v1_1_3.

* `bdtweight_series2_june1_run[1/3].pickle`: from Davio for 1m1p. Used in v1_1_3.


## Vars for 1m1p BDTs

You can find these models being used in `ubdlana/bdt_util.py`.

```
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
```

## Vars for 1e1p BDTs

You can find these models being used in `ubdlana/bdt_util.py`.

```
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
```