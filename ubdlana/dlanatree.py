from array import array
import ROOT

from maketreebranch import MakeTreeBranch

# --------------------------------------------- #
# DLanaTree
# --------------------------------------------- #
class DLanaTree:

    def __init__(self):
        self.outTree = None
        self.makeDLanaTree()


    # --------------------------------------------- #
    # Make DL ANA Tree
    # --------------------------------------------- #
    def makeDLanaTree(self):
        """ define output ana tree """
        if self.outTree is not None:
            raise RuntimeError("Already created tree and branches")

        outTree = ROOT.TTree('FinalVertexVariables','Final Vertex Variable Tree')
        ## lepton agnostic
        self._run = MakeTreeBranch(outTree,'run','int')
        self._subrun = MakeTreeBranch(outTree,'subrun','int')
        self._event = MakeTreeBranch(outTree,'event','int')
        self._vtxid = MakeTreeBranch(outTree,'vtxid','int')
        self._x = MakeTreeBranch(outTree,'Xreco','float')
        self._y = MakeTreeBranch(outTree,'Yreco','float')
        self._z = MakeTreeBranch(outTree,'Zreco','float')
        self._infiducial = MakeTreeBranch(outTree,'InFiducial','int')
        self._anyReco = MakeTreeBranch(outTree,'AnyReco','int')
        self._ntracks = MakeTreeBranch(outTree,'NTracks','int')
        self._n5tracks = MakeTreeBranch(outTree,'N5cmTracks','int')
        self._passCuts = MakeTreeBranch(outTree,'PassSimpleCuts','int')
        self._passShowerReco = MakeTreeBranch(outTree,'PassShowerReco','int')
        self._passSecShr = MakeTreeBranch(outTree,'PassSecondShower','int')
        self._failBoost = MakeTreeBranch(outTree,'FailedBoost','int')
        self._failBoost_1m1p = MakeTreeBranch(outTree,'FailedBoost_1m1p','int')
        self._failBoost_1e1p = MakeTreeBranch(outTree,'FailedBoost_1e1p','int')
        self._good3DReco = MakeTreeBranch(outTree,'Good3DReco','int')
        self._eta = MakeTreeBranch(outTree,'Eta','float')
        self._openAng = MakeTreeBranch(outTree,'OpenAng','float')
        self._thetas = MakeTreeBranch(outTree,'Thetas','float')
        self._phis = MakeTreeBranch(outTree,'Phis','float')
        self._qcorrection_factor_vtx = MakeTreeBranch(outTree,'QCorrectionFactorVertex','float')
        self._charge_near_trunk  = MakeTreeBranch(outTree,'ChargeNearTrunk','float')
        self._longtracklen = MakeTreeBranch(outTree,'LongTrackLen','float')
        self._shorttracklen = MakeTreeBranch(outTree,'ShortTrackLen','float')
        self._maxshrFrac = MakeTreeBranch(outTree,'MaxShrFrac','float')
        self._minshrFrac = MakeTreeBranch(outTree,'MinShrFrac','float')

        # 1mu1p stuff
        self._CCQE_energy_shift_1m1p = MakeTreeBranch(outTree,'CCQEEnergyShift_1m1p','float')
        self._enu_1m1p = MakeTreeBranch(outTree,'Enu_1m1p','float')
        self._phiT_1m1p = MakeTreeBranch(outTree,'PhiT_1m1p','float')
        self._alphaT_1m1p = MakeTreeBranch(outTree,'AlphaT_1m1p','float')
        self._pT_1m1p = MakeTreeBranch(outTree,'PT_1m1p','float')
        self._pTRat_1m1p = MakeTreeBranch(outTree,'PTRat_1m1p','float')
        self._bjX_1m1p = MakeTreeBranch(outTree,'BjX_1m1p','float')
        self._bjY_1m1p = MakeTreeBranch(outTree,'BjY_1m1p','float')
        self._q2_1m1p = MakeTreeBranch(outTree,'Q2_1m1p','float')
        self._sph_1m1p = MakeTreeBranch(outTree,'Sph_1m1p','float')
        self._pzEnu_1m1p = MakeTreeBranch(outTree,'PzEnu_1m1p','float')
        self._q0_1m1p = MakeTreeBranch(outTree,'Q0_1m1p','float')
        self._q3_1m1p = MakeTreeBranch(outTree,'Q3_1m1p','float')
        self._openAngB_1m1p = MakeTreeBranch(outTree,'OpenAngB_1m1p','float')
        self._thetasB_1m1p = MakeTreeBranch(outTree,'ThetasB_1m1p','float')
        self._phisB_1m1p = MakeTreeBranch(outTree,'PhisB_1m1p','float')
        self._phiTB_1m1p = MakeTreeBranch(outTree,'PhiTB_1m1p','float')
        self._alphaTB_1m1p = MakeTreeBranch(outTree,'AlphaTB_1m1p','float')
        self._pTB_1m1p = MakeTreeBranch(outTree,'PTB_1m1p','float')
        self._bjXB_1m1p = MakeTreeBranch(outTree,'BjXB_1m1p','float')
        self._bjYB_1m1p = MakeTreeBranch(outTree,'BjYB_1m1p','float')
        self._q2B_1m1p = MakeTreeBranch(outTree,'Q2B_1m1p','float')
        self._sphB_1m1p = MakeTreeBranch(outTree,'SphB_1m1p','float')

        # same stuff, but for 1e1p
        self._CCQE_energy_shift_1e1p = MakeTreeBranch(outTree,'CCQEEnergyShift_1e1p','float')
        self._enu_1e1p = MakeTreeBranch(outTree,'Enu_1e1p','float')
        self._phiT_1e1p = MakeTreeBranch(outTree,'PhiT_1e1p','float')
        self._alphaT_1e1p = MakeTreeBranch(outTree,'AlphaT_1e1p','float')
        self._pT_1e1p = MakeTreeBranch(outTree,'PT_1e1p','float')
        self._pTRat_1e1p = MakeTreeBranch(outTree,'PTRat_1e1p','float')
        self._bjX_1e1p = MakeTreeBranch(outTree,'BjX_1e1p','float')
        self._bjY_1e1p = MakeTreeBranch(outTree,'BjY_1e1p','float')
        self._q2_1e1p = MakeTreeBranch(outTree,'Q2_1e1p','float')
        self._sph_1e1p = MakeTreeBranch(outTree,'Sph_1e1p','float')
        self._pzEnu_1e1p = MakeTreeBranch(outTree,'PzEnu_1e1p','float')
        self._q0_1e1p = MakeTreeBranch(outTree,'Q0_1e1p','float')
        self._q3_1e1p = MakeTreeBranch(outTree,'Q3_1e1p','float')
        self._openAngB_1e1p = MakeTreeBranch(outTree,'OpenAngB_1e1p','float')
        self._thetasB_1e1p = MakeTreeBranch(outTree,'ThetasB_1e1p','float')
        self._phisB_1e1p = MakeTreeBranch(outTree,'PhisB_1e1p','float')
        self._phiTB_1e1p = MakeTreeBranch(outTree,'PhiTB_1e1p','float')
        self._alphaTB_1e1p = MakeTreeBranch(outTree,'AlphaTB_1e1p','float')
        self._pTB_1e1p = MakeTreeBranch(outTree,'PTB_1e1p','float')
        self._bjXB_1e1p = MakeTreeBranch(outTree,'BjXB_1e1p','float')
        self._bjYB_1e1p = MakeTreeBranch(outTree,'BjYB_1e1p','float')
        self._q2B_1e1p = MakeTreeBranch(outTree,'Q2B_1e1p','float')
        self._sphB_1e1p = MakeTreeBranch(outTree,'SphB_1e1p','float')

        self._lepton_id = MakeTreeBranch(outTree,'Lepton_ID','int')
        self._lepton_phi = MakeTreeBranch(outTree,'Lepton_PhiReco','float')
        self._lepton_theta = MakeTreeBranch(outTree,'Lepton_ThetaReco','float')
        self._lepton_length = MakeTreeBranch(outTree,'Lepton_TrackLength','float')
        self._lepton_dqdx_uncalibrated = MakeTreeBranch(outTree,'Lepton_dQdx_uncalibrated','float')
        self._lepton_dqdx_calibrated = MakeTreeBranch(outTree,'Lepton_dQdx','float')
        self._muon_E = MakeTreeBranch(outTree,'Muon_Edep','float')
        self._electron_E = MakeTreeBranch(outTree,'Electron_Edep','float')
        self._lepton_edge_dist = MakeTreeBranch(outTree,'Lepton_EdgeDist','float')
        self._muon_phiB_1m1p = MakeTreeBranch(outTree,'Muon_PhiRecoB_1m1p','float')
        self._muon_thetaB_1m1p = MakeTreeBranch(outTree,'Muon_ThetaRecoB_1m1p','float')
        self._muon_EB_1m1p = MakeTreeBranch(outTree,'Muon_EdepB_1m1p','float')
        self._electron_phiB_1e1p = MakeTreeBranch(outTree,'Electron_PhiRecoB_1e1p','float')
        self._electron_thetaB_1e1p = MakeTreeBranch(outTree,'Electron_ThetaRecoB_e1ep','float')
        self._electron_EB_1e1p = MakeTreeBranch(outTree,'Electron_EdepB_1e1p','float')
        self._proton_id = MakeTreeBranch(outTree,'Proton_ID','float')
        self._proton_phi = MakeTreeBranch(outTree,'Proton_PhiReco','float')
        self._proton_theta = MakeTreeBranch(outTree,'Proton_ThetaReco','float')
        self._proton_length = MakeTreeBranch(outTree,'Proton_TrackLength','float')
        self._proton_dqdx_uncalibrated = MakeTreeBranch(outTree,'Proton_dQdx_uncalibrated','float')
        self._proton_dqdx_calibrated = MakeTreeBranch(outTree,'Proton_dQdx','float')
        self._proton_E = MakeTreeBranch(outTree,'Proton_Edep','float')
        self._proton_edge_dist = MakeTreeBranch(outTree,'Proton_EdgeDist','float')
        self._proton_phiB_1m1p = MakeTreeBranch(outTree,'Proton_PhiRecoB_1m1p','float')
        self._proton_thetaB_1m1p = MakeTreeBranch(outTree,'Proton_ThetaRecoB_1m1p','float')
        self._proton_EB_1m1p = MakeTreeBranch(outTree,'Proton_EdepB_1m1p','float')
        self._proton_phiB_1e1p = MakeTreeBranch(outTree,'Proton_PhiRecoB_1e1p','float')
        self._proton_thetaB_1e1p = MakeTreeBranch(outTree,'Proton_ThetaRecoB_1e1p','float')
        self._proton_EB_1e1p = MakeTreeBranch(outTree,'Proton_EdepB_1e1p','float')

        # Shower variables
        self._shower1_E_U = MakeTreeBranch(outTree,'shower1_E_U','float')
        self._shower1_E_V = MakeTreeBranch(outTree,'shower1_E_V','float')
        self._shower1_E_Y = MakeTreeBranch(outTree,'shower1_E_Y','float')
        self._shower1_gap_U = MakeTreeBranch(outTree,'shower1_gap_U','int')
        self._shower1_gap_V = MakeTreeBranch(outTree,'shower1_gap_V','int')
        self._shower1_gap_Y = MakeTreeBranch(outTree,'shower1_gap_Y','int')
        self._shower1_dir_3d_X = MakeTreeBranch(outTree,'shower1_dir_3d_X','float')
        self._shower1_dir_3d_Y = MakeTreeBranch(outTree,'shower1_dir_3d_Y','float')
        self._shower1_dir_3d_Z = MakeTreeBranch(outTree,'shower1_dir_3d_Z','float')
        self._shower1_dir_2d_U = MakeTreeBranch(outTree,'shower1_dir_2d_U','float')
        self._shower1_dir_2d_V = MakeTreeBranch(outTree,'shower1_dir_2d_V','float')
        self._shower1_dir_2d_Y = MakeTreeBranch(outTree,'shower1_dir_2d_Y','float')
        self._shower1_op_2d_U = MakeTreeBranch(outTree,'shower1_op_2d_U','float')
        self._shower1_op_2d_V = MakeTreeBranch(outTree,'shower1_op_2d_V','float')
        self._shower1_op_2d_Y = MakeTreeBranch(outTree,'shower1_op_2d_Y','float')
        self._shower1_start_2d_U_X = MakeTreeBranch(outTree,'shower1_start_2d_U_X','int')
        self._shower1_start_2d_U_Y = MakeTreeBranch(outTree,'shower1_start_2d_U_Y','int')
        self._shower1_start_2d_U_Z = MakeTreeBranch(outTree,'shower1_start_2d_U_Z','int')
        self._shower1_start_2d_V_X = MakeTreeBranch(outTree,'shower1_start_2d_V_X','int')
        self._shower1_start_2d_V_Y = MakeTreeBranch(outTree,'shower1_start_2d_V_Y','int')
        self._shower1_start_2d_V_Z = MakeTreeBranch(outTree,'shower1_start_2d_V_Z','int')
        self._shower1_start_2d_Y_X = MakeTreeBranch(outTree,'shower1_start_2d_Y_X','int')
        self._shower1_start_2d_Y_Y = MakeTreeBranch(outTree,'shower1_start_2d_Y_Y','int')
        self._shower1_start_2d_Y_Z = MakeTreeBranch(outTree,'shower1_start_2d_Y_Z','int')
        self._shower1_impact = MakeTreeBranch(outTree,'_shower1_impact','float')
        self._shower1_smallq_U = MakeTreeBranch(outTree,'shower1_smallQ_U','float')
        self._shower1_smallq_V = MakeTreeBranch(outTree,'shower1_smallQ_V','float')
        self._shower1_smallq_Y = MakeTreeBranch(outTree,'shower1_smallQ_Y','float')
        self._shower1_sumq_U = MakeTreeBranch(outTree,'shower1_sumQ_U','float')
        self._shower1_sumq_V = MakeTreeBranch(outTree,'shower1_sumQ_V','float')
        self._shower1_sumq_Y = MakeTreeBranch(outTree,'shower1_sumQ_Y','float')

        self._shower2_E_U = MakeTreeBranch(outTree,'shower2_E_U','float')
        self._shower2_E_V = MakeTreeBranch(outTree,'shower2_E_V','float')
        self._shower2_E_Y = MakeTreeBranch(outTree,'shower2_E_Y','float')
        self._shower2_gap_U = MakeTreeBranch(outTree,'shower2_gap_U','int')
        self._shower2_gap_V = MakeTreeBranch(outTree,'shower2_gap_V','int')
        self._shower2_gap_Y = MakeTreeBranch(outTree,'shower2_gap_Y','int')
        self._shower2_dir_3d_X = MakeTreeBranch(outTree,'shower2_dir_3d_X','float')
        self._shower2_dir_3d_Y = MakeTreeBranch(outTree,'shower2_dir_3d_Y','float')
        self._shower2_dir_3d_Z = MakeTreeBranch(outTree,'shower2_dir_3d_Z','float')
        self._shower2_dir_2d_U = MakeTreeBranch(outTree,'shower2_dir_2d_U','float')
        self._shower2_dir_2d_V = MakeTreeBranch(outTree,'shower2_dir_2d_V','float')
        self._shower2_dir_2d_Y = MakeTreeBranch(outTree,'shower2_dir_2d_Y','float')
        self._shower2_op_2d_U = MakeTreeBranch(outTree,'shower2_op_2d_U','float')
        self._shower2_op_2d_V = MakeTreeBranch(outTree,'shower2_op_2d_V','float')
        self._shower2_op_2d_Y = MakeTreeBranch(outTree,'shower2_op_2d_Y','float')
        self._shower2_start_2d_U_X = MakeTreeBranch(outTree,'shower2_start_2d_U_X','int')
        self._shower2_start_2d_U_Y = MakeTreeBranch(outTree,'shower2_start_2d_U_Y','int')
        self._shower2_start_2d_U_Z = MakeTreeBranch(outTree,'shower2_start_2d_U_Z','int')
        self._shower2_start_2d_V_X = MakeTreeBranch(outTree,'shower2_start_2d_V_X','int')
        self._shower2_start_2d_V_Y = MakeTreeBranch(outTree,'shower2_start_2d_V_Y','int')
        self._shower2_start_2d_V_Z = MakeTreeBranch(outTree,'shower2_start_2d_V_Z','int')
        self._shower2_start_2d_Y_X = MakeTreeBranch(outTree,'shower2_start_2d_Y_X','int')
        self._shower2_start_2d_Y_Y = MakeTreeBranch(outTree,'shower2_start_2d_Y_Y','int')
        self._shower2_start_2d_Y_Z = MakeTreeBranch(outTree,'shower2_start_2d_Y_Z','int')
        self._shower2_impact = MakeTreeBranch(outTree,'_shower2_impact','float')
        self._shower_alpha = MakeTreeBranch(outTree,'_shower_alpha','float')
        self._pi0mass = MakeTreeBranch(outTree,'_pi0mass','float')
        self._shower2_smallq_U = MakeTreeBranch(outTree,'shower2_smallQ_U','float')
        self._shower2_smallq_V = MakeTreeBranch(outTree,'shower2_smallQ_V','float')
        self._shower2_smallq_Y = MakeTreeBranch(outTree,'shower2_smallQ_Y','float')
        self._shower2_sumq_U = MakeTreeBranch(outTree,'shower2_sumQ_U','float')
        self._shower2_sumq_V = MakeTreeBranch(outTree,'shower2_sumQ_V','float')
        self._shower2_sumq_Y = MakeTreeBranch(outTree,'shower2_sumQ_Y','float')

        # CNN Shower variables
        self._cnn_shower_energy = MakeTreeBranch(outTree,'cnn_shower_energy','float')

        # mc shower Variables
        self._haspi0 = MakeTreeBranch(outTree,'haspi0','int')
        self._ccnc = MakeTreeBranch(outTree,'ccnc','int')
        self._truefid = MakeTreeBranch(outTree,'truefid','int')
        self._numshowers = MakeTreeBranch(outTree,'numshowers_true','int')
        self._seconddirection_true_X = MakeTreeBranch(outTree,'second_direction_true_X','float')
        self._seconddirection_true_Y = MakeTreeBranch(outTree,'second_direction_true_Y','float')
        self._seconddirection_true_Z = MakeTreeBranch(outTree,'second_direction_true_Z','float')
        self._firstdirection_true_X = MakeTreeBranch(outTree,'first_direction_true_X','float')
        self._firstdirection_true_Y = MakeTreeBranch(outTree,'first_direction_true_Y','float')
        self._firstdirection_true_Z = MakeTreeBranch(outTree,'first_direction_true_Z','float')
        self._shower_start_2d_true_row = MakeTreeBranch(outTree,'shower_start_2d_true_row','int')
        self._shower_start_2d_true_colU = MakeTreeBranch(outTree,'shower_start_2d_true_colU','int')
        self._shower_start_2d_true_colV = MakeTreeBranch(outTree,'shower_start_2d_true_colV','int')
        self._shower_start_2d_true_colY = MakeTreeBranch(outTree,'shower_start_2d_true_colY','int')
        self._secondshower_start_2d_true_row = MakeTreeBranch(outTree,'second_shower_start_2d_true_row','int')
        self._secondshower_start_2d_true_colU = MakeTreeBranch(outTree,'second_shower_start_2d_true_colU','int')
        self._secondshower_start_2d_true_colV = MakeTreeBranch(outTree,'second_shower_start_2d_true_colV','int')
        self._secondshower_start_2d_true_colY = MakeTreeBranch(outTree,'second_shower_start_2d_true_colY','int')
        self._shower_energy_true = MakeTreeBranch(outTree,'shower_energy_true','float')
        self._secondshower_energy_true = MakeTreeBranch(outTree,'secondshower_energy_true','float')
        self._shower_recotrue_dist = MakeTreeBranch(outTree,'shower_recotrue_dist','float')
        self._secondshower_recotrue_dist = MakeTreeBranch(outTree,'secondshower_recotrue_dist','float')


        # Precut stuff
        self._totPE = MakeTreeBranch(outTree,'TotPE','float')
        self._porchTotPE = MakeTreeBranch(outTree,'PorchTotPE','float')
        self._maxPEFrac = MakeTreeBranch(outTree,'MaxPEFrac','float')
        self._passPMTPrecut = MakeTreeBranch(outTree,'PassPMTPrecut','int')
        self._precutBeamFirstTick = MakeTreeBranch(outTree,'PrecutBeamFirstTick','int')
        self._precutVetoFirstTick = MakeTreeBranch(outTree,'PrecutVetoFirstTick','int')


        # BDT output
        self._bdtscore_1e1p         = MakeTreeBranch(outTree,'BDTscore_1e1p','float')
        self._bdtscore_1mu1p_cosmic = MakeTreeBranch(outTree,'BDTscore_1mu1p_cosmic','float')
        self._bdtscore_1mu1p_nu     = MakeTreeBranch(outTree,'BDTscore_1mu1p_nu','float')

        # MC stuff
        self._parentPDG = MakeTreeBranch(outTree,'MC_parentPDG','int')
        self._energyInit = MakeTreeBranch(outTree,'MC_energyInit','float')
        self._parentX = MakeTreeBranch(outTree,'MC_parentX','float')
        self._parentY = MakeTreeBranch(outTree,'MC_parentY','float')
        self._parentZ = MakeTreeBranch(outTree,'MC_parentZ','float')
        self._nproton = MakeTreeBranch(outTree,'MC_nproton','int')
        self._nlepton = MakeTreeBranch(outTree,'MC_nlepton','int')
        self._parentSCEX = MakeTreeBranch(outTree,'MC_parentSCEX','float')
        self._parentSCEY = MakeTreeBranch(outTree,'MC_parentSCEY','float')
        self._parentSCEZ = MakeTreeBranch(outTree,'MC_parentSCEZ','float')
        self._scedr = MakeTreeBranch(outTree,'MC_scedr','float')

        # MPID stuff
        self._eminusPID_int_v = MakeTreeBranch(outTree,'EminusPID_int_v','tvector')
        self._muonPID_int_v = MakeTreeBranch(outTree,'MuonPID_int_v','tvector')
        self._protonPID_int_v = MakeTreeBranch(outTree,'ProtonPID_int_v','tvector')
        self._gammaPID_int_v = MakeTreeBranch(outTree,'GammaPID_int_v','tvector')
        self._pionPID_int_v = MakeTreeBranch(outTree,'PionPID_int_v','tvector')
        self._eminusPID_pix_v = MakeTreeBranch(outTree,'EminusPID_pix_v','tvector')
        self._muonPID_pix_v = MakeTreeBranch(outTree,'MuonPID_pix_v','tvector')
        self._protonPID_pix_v = MakeTreeBranch(outTree,'ProtonPID_pix_v','tvector')
        self._gammaPID_pix_v = MakeTreeBranch(outTree,'GammaPID_pix_v','tvector')
        self._pionPID_pix_v = MakeTreeBranch(outTree,'PionPID_pix_v','tvector')

        self.outTree = outTree

    def clear_vertex(self):
        """ clear values. also define when run for first time """
        self._lepton_id        = int(-9999)
        self._lepton_phi       = float(-9999)
        self._lepton_theta     = float(-9999)
        self._lepton_length    = float(-9999)
        self._lepton_dqdx_uncalibrated      = float(-9999)
        self._lepton_dqdx_calibrated      = float(-9999)
        self._lepton_edge_dist = float(-9999)
        self._muon_E         = float(-9999)
        self._muon_phiB_1m1p      = float(-9999)
        self._muon_thetaB_1m1p    = float(-9999)
        self._muon_EB_1m1p        = float(-9999)
        self._electron_E          = float(-9999)
        self._electron_phiB_1e1p      = float(-9999)
        self._electron_thetaB_1e1p    = float(-9999)
        self._electron_EB_1e1p        = float(-9999)

        self._proton_id        = int(-9999)
        self._proton_phi       = float(-9999)
        self._proton_theta     = float(-9999)
        self._proton_length    = float(-9999)
        self._proton_dqdx_uncalibrated      = float(-9999)
        self._proton_dqdx_calibrated      = float(-9999)
        self._proton_E         = float(-9999)
        self._proton_edge_dist = float(-9999)
        self._proton_phiB_1m1p      = float(-9999)
        self._proton_thetaB_1m1p    = float(-9999)
        self._proton_EB_1m1p        = float(-9999)
        self._proton_phiB_1e1p      = float(-9999)
        self._proton_thetaB_1e1p    = float(-9999)
        self._proton_EB_1e1p        = float(-9999)

        self._shower1_E_U = float(-9999)
        self._shower1_E_V = float(-9999)
        self._shower1_E_Y = float(-9999)
        self._shower2_E_U = float(-9999)
        self._shower2_E_V = float(-9999)
        self._shower2_E_Y = float(-9999)
        self._shower1_gap_U = int(-9999)
        self._shower1_gap_V = int(-9999)
        self._shower1_gap_Y = int(-9999)
        self._shower1_dir_3d_X = float(-9999)
        self._shower1_dir_3d_Y  = float(-9999)
        self._shower1_dir_3d_Z = float(-9999)
        self._shower1_dir_2d_U = float(-9999)
        self._shower1_dir_2d_V = float(-9999)
        self._shower1_dir_2d_Y = float(-9999)
        self._shower1_op_2d_U = float(-9999)
        self._shower1_op_2d_V = float(-9999)
        self._shower1_op_2d_Y = float(-9999)
        self._shower1_start_2d_U_X = int(-9999)
        self._shower1_start_2d_U_Y = int(-9999)
        self._shower1_start_2d_U_Z = int(-9999)
        self._shower1_start_2d_V_X  = int(-9999)
        self._shower1_start_2d_V_Y = int(-9999)
        self._shower1_start_2d_V_Z  = int(-9999)
        self._shower1_start_2d_Y_X = int(-9999)
        self._shower1_start_2d_Y_Y = int(-9999)
        self._shower1_start_2d_Y_Z = int(-9999)
        self._shower1_impact = float(-9999)
        self._shower1_smallq_U = float(-9999)
        self._shower1_smallq_V = float(-9999)
        self._shower1_smallq_Y = float(-9999)
        self._shower2_gap_U = int(-9999)
        self._shower2_gap_V = int(-9999)
        self._shower2_gap_Y = int(-9999)
        self._shower2_dir_3d_X = float(-9999)
        self._shower2_dir_3d_Y = float(-9999)
        self._shower2_dir_3d_Z  = float(-9999)
        self._shower2_dir_2d_U = float(-9999)
        self._shower2_dir_2d_V = float(-9999)
        self._shower2_dir_2d_Y = float(-9999)
        self._shower2_op_2d_U  = float(-9999)
        self._shower2_op_2d_V = float(-9999)
        self._shower2_op_2d_Y = float(-9999)
        self._shower2_start_2d_U_X = int(-9999)
        self._shower2_start_2d_U_Y = int(-9999)
        self._shower2_start_2d_U_Z = int(-9999)
        self._shower2_start_2d_V_X = int(-9999)
        self._shower2_start_2d_V_Y = int(-9999)
        self._shower2_start_2d_V_Z = int(-9999)
        self._shower2_start_2d_Y_X = int(-9999)
        self._shower2_start_2d_Y_Y = int(-9999)
        self._shower2_start_2d_Y_Z = int(-9999)
        self._shower2_impact = float(-9999)
        self._shower2_smallq_U = float(-9999)
        self._shower2_smallq_V = float(-9999)
        self._shower2_smallq_Y = float(-9999)
        self._shower_alpha = float(-9999)
        self._pi0mass = float(-9999)
        self._haspi0= int(-9999)
        self._ccnc= int(-9999)
        self._truefid= int(-9999)
        self._numshowers= int(-9999)
        self._seconddirection_true_X= float(-9999)
        self._seconddirection_true_Y= float(-9999)
        self._seconddirection_true_Z= float(-9999)
        self._firstdirection_true_X= float(-9999)
        self._firstdirection_true_Y= float(-9999)
        self._firstdirection_true_Z= float(-9999)
        self._shower_start_2d_true_row= int(-9999)
        self._shower_start_2d_true_colU= int(-9999)
        self._shower_start_2d_true_colV= int(-9999)
        self._shower_start_2d_true_colY= int(-9999)
        self._secondshower_start_2d_true_row= int(-9999)
        self._secondshower_start_2d_true_colU= int(-9999)
        self._secondshower_start_2d_true_colV= int(-9999)
        self._secondshower_start_2d_true_colY= int(-9999)
        self._shower_energy_true= float(-9999)
        self._secondshower_energy_true= float(-9999)
        self._shower_recotrue_dist= float(-9999)
        self._secondshower_recotrue_dist= float(-9999)



        return