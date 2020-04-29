import os,sys
import numpy as np
import torch
from collections import OrderedDict
import ShowerEstimatorCNN
from larlitecv import larlitecv

def load_showercnn_model( weightfile_name ):
    """ model is loaded from a pickle """
    weight_path = os.environ["UBDLANA_DIR"]+"/showercnn_models/"+weightfile_name
    state_dict = torch.load(weight_path, map_location="cpu")
    new_state_dict = OrderedDict()
    for k,v in state_dict.items():
        name = k[7:] #removing 'module.' for each weight tensor
        new_state_dict[name] = v

    model = ShowerEstimatorCNN.resnet18(in_channels = 6, n_classes = 1)
    model.load_state_dict( new_state_dict )
    model.eval()

    return model

def run_showercnn_event( iolcv, model  ):
    shrutil = larlitecv.ssnetshowerreco.ShowerRecoUtil()
    vtx_img_v = shrutil.process_event_get_ndarray( iolcv )
    vtx_shr_energy_mev_v = []
    for vtx_id, img_dict in enumerate(vtx_img_v):

        if len(img_dict["adc"])!=3:
            print "[showercnnutil::run_showercnn_event] no image crops for vtxid[",vtx_id,"]"
            vtx_shr_energy_mev_v.append(-9999.0)
            continue
        
        img_arrU = np.log(img_dict["adc"][0] + 1)
        img_arrV = np.log(img_dict["adc"][1] + 1)
        img_arrY = np.log(img_dict["adc"][2] + 1)
        chs_arrU = img_dict["status"][0]
        chs_arrV = img_dict["status"][1]
        chs_arrY = img_dict["status"][2]


        x = np.array((img_arrU,img_arrV,img_arrY,chs_arrU,chs_arrV,chs_arrY))
        x = x.reshape( (1,x.shape[0],x.shape[1],x.shape[2]) )
        print "[showercnnutil::run_showercnn_event] crop for vtxid[",vtx_id,"] prepared. input.shape=",x.shape
        
        with torch.set_grad_enabled(False):
            x = torch.from_numpy(x).float().to("cpu")
            outputs = model(x).cpu().detach().numpy()
            vtx_shr_energy_mev_v.append(outputs[0][0]*2000.0)
            print "[showercnnutil::run_showercnn_event] result for vtxid[",vtx_id,"] ",vtx_shr_energy_mev_v[-1]," MeV"
    return vtx_shr_energy_mev_v
