import os,sys
import torch
from collections import OrderedDict
import ShowerEstimatorCNN

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
