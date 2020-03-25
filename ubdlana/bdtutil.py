import os,sys,pickle
import scipy
import xgboost

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

def apply_1e1p_model( dlanavars, model ):
    pass

def apply_1m1p_model( dlanavars, model ):
    pass
