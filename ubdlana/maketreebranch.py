import os,sys
from array import array
import ROOT

# --------------------------------------------- #
# MakeTreeBranch
# --------------------------------------------- #
def MakeTreeBranch(ttree,s_name,s_type):
    """ utility function to setup tree branch """
    if s_type == 'int':
        _myvar = array('i',[-9999])
        ttree.Branch(s_name,_myvar,s_name+'/I')
        return _myvar

    if s_type == 'double':
        _myvar = array('d',[-9999])
        ttree.Branch(s_name,_myvar,s_name+'/D')
	return _myvar

    if s_type == 'float':
        _myvar = array('f',[-9999])
        ttree.Branch(s_name,_myvar,s_name+'/F')
        return _myvar

    if s_type == 'tvector':
        _myvar = ROOT.vector('double')()
        ttree.Branch(s_name,_myvar)
        return _myvar
