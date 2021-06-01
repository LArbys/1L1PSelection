import os,sys
from event_indexed_trees_util import is_tree_event_indexed

def get_tree_lists(input_file,IS_MC=False):
    """
    Get list of trees in the input (DLANA) file. Split them up into certain categories
    * vertex-indexed
    * event-indexed
    * other trees
    * pot summary tree
    * THE final vertex variable tree
    """
    vertex_indexed_trees = []
    event_indexed_trees  = []
    other_trees = []
    pot_sum_tree = None
    fvv_tree = None
    
    

    dirlist = [None]
    dirdict = {}
        
    while len(dirlist)>0:
        dirname = dirlist.pop(0)
        print "searching dir: ",dirname
        tdir = input_file
        if dirname is not None:
            tdir = input_file.Get(dirname)
        tdir.cd()
        tree_keys = tdir.GetListOfKeys()
        nkeys = tree_keys.GetEntries()
        for ikey in xrange(nkeys):
            basetreename = tree_keys.At(ikey).GetName()
            if dirname is not None:
                treename = dirname+"/"+basetreename
                dirdict[basetreename] = dirname
            else:
                treename = basetreename            
                        
            atree = input_file.Get(treename)
            if treename=="dlana/FinalVertexVariables":
                fvv_tree = atree

            if atree is not None and atree.ClassName()=="TTree":
                nentries = atree.GetEntries()
                print "Tree: ",treename," ",atree.ClassName()," nentries=",nentries
                if treename=="potsummary_generator_tree":
                    print "Found potsummary_generator_tree. special case."
                    pot_sum_tree = atree
                    continue

                if is_tree_event_indexed( basetreename, is_mc=IS_MC ):
                    event_indexed_trees.append(atree)
                else:
                    vertex_indexed_trees.append(atree)

            elif atree is not None and atree.ClassName()=="TDirectoryFile":
                dirlist.append(treename)
            else:
                print "unrecognized: ",atree.ClassName()

        print "directories remaining: ",len(dirlist),dirlist

    return vertex_indexed_trees,event_indexed_trees,other_trees,fvv_tree,pot_sum_tree,dirdict
