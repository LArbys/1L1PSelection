import ROOT
import numpy as np
from larlite import larlite,larutil


def PerformPMTPrecuts(s_dlmerged,
                      ophittree = "ophitBeam",
                      PP_WINDOW_LENGTH = 130,
                      PP_COINC_WINDOW = 6,
                      PP_PE_THRESH = 20.0,
                      PP_PMT_MAX_FRAC_CUTOFF = 0.60,
                      PP_WIN_START = 190,
                      PP_WIN_END = 320,
                      PP_PORCH_WIN_START = 60,
                      PP_PORCH_WIN_END =190,
                      PP_TICK_SIZE = 0.015625 ):
    """
    calculate PMT precuts
    TICK_SIZE units in 15.625 ns per tick
    """

    # make config
    from ROOT import fcllite
    fcllite.PSet
    fcllite.CreatePSetFromFile

    cfg = """OpHitProducer:\"%s\"
         BinTickWidth:%d
         WinStartTick:%d
         WinEndTick:%d
         PEThreshold:%.2f
         VetoPEThreshold:%.2f
         MaxVetoPE:%.2f
         VetoStartTick:%.2f
         VetoEndTick:%.2f
         PMTMaxFrac:%.2f
        """%(ophittree,
            PP_COINC_WINDOW,
            PP_WIN_START,
            PP_WIN_END,
            PP_PE_THRESH,
            PP_PE_THRESH,
            PP_PE_THRESH,
            PP_PORCH_WIN_START,
            PP_PORCH_WIN_END,
            PP_PMT_MAX_FRAC_CUTOFF)

    fout = open("precut.cfg",'w')
    print>>fout,cfg
    fout.close()

    pset = fcllite.CreatePSetFromFile("precut.cfg","PMTPreCut")
    #print pset.dump()

    pset = fcllite.PSet("PMTPreCut",cfg)
    precutalgo = larlite.LEEPreCut()
    precutalgo.configure(pset)

    print "PMT Precut algo configured"

    # ---------------------------------------------- #
    #print('<EVID: %s> -- First, we will figure out the PMT Precut info.'%_tag)  #gotta do this first for io reasons
    PMTPrecutDict = {}

    # Load up larlite
    ll_manager = larlite.storage_manager()
    ll_manager.set_io_mode(ll_manager.kREAD)
    ll_manager.add_in_filename(s_dlmerged)
    ll_manager.set_in_rootdir("")
    ll_manager.open()

    while ll_manager.next_event():

        id_rse = tuple((ll_manager.run_id(),ll_manager.subrun_id(),ll_manager.event_id()))

        ev_ophits   = ll_manager.get_data(larlite.data.kOpHit,ophittree)
        print "[event {}] number of ophits: {}".format( id_rse, ev_ophits.size() )
        passcuts = precutalgo.apply( ev_ophits )
        PMTPrecutDict[id_rse] = dict(_totpe=precutalgo.beamPE(),
                                     _porchtotpe=precutalgo.vetoPE(),
                                     _maxpefrac=precutalgo.maxFrac(),
                                     _passpmtprecut=passcuts,
                                     _beamFirstTick=precutalgo.beamFirstTick(),
                                     _vetoFirstTick=precutalgo.vetoFirstTick() )

    return PMTPrecutDict

def GetPMTPrecutDict(s_dlmerged):
    """ get results from previously applied pmt precuts """
    PMTPrecutDict = {}

    # Load up larlite
    ll_manager = larlite.storage_manager()
    ll_manager.set_io_mode(ll_manager.kREAD)
    ll_manager.add_in_filename(s_dlmerged)
    ll_manager.set_data_to_read(larlite.data.kUserInfo,'precutresults')
    ll_manager.set_in_rootdir("")
    ll_manager.open()

    while ll_manager.next_event():
        id_rse = tuple((ll_manager.run_id(),ll_manager.subrun_id(),ll_manager.event_id()))

        # ---------- Grab vectors containing particle info --------- #
        precutresults   = ll_manager.get_data(larlite.data.kUserInfo,"precutresults")
        precutresult = precutresults[0]
        # ---------------------------------------------------------- #
        PMTPrecutDict[id_rse] = dict(_totpe=precutresult.get_double("beamPE"),
                                     _porchtotpe=precutresult.get_double("vetoPE"),
                                     _maxpefrac=precutresult.get_double("maxFrac"),
                                     _passpmtprecut=precutresult.get_int("pass"))

    ll_manager.close()
    return PMTPrecutDict
