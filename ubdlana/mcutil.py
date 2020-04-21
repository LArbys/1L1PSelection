import os,sys

class MCUtil:

    def __init__(self):
        self.mcdata = {}

    def process( self, ioll ):
        ev_mctruth  = ioll.get_data( larlite.data.MCTruth,  "generator" )
        ev_mctrack  = ioll.get_data( larlite.data.MCTrack,  "mcreco" )
        ev_mcshower = ioll.get_data( larlite.data.MCShower, "mcshower" )

        event_data = {"mode":None,
                      "select1L1P0X":None,
                      "Enu":None,
                      "flavor":None,
                      "CCNC":None,
                      "leptonE":None,
                      "protonE":None,
                      "nproton_ke65MeV":None,
                      "npion_ke35MeV":None}

        
