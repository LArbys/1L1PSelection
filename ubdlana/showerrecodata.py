import os,sys

class ShowerRecoData:
    """ provide interface into showerreco json structure """
    def __init__( self, json ):
        self.json = json
        # index the dictionary
        self.index = {}
        for index,entry in enumerate(json["entries"]):
            self.index[ (entry["run"],entry["subrun"],entry["event"]) ] = index
    def get_entry(self,run,subrun,event):
        if (run,subrun,event) not in self.index:
            return None
        return self.json["entries"][self.index[(run,subrun,event)]]
    def get_vertex(self,run,subrun,event,vtxid):
        entrydata = self.get_entry( run, subrun, event )
        if entrydata is None:
            return None

        vertexinfo = {}
        if vtxid<0 or vtxid>=len( entrydata["shower_energies"] ):
            return None

        for k,data in entrydata.items():
            if k in ["run","subrun","entry","event"]:
                continue
            #print type(data),data
            vertexinfo[k] = data[vtxid]

        return vertexinfo
