import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
from numpy import mean
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from numpy import asarray
from root_pandas import read_root
from matplotlib import gridspec
from scipy import stats,signal
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from textwrap import wrap
import seaborn as sns

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os


class distVar:

    def __init__(self,s_name,n_range,s_label='',s_cov=''):
        self.myname = s_name
        self.myrange = n_range
        self.myscov = s_cov
        if s_label == '':
            self.mylabel = s_name
        else:
            self.mylabel = s_label

class FullCov:
    def __init__(self,s_detsyspickle):
        
        with open(s_detsyspickle,'rb') as handle: (self.a_overlap_sys_numu,self.a_cv_sys_numu,self.a_overlap_sys_nue,self.a_cv_sys_nue,self.s_detsyslist) = pickle.load(handle)
        self.SetDefaultCuts()
        self.a_detSysEnabled = np.ones(len(self.s_detsyslist))
    
    def SetDefaultCuts(self):     # the current selection
        self.s_cuts_numu = 'bkgBDT_run3 < .4 and Enu_1m1p < 1000 and Proton_CosTheta > 0'
        self.s_cuts_cv_numu = 'bkgBDT_run3_cv < .4 and Enu_1m1p_cv < 1000 and Proton_CosTheta_cv > 0'
        self.s_cuts_nue = 'PassFinalSelection1e1p >= 0'
        self.s_cuts_cv_nue = 'PassFinalSelection1e1p_cv >= 0'
    
    def SetCuts(self,s_cuts_numu,s_cuts_cv_numu,s_cuts_nue,s_cuts_cv_nue):
        self.s_cuts_numu = s_cuts_numu
        self.s_cuts_cv_numu = s_cuts_cv_numu
        self.s_cuts_nue = s_cuts_nue
        self.s_cuts_cv_nue = s_cuts_cv_nue
        
    def ListVariables(self):
        varlist = list(self.a_overlap_sys_numu[0])
        s_varlist = ''
        for s_var in varlist:
            if s_var[-3:]=='_cv':
                continue
            else:
                s_varlist+=s_var+'\n'

        print(s_varlist)

    def ListDetectorSystematics(self):
        print('Detector Systematics:',self.s_detsyslist)
        print('Enabled:',self.a_detSysEnabled)
        
    def EnableDetectorSystematics(self,a_enabled):  
        # here, we can turn different systematics on or off with an array of ones and zeros
        # for each systematic, 1 - on, 0 - off. the default is all ones upon startup.
        
        
        # make sure it's the right length
        if len(a_enabled) != len(self.a_detSysEnabled):
            print('Failed. Array must be of length %i'%len(self.a_detSysEnabled))
            return
            
        # make sure it's just ones and zeros
        for val in a_enabled:
            if val!=0 and val!=1:
                print('Failed. Only ones and zeros please. PLEASE')
                return
        
        # can't just turn everything off...
        if np.asarray(a_enabled).sum() == 0:
            print('Failed. You cannot just turn everything off.')
            return
        
        self.a_detSysEnabled = a_enabled.copy()
        for i in range(len(self.a_detSysEnabled)):
            if self.a_detSysEnabled[i]==1:
                print('%s is enabled'%(self.s_detsyslist[i]))
        print('Enabled:',self.a_detSysEnabled)
        
    def Nominal(self,dvar_numu,dvar_nue,a_nbins,draw=False):

        nbins_numu = a_nbins[0]
        nbins_nue = a_nbins[1]
        
        cov_tru = np.zeros((nbins_numu+nbins_nue,nbins_numu+nbins_nue))
        for sysi in range(len(self.a_overlap_sys_numu)):
            if self.a_detSysEnabled[sysi]==0:
                continue
                
            myvardf_numu = self.a_overlap_sys_numu[sysi].query(self.s_cuts_numu)
            myvarcv_numu = self.a_cv_sys_numu[sysi].query(self.s_cuts_cv_numu)
            myvardf_nue = self.a_overlap_sys_nue[sysi].query(self.s_cuts_nue)
            myvarcv_nue = self.a_cv_sys_nue[sysi].query(self.s_cuts_cv_nue)
            
            var_sys_numu = myvardf_numu[dvar_numu.myname]
            var_cv_numu = myvarcv_numu[dvar_numu.myname+'_cv']
            var_sys_nue = myvardf_nue[dvar_nue.myname]
            var_cv_nue = myvarcv_nue[dvar_nue.myname+'_cv']
            
            hCV_numu,binedges_numu = np.histogram(var_cv_numu,bins=nbins_numu,range=dvar_numu.myrange)
            h0_numu,_ = np.histogram(var_sys_numu,bins=nbins_numu,range=dvar_numu.myrange)
            hCV_nue,binedges_nue = np.histogram(var_cv_nue,bins=nbins_nue,range=dvar_nue.myrange)
            h0_nue,_ = np.histogram(var_sys_nue,bins=nbins_nue,range=dvar_nue.myrange)

            hCV = np.concatenate((hCV_numu,hCV_nue))
            h0 = np.concatenate((h0_numu,h0_nue))
            
            for i in range(nbins_numu+nbins_nue):
                for j in range(nbins_numu+nbins_nue):
                    cov_tru[i][j] += (h0[i]-hCV[i])*(h0[j]-hCV[j])/(hCV[i]*hCV[j])
    
        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(range(nbins_nue+nbins_numu+1),range(nbins_numu+nbins_nue+1))
            crat_tru = ax.pcolormesh(X, Y,cov_tru,cmap='cool')
            cbar = fig.colorbar(crat_tru)
            ax.set_xlabel(dvar_numu.mylabel+' '+dvar_nue.mylabel,fontsize=20)
            ax.set_ylabel(dvar_numu.mylabel+' '+dvar_nue.mylabel,fontsize=20)
            ax.axvline(nbins_numu,color='red')
            ax.axhline(nbins_numu,color='red')
    
        return cov_tru
        
    def Flat(self,dvar_numu,dvar_nue,a_nbins,draw=False):
        
        nbins_numu = a_nbins[0]
        nbins_nue = a_nbins[1]
        
        flatsys = 0
        for sysi in range(len(self.a_overlap_sys_numu)):
            if self.a_detSysEnabled[sysi]==0:
                continue
                
            myvardf_numu = self.a_overlap_sys_numu[sysi].query(self.s_cuts_numu)
            myvarcv_numu = self.a_cv_sys_numu[sysi].query(self.s_cuts_cv_numu)
            myvardf_nue = self.a_overlap_sys_nue[sysi].query(self.s_cuts_nue)
            myvarcv_nue = self.a_cv_sys_nue[sysi].query(self.s_cuts_cv_nue)
            
            var_sys_numu = myvardf_numu[dvar_numu.myname]
            var_cv_numu = myvarcv_numu[dvar_numu.myname+'_cv']
            var_sys_nue = myvardf_nue[dvar_nue.myname]
            var_cv_nue = myvarcv_nue[dvar_nue.myname+'_cv']
            
            hCV_numu,binedges_numu = np.histogram(var_cv_numu,bins=nbins_numu,range=dvar_numu.myrange)
            h0_numu,_ = np.histogram(var_sys_numu,bins=nbins_numu,range=dvar_numu.myrange)
            hCV_nue,binedges_nue = np.histogram(var_cv_nue,bins=nbins_nue,range=dvar_nue.myrange)
            h0_nue,_ = np.histogram(var_sys_nue,bins=nbins_nue,range=dvar_nue.myrange)

            hCV = np.concatenate((hCV_numu,hCV_nue))
            h0 = np.concatenate((h0_numu,h0_nue))
        
            flatsys += (h0.sum()-hCV.sum())*(h0.sum()-hCV.sum())/(hCV.sum()*hCV.sum())
        
        if(draw):
            print('Flatsys:',flatsys,' - (',np.sqrt(flatsys),' frac error)')        
        return np.ones((nbins_numu+nbins_nue,nbins_numu+nbins_nue))*flatsys
    
    def Polyfit(self,dvar_numu,dvar_nue,a_nbins,draw=False):
        
        nbins_numu = a_nbins[0]
        nbins_nue = a_nbins[1]
        
        cov_poly = np.zeros((nbins_numu+nbins_nue,nbins_numu+nbins_nue))
        for sysi in range(len(self.a_overlap_sys_numu)):
            if self.a_detSysEnabled[sysi]==0:
                continue
            
            myvardf_numu = self.a_overlap_sys_numu[sysi].query(self.s_cuts_numu)
            myvarcv_numu = self.a_cv_sys_numu[sysi].query(self.s_cuts_cv_numu)
            myvardf_nue = self.a_overlap_sys_nue[sysi].query(self.s_cuts_nue)
            myvarcv_nue = self.a_cv_sys_nue[sysi].query(self.s_cuts_cv_nue)
            
            var_sys_numu = myvardf_numu[dvar_numu.myname]
            var_cv_numu = myvarcv_numu[dvar_numu.myname+'_cv']
            var_sys_nue = myvardf_nue[dvar_nue.myname]
            var_cv_nue = myvarcv_nue[dvar_nue.myname+'_cv']
            
            hCV_numu,binedges_numu = np.histogram(var_cv_numu,bins=nbins_numu,range=dvar_numu.myrange)
            h0_numu,_ = np.histogram(var_sys_numu,bins=nbins_numu,range=dvar_numu.myrange)
            hCV_nue,binedges_nue = np.histogram(var_cv_nue,bins=nbins_nue,range=dvar_nue.myrange)
            h0_nue,_ = np.histogram(var_sys_nue,bins=nbins_nue,range=dvar_nue.myrange)
            bincenters_numu = np.diff(binedges_numu)/2 + binedges_numu[:-1]    
            bincenters_nue = np.diff(binedges_nue)/2 + binedges_nue[:-1]    
            
            
            # get polyfit degs for numu ----
            truRat_numu = np.true_divide(h0_numu,hCV_numu,out=np.ones_like(bincenters_numu),where=hCV_numu!=0)
            aics = []
            degs = []
            for deg in range(min(nbins_numu-2,int(nbins_numu/2))):
                params = deg + 1
                polyRat = np.polyfit(bincenters_numu, truRat_numu, deg)
                fRat = np.poly1d(polyRat)
                
                # now calculate chi2 for fit
                yerr_rat = np.true_divide(np.sqrt(fRat(bincenters_numu)*hCV_numu),hCV_numu,out=np.zeros_like(bincenters_numu),where=hCV_numu!=0)
                chi2_fit = np.power(np.true_divide(fRat(bincenters_numu)-truRat_numu,yerr_rat),2).sum()
                aic = chi2_fit + 2*params + 2*params*(params+1)/float(nbins_numu-params-1)
                aics.append(aic)
                degs.append(deg)
            polyterms = degs[np.argmin(aics)]
            print(self.s_detsyslist[sysi],'Numu Polyfit Degrees:',polyterms,aics[np.argmin(aics)])  

            polyRat_numu = np.polyfit(bincenters_numu, np.true_divide(h0_numu,hCV_numu,where=hCV_numu!=0), polyterms)
            fRat_numu = np.poly1d(polyRat_numu) 
            h0_fit_numu = fRat_numu(bincenters_numu)*hCV_numu
            
            # get polyfit degs for nue ----
            truRat_nue = np.true_divide(h0_nue,hCV_nue,out=np.ones_like(bincenters_nue),where=hCV_nue!=0)
            aics = []
            degs = []
            for deg in range(min(nbins_nue-2,int(nbins_nue/2))):
                params = deg + 1
                polyRat = np.polyfit(bincenters_nue, truRat_nue, deg)
                fRat = np.poly1d(polyRat)
                
                # now calculate chi2 for fit
                yerr_rat = np.true_divide(np.sqrt(fRat(bincenters_nue)*hCV_nue),hCV_nue,out=np.zeros_like(bincenters_nue),where=hCV_nue!=0)
                chi2_fit = np.power(np.true_divide(fRat(bincenters_nue)-truRat_nue,yerr_rat),2).sum()
                aic = chi2_fit + 2*params + 2*params*(params+1)/float(nbins_nue-params-1)
                aics.append(aic)
                degs.append(deg)
            polyterms = degs[np.argmin(aics)]
            print(self.s_detsyslist[sysi],'Nue Polyfit Degrees:',polyterms,aics[np.argmin(aics)])  

            polyRat = np.polyfit(bincenters_nue, np.true_divide(h0_nue,hCV_nue,where=hCV_nue!=0), polyterms)
            fRat = np.poly1d(polyRat) 
            h0_fit_nue = fRat(bincenters_nue)*hCV_nue
            
            
            hCV = np.concatenate((hCV_numu,hCV_nue))
            h0_fit = np.concatenate((h0_fit_numu,h0_fit_nue))
                 
            for i in range(nbins_nue+nbins_numu):
                for j in range(nbins_nue+nbins_numu):
                    cov_poly[i][j] += (h0_fit[i]-hCV[i])*(h0_fit[j]-hCV[j])/(hCV[i]*hCV[j])
            
        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(range(nbins_nue+nbins_numu+1),range(nbins_numu+nbins_nue+1))
            crat_poly = ax.pcolormesh(X, Y,cov_poly,cmap='cool')
            cbar = fig.colorbar(crat_poly)
            ax.set_xlabel(dvar_numu.mylabel+' '+dvar_nue.mylabel,fontsize=20)
            ax.set_ylabel(dvar_numu.mylabel+' '+dvar_nue.mylabel,fontsize=20)
            ax.axvline(nbins_numu,color='red')
            ax.axhline(nbins_numu,color='red')
    
        return cov_poly
    
    
    def DetectorSystematicDiagnosticsNumu(self,dvar,nbins):
        self.__DetectorSystematicDiagnostics(dvar,self.s_cuts_numu,self.s_cuts_cv_numu,nbins,self.a_overlap_sys_numu,self.a_overlap_sys_numu)
    
    def DetectorSystematicDiagnosticsNue(self,dvar,nbins):
        self.__DetectorSystematicDiagnostics(dvar,self.s_cuts_nue,self.s_cuts_cv_nue,nbins,self.a_overlap_sys_nue,self.a_overlap_sys_nue)
    
    
    def __DetectorSystematicDiagnostics(self,dvar,s_cuts,s_cuts_cv,nbins,a_overlap_sys,a_cv_sys):
    
        flatsys = 0.0

        for sysi in range(len(self.a_detSysEnabled)):
            print(self.s_detsyslist[sysi])
            myvardf = a_overlap_sys[sysi].query(s_cuts)
            myvarcv = a_cv_sys[sysi].query(s_cuts_cv)
        
            var_sys = myvardf[dvar.myname]
            var_cv = myvarcv[dvar.myname+'_cv']
          
            hCV,binedges = np.histogram(var_cv,bins=nbins,range=dvar.myrange)
            h0,_ = np.histogram(var_sys,bins=nbins,range=dvar.myrange)
            bincenters = np.diff(binedges)/2 + binedges[:-1]

            truRat = np.true_divide(h0,hCV,out=np.ones_like(bincenters),where=hCV!=0)
      
            # get polyfit degs
            aics = []
            degs = []
            for deg in range(min(nbins-2,int(nbins/2))):
                params = deg + 1
                polyRat = np.polyfit(bincenters, truRat, deg)
                fRat = np.poly1d(polyRat)
    
                # now calculate chi2 for fit
                yerr_rat = np.true_divide(np.sqrt(fRat(bincenters)*hCV),hCV,out=np.zeros_like(bincenters),where=hCV!=0)
                chi2_fit = np.power(np.true_divide(fRat(bincenters)-truRat,yerr_rat),2).sum()
                aic = chi2_fit + 2*params + 2*params*(params+1)/float(nbins-params-1)
                aics.append(aic)
                degs.append(deg)
      
            polyterms = degs[np.argmin(aics)]
            print('polyfit degrees:',polyterms,aics[np.argmin(aics)])  
            polyRat = np.polyfit(bincenters, np.true_divide(h0,hCV,where=hCV!=0), polyterms)
            fRat = np.poly1d(polyRat) 
            h0_fit = fRat(bincenters)*hCV
    
            fig,ax = plt.subplots(figsize=(27,11))    
            gs = gridspec.GridSpec(1, 2)
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1])
    
            dvarLinspace = np.linspace(dvar.myrange[0],dvar.myrange[1],40)
            ax0.scatter(bincenters,np.true_divide(h0,hCV,where=hCV!=0),label='Variation/CV',s=100)
            ax0.plot(dvarLinspace,fRat(dvarLinspace),label='PolyFit',color='green')
            
            ax0.set_title(self.s_detsyslist[sysi],fontsize=30)
            ax0.set_xlabel(dvar.mylabel,fontsize=20)
            ax0.legend(fontsize=15)
    
            ax1.hist(var_cv,nbins,range=dvar.myrange,histtype='step',linewidth=2,label='CV Raw')
            ax1.hist(var_sys,nbins,range=dvar.myrange,histtype='step',linewidth=3,linestyle='--',label='Variation Raw')
            ax1.plot(bincenters,h0_fit,label='PolyFit Ratio x CV',c='green')
            ax1.legend(fontsize=15)
            ax1.set_xlim(dvar.myrange)
    
            # chi2s
            chisq_nom = 0
            chisq_polyfit = 0
            dof = 0
            for i in range(nbins):    
                if(hCV[i]!=0):
                    dof += 1
                    chisq_nom += np.power(h0[i]-hCV[i],2)/hCV[i]
                    chisq_polyfit += np.power(h0_fit[i]-hCV[i],2)/hCV[i]
            print('Chisq/(%i dof) [CV vs Nominal Variation]:%.3f'%(dof,chisq_nom/float(dof)))
            print('Chisq/(%i dof) [CV vs Polyfit Variation]:%.3f'%(dof,chisq_polyfit/float(dof)))
