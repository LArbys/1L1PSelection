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

# from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    
    # self.a_detSysEnabled, self.a_detSysName, self.a_OverlapSysDf, self.a_OverlapCVDf
    # self.a_channelSysSCut,  self.a_channelCVSCut, self.a_channelName
    
    def __init__(self,s_detsyspickle):
        
        self.a_OverlapSysDf = []
        self.a_OverlapCVDf = []
        
        self.a_channelSCut = []
        self.a_channelName = []
        self.a_channelSampleID = []
        
        with open(s_detsyspickle,'rb') as handle: (a_overlap_sys_numu,a_cv_sys_numu,a_overlap_sys_nue,a_cv_sys_nue,s_detSysName) = pickle.load(handle)
        
        self.a_OverlapSysDf.append(a_overlap_sys_numu.copy())
        self.a_OverlapSysDf.append(a_overlap_sys_nue.copy())
        self.a_OverlapCVDf.append(a_cv_sys_numu.copy())
        self.a_OverlapCVDf.append(a_cv_sys_nue.copy())
        self.a_detSysName = s_detSysName.copy()
            
        self.a_detSysEnabled = np.ones(len(self.a_detSysName))
        
    def ClearChannels(self):
        self.a_channelSCut = []
        self.a_channelName = []
        self.a_channelSampleID = []
    
    def AddChannel(self,s_cut,i_sample,s_channelName):
        self.a_channelSCut.append(s_cut)
        self.a_channelSampleID.append(i_sample)
        self.a_channelName.append(s_channelName)
    
    def ListChannels(self):
        if len(self.a_channelName)==0:
            print('Need at least one channel, jabroni.')
            return
        else:
            for i in range(len(self.a_channelName)):
                print('Channel ',i,':',self.a_channelName[i])
                print('|',self.a_channelSCut[i],'|')
                print('-----------------------------')
    
    def ListVariables(self):
        varlist = list(self.a_OverlapCVDf[0][0])
        s_varlist = ''
        for s_var in varlist:
            s_varlist+=s_var+'\n'
        print(s_varlist)

    def ListDetectorSystematics(self):
        print('Detector Systematics:',self.a_detSysName)
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
                print('%s is enabled'%(self.a_detSysName[i]))
        print('Enabled:',self.a_detSysEnabled)
    
    
    def Nominal(self,a_dvar,_a_nbins,draw=False):
        a_nbins = np.asarray(_a_nbins)
        if len(self.a_channelName)==0:
            print('Need at least one channel, jabroni.')
            return
        if len(a_dvar) != len(self.a_channelName) or len(a_nbins) != len(self.a_channelName):
            print('variable and bin arrays must match number of channels: ',len(self.a_channelName))
        
        
        cov_tru = np.zeros((a_nbins.sum(),a_nbins.sum()))
        for sysi in range(len(self.a_detSysName)):   # for each systematic
            if self.a_detSysEnabled[sysi]==0:
                continue
            
            hCV = np.zeros(a_nbins.sum())
            hSys = np.zeros(a_nbins.sum())
            
            for chi in range(len(self.a_channelSCut)):  # for each channel
                # get the list of values for the variable assigned to this channel, in the sample assigned to this channel, after applying the cut assigned to this channel
                var_sys = self.a_OverlapSysDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[a_dvar[chi].myname]
                var_cv = self.a_OverlapCVDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[a_dvar[chi].myname]
                    
                h0,_ = np.histogram(var_cv,bins=a_nbins[chi],range=a_dvar[chi].myrange)
                h1,_ = np.histogram(var_sys,bins=a_nbins[chi],range=a_dvar[chi].myrange)
                              
                # fill in the corresponding parts of our grand histogram
                hCV[a_nbins[:chi].sum():a_nbins[:chi+1].sum()] += h0
                hSys[a_nbins[:chi].sum():a_nbins[:chi+1].sum()] += h1
            
            for i in range(a_nbins.sum()):
                for j in range(a_nbins.sum()):
                    cov_tru[i,j] += (hSys[i]-hCV[i])*(hSys[j]-hCV[j])/(hCV[i]*hCV[j])
        
        
        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(range(a_nbins.sum()+1),range(a_nbins.sum()+1))
            crat_tru = ax.pcolormesh(X, Y,cov_tru,cmap='cool')
            cbar = fig.colorbar(crat_tru)
            
            s_lab  = ''
            for chi in range(len(self.a_channelName)):
                s_lab += a_dvar[chi].mylabel+' '
                ax.axvline(a_nbins[:chi].sum(),color='red')
                ax.axhline(a_nbins[:chi].sum(),color='red')   
            ax.set_xlabel(s_lab,fontsize=20)
            ax.set_ylabel(s_lab,fontsize=20)
            
        return cov_tru
        
    def Flat(self,a_dvar,_a_nbins,draw=False):
        a_nbins = np.asarray(_a_nbins)
        if len(self.a_channelName)==0:
            print('Need at least one channel, jabroni.')
            return
        if len(a_dvar) != len(self.a_channelName) or len(a_nbins) != len(self.a_channelName):
            print('variable and bin arrays must match number of channels: ',len(self.a_channelName))
        
        cov_flat = np.zeros((a_nbins.sum(),a_nbins.sum()))
        for sysi in range(len(self.a_detSysName)):   # for each systematic
            if self.a_detSysEnabled[sysi]==0:
                continue
            
            hCV = np.zeros(a_nbins.sum())
            hSys = np.zeros(a_nbins.sum())
            
            for chi in range(len(self.a_channelSCut)):  # for each channel
                # get the list of values for the variable assigned to this channel, in the sample assigned to this channel, after applying the cut assigned to this channel
                var_sys = self.a_OverlapSysDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[a_dvar[chi].myname]
                var_cv = self.a_OverlapCVDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[a_dvar[chi].myname]
                    
                h0,_ = np.histogram(var_cv,bins=1,range=a_dvar[chi].myrange)
                h1,_ = np.histogram(var_sys,bins=1,range=a_dvar[chi].myrange)
                              
                # fill in the corresponding parts of our grand histogram
                hCV[a_nbins[:chi].sum():a_nbins[:chi+1].sum()] += h0[0]
                hSys[a_nbins[:chi].sum():a_nbins[:chi+1].sum()] += h1[0]
        
        for i in range(a_nbins.sum()):
                for j in range(a_nbins.sum()):
                    cov_flat[i,j] += (hSys[i]-hCV[i])*(hSys[j]-hCV[j])/(hCV[i]*hCV[j])
    
        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(range(a_nbins.sum()+1),range(a_nbins.sum()+1))
            crat_tru = ax.pcolormesh(X, Y,cov_flat,cmap='cool')
            cbar = fig.colorbar(crat_tru)
            
            s_lab  = ''
            for chi in range(len(self.a_channelName)):
                s_lab += a_dvar[chi].mylabel+' '
                ax.axvline(a_nbins[:chi].sum(),color='red')
                ax.axhline(a_nbins[:chi].sum(),color='red')   
            ax.set_xlabel(s_lab,fontsize=20)
            ax.set_ylabel(s_lab,fontsize=20)
        return cov_flat
    
    def Polyfit(self,a_dvar,_a_nbins,draw=False):
        a_nbins = np.asarray(_a_nbins)
        if len(self.a_channelName)==0:
            print('Need at least one channel, jabroni.')
            return
        if len(a_dvar) != len(self.a_channelName) or len(a_nbins) != len(self.a_channelName):
            print('variable and bin arrays must match number of channels: ',len(self.a_channelName))
        
        
        cov_poly = np.zeros((a_nbins.sum(),a_nbins.sum()))
        
        for sysi in range(len(self.a_detSysName)):   # for each systematic
            if self.a_detSysEnabled[sysi]==0:
                continue
            
            hCV = np.zeros(a_nbins.sum())
            hSys = np.zeros(a_nbins.sum())
            
            for chi in range(len(self.a_channelSCut)):  # for each channel
                # get the list of values for the variable assigned to this channel, in the sample assigned to this channel, after applying the cut assigned to this channel
                var_sys = self.a_OverlapSysDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[a_dvar[chi].myname]
                var_cv = self.a_OverlapCVDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[a_dvar[chi].myname]
                    
                h0,binedges = np.histogram(var_cv,bins=a_nbins[chi],range=a_dvar[chi].myrange)
                h1,_ = np.histogram(var_sys,bins=a_nbins[chi],range=a_dvar[chi].myrange)
                bincenters = np.diff(binedges)/2 + binedges[:-1]
                
                truRat = np.true_divide(h1,h0,out=np.ones_like(bincenters),where=h0!=0)
                aics = []
                degs = []
                for deg in range(min(np.int(a_nbins[chi]/2),3)):
                    params = deg + 1
                    polyRat = np.polyfit(bincenters, truRat, deg)
                    fRat = np.poly1d(polyRat)
                
                    # now calculate chi2 for fit
                    yerr_rat = np.true_divide(np.sqrt(fRat(bincenters)*h0),h0,out=np.zeros_like(bincenters),where=h0!=0)
                    chi2_fit = np.power(np.true_divide(fRat(bincenters)-truRat,yerr_rat,out=np.zeros_like(bincenters),where=yerr_rat!=0),2).sum()
                    aic = chi2_fit + 2*params + 2*params*(params+1)/float(a_nbins[chi]-params-1)
                    aics.append(aic)
                    degs.append(deg)
                polyterms = degs[np.argmin(aics)]
                print(self.a_detSysName[sysi],self.a_channelName[chi],'Polyfit Degrees:',polyterms,aics[np.argmin(aics)])  

                polyRat = np.polyfit(bincenters, np.true_divide(h1,h0,where=h0!=0), polyterms)
                fRat = np.poly1d(polyRat) 
                h1_fit = fRat(bincenters)*h0
                
                hCV[a_nbins[:chi].sum():a_nbins[:chi+1].sum()] += h0
                hSys[a_nbins[:chi].sum():a_nbins[:chi+1].sum()] += h1_fit
            
        for i in range(a_nbins.sum()):
                for j in range(a_nbins.sum()):
                    cov_poly[i,j] += (hSys[i]-hCV[i])*(hSys[j]-hCV[j])/(hCV[i]*hCV[j])
    
        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(range(a_nbins.sum()+1),range(a_nbins.sum()+1))
            crat_tru = ax.pcolormesh(X, Y,cov_poly,cmap='cool')
            cbar = fig.colorbar(crat_tru)
            
            s_lab  = ''
            for chi in range(len(self.a_channelName)):
                s_lab += a_dvar[chi].mylabel+' '
                ax.axvline(a_nbins[:chi].sum(),color='red')
                ax.axhline(a_nbins[:chi].sum(),color='red')   
            ax.set_xlabel(s_lab,fontsize=20)
            ax.set_ylabel(s_lab,fontsize=20)
    
        return cov_poly
    
    def DetectorSystematicDiagnostics(self,dvar,nbins,chi):
    
        flatsys = 0.0
        for sysi in range(len(self.a_detSysName)):   # for each systematic
            if self.a_detSysEnabled[sysi]==0:
                continue
            
            var_sys = self.a_OverlapSysDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[dvar.myname]
            print(dvar.myname)
            var_cv = self.a_OverlapCVDf[self.a_channelSampleID[chi]][sysi].query(self.a_channelSCut[chi])[dvar.myname]
                    
            hCV,binedges = np.histogram(var_cv,bins=nbins,range=dvar.myrange)
            hSys,_ = np.histogram(var_sys,bins=nbins,range=dvar.myrange)
            bincenters = np.diff(binedges)/2 + binedges[:-1]
        
            truRat = np.true_divide(hSys,hCV,out=np.ones_like(bincenters),where=hCV!=0)
            aics = []
            degs = []
            for deg in range(min(np.int(nbins/2),3)):
                params = deg + 1
                polyRat = np.polyfit(bincenters, truRat, deg)
                fRat = np.poly1d(polyRat)
                
                # now calculate chi2 for fit
                yerr_rat = np.true_divide(np.sqrt(fRat(bincenters)*hCV),hCV,out=np.zeros_like(bincenters),where=hCV!=0)
                chi2_fit = np.power(np.true_divide(fRat(bincenters)-truRat,yerr_rat,out=np.zeros_like(bincenters),where=yerr_rat!=0),2).sum()
                aic = chi2_fit + 2*params + 2*params*(params+1)/float(nbins-params-1)
                aics.append(aic)
                degs.append(deg)
            polyterms = degs[np.argmin(aics)]
            print(self.a_detSysName[sysi],self.a_channelName[chi],'Polyfit Degrees:',polyterms,aics[np.argmin(aics)])  
            polyRat = np.polyfit(bincenters, np.true_divide(hSys,hCV,where=hCV!=0), polyterms)
            fRat = np.poly1d(polyRat) 
            hSys_fit = fRat(bincenters)*hCV
    
            fig,ax = plt.subplots(figsize=(27,11))    
            gs = gridspec.GridSpec(1, 2)
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1])
    
            dvarLinspace = np.linspace(dvar.myrange[0],dvar.myrange[1],40)
            ax0.scatter(bincenters,np.true_divide(hSys,hCV,where=hCV!=0),label='Variation/CV',s=100)
            ax0.plot(dvarLinspace,fRat(dvarLinspace),label='PolyFit',color='green')
            
            ax0.set_title(self.a_detSysName[sysi],fontsize=30)
            ax0.set_xlabel(dvar.mylabel,fontsize=20)
            ax0.set_xlim(dvar.myrange)
            ax0.set_ylim(0,2)
            ax0.legend(fontsize=15)
    
            ax1.hist(var_cv,nbins,range=dvar.myrange,histtype='step',linewidth=2,label='CV Raw')
            ax1.hist(var_sys,nbins,range=dvar.myrange,histtype='step',linewidth=3,linestyle='--',label='Variation Raw')
            ax1.plot(bincenters,hSys_fit,label='PolyFit Ratio x CV',c='green')
            ax1.legend(fontsize=15)
            ax1.set_xlim(dvar.myrange)
