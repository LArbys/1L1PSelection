import platform
from ROOT import TFile,vector
import ROOT
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
from numpy import mean
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from numpy import asarray
from root_pandas import read_root
from matplotlib import gridspec
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import math


import pickle
import pandas as pd
import numpy as np
from numpy import mean
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from numpy import asarray
from root_pandas import read_root
from matplotlib import gridspec
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import quad
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from textwrap import wrap
import copy

import seaborn as sns

def distplot(nbins,myrange,_hists,data,data_weight,stxcoord,ymax):

    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)

    for i in xrange(0,len(_hists)):
        v1,binedges = np.histogram(_hists[i],nbins,range=myrange,weights=a_weights[i])
        vals_mc += v1

    vals_data,_ = np.histogram(var_data,nbins,range=myrange,weights=data_weight)

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    plt.ylim(0,ymax)
    plt.xlim(myrange)
    ax.set_xlabel(stxcoord,fontsize=20)
    ax.set_ylabel('Events in 7e20 POT',fontsize=20)
    ax.set_title('MCC9 Data/MC',fontsize=30)

    plt.hist(_hists,nbins,range=myrange,weights=a_weights,color=a_colors,stacked=True,linewidth=2.5,label=a_labels,edgecolor='None')
    ax.errorbar(bincenters,vals_data,fmt='.',yerr=np.sqrt(vals_data/DATA_POT*POT_SCALE_FACTOR),color='black',capsize=5,label='bnb5e19')

    ax.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    print ('Events:',vals_data.sum())
    print ('MC: ',vals_mc.sum())

def distplot_normed(nbins,myrange,_hists,data,data_weight,stxcoord,ymax):

    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)


    for i in xrange(0,len(_hists)):
        v1,binedges = np.histogram(_hists[i],nbins,range=myrange,weights=a_weights[i])
        vals_mc += v1

    vals_data,_ = np.histogram(var_data,nbins,range=myrange,weights=data_weight,density=1)

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    plt.ylim(0,ymax)
    plt.xlim(myrange)
    ax.set_xlabel(stxcoord,fontsize=20)
    ax.set_ylabel('Area Normalized',fontsize=20)
    ax.set_title('MCC9 Data/MC [AREA NORMED]',fontsize=30)

    plt.hist(_hists,nbins,range=myrange,weights=a_weights,color=a_colors,stacked=True,linewidth=2.5,label=a_labels,density=1,edgecolor='None')
    ax.errorbar(bincenters,vals_data,fmt='.',yerr=np.sqrt(vals_data/DATA_POT*POT_SCALE_FACTOR),color='black',capsize=5,label='bnb5e19')

    ax.legend(loc='upper right',fontsize=15,frameon=False)
    print ('Events:',vals_data.sum())

def distplot_wratio(nbins,myrange,_hists,data,data_weight,stxcoord,ymax,writepval=True):

    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)
    yerr_mc = np.zeros(nbins)
    vals_mc_raw = np.zeros(nbins)

    a_labels_evts = []

    print(len(_hists),len(a_weights))

    for i in range(0,len(_hists)):
        print(len(_hists[i]),len(a_weights[i]))
        v1,binedges = np.histogram(_hists[i],nbins,range=myrange,weights=a_weights[i])
        yerr_mc += v1*a_scale[i]
        vals_mc += v1
        v1_raw,binedges_raw = np.histogram(_hists[i],nbins,range=myrange)
        vals_mc_raw+=v1_raw
        a_labels_evts.append(a_labels[i]+' (%.2f)'%v1.sum())

    yerr_mc = np.sqrt(yerr_mc)

#     vals_data_err,_err= np.histogram(data,nbins,range=myrange)
    vals_data,_ = np.histogram(data,nbins,range=myrange,weights=data_weight)
    vals_data_err = vals_data*DATA_SCALE
    bincenters = np.diff(binedges)/2 + binedges[:-1]

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    ax0.set_ylim(0,ymax)
    ax0.set_xlim(myrange)
    ax1.set_ylim(0,2)
    ax1.set_xlim(myrange)
    ax1.set_xlabel(stxcoord,fontsize=20)
    ax0.set_ylabel('Events in 5e19 POT',fontsize=20)
    ax1.set_ylabel('Data/MC',fontsize=20)
    ax0.set_title('MCC9 Data/MC',fontsize=30)

    ax0.hist(_hists,nbins,range=myrange,weights=a_weights,color=a_colors,stacked=True,linewidth=0,label=a_labels_evts,edgecolor=None)
    ax0.errorbar(bincenters,vals_data,fmt='.',yerr=np.sqrt(vals_data/DATA_POT*POT_SCALE_FACTOR),color='black',capsize=5,label='bnb5e19(scaled) (%i)'%vals_data.sum())

    ax0.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    er_rat = np.sqrt(np.power(np.divide(yerr_mc,vals_mc),2) + np.power(np.divide(np.sqrt(vals_data_err),vals_data),2)) * np.true_divide(vals_data,vals_mc)

    chisq = 0.0
    m_cov = np.zeros((nbins,nbins))
    ndof = 0
    normshift =1
    yerrsq_data = np.zeros(nbins)
    for i in range(nbins):
        if vals_mc_raw[i] > 0:
            ndof += 1
            if vals_data[i] > 0:
                yerrsq_data[i] += (3.0*vals_data[i]*vals_mc[i]*normshift)/(vals_mc[i]*normshift+2.0*vals_data[i])
#                 print(yerrsq_data[i])
            else:
                yerrsq_data[i] += vals_mc[i]*normshift/2.0
            m_cov[i][i] += yerrsq_data[i]
        else:
            m_cov[i][i] += 999

    yerr_data = np.sqrt(yerrsq_data)

#     print(yerr_data)
#     print(m_cov)

    invcov = np.linalg.inv(m_cov)
#     print(invcov)
#     invcov= m_cov
    # calc chi2
    for i in range(nbins):
        for j in range(nbins):
            if vals_mc_raw[i] > 0 and vals_mc_raw[j] > 0:
#                 print((vals_mc[i]*normshift - vals_data[i]) * (vals_mc[j]*normshift - vals_data[j]) * invcov[i][j])
                chisq += (vals_mc[i]*normshift - vals_data[i]) * (vals_mc[j]*normshift - vals_data[j]) * invcov[i][j]
    pval = 1 - stats.chi2.cdf(chisq, ndof)
    print ("stats:",chisq,pval)


    ax1.errorbar(bincenters,np.true_divide(vals_data,vals_mc),yerr=er_rat,fmt='o',color='maroon',capsize=0)
    ax1.axhline(1,color='maroon',linestyle=':')
    ax0.annotate(r'$\sum$data/$\sum$MC = %.2f'%(vals_data.sum()/float(vals_mc.sum())),xy=(.01,.92),xycoords='axes fraction',fontsize=20,bbox=dict(boxstyle="square", fc="ghostwhite",alpha=.8))

    ax1.annotate('p-value: %.3f\n'%pval+r'$\chi^2_{CNP}/%i  (dof)$: %.3f'%(ndof,chisq/float(ndof)),xy=(.85,.7), xycoords='axes fraction',fontsize=15,bbox=dict(boxstyle="round4", fc="w",alpha=.9))

    plt.tight_layout()
    print ('Events:',vals_data.sum())
    print ('MC: ',vals_mc.sum())
    return ax0


def Tune1(df):
    return df['xsec_tune1_weight'].values

def CV(df):
    wgt = df['xsec_corr_weight'].values
    wgt[wgt == np.inf] = 1
    wgt = np.nan_to_num(wgt)
    wgt[wgt <= 0] = 1
    return wgt

def Spline(df):
    wgt = df['spline_weight'].values
    wgt[wgt == np.inf] = 1
    wgt = np.nan_to_num(wgt)
    wgt[wgt <= 0] = 1
    return wgt


class sampHist:

    def __init__(self,samp_df,samp_l,samp_c,samp_wind,samp_s):
        self._df = samp_df.copy()
        self._label = samp_l
        self._color = samp_c
        self._scale = samp_s
        self._wi = samp_wind
        if samp_wind == 0:
            self._wgt = np.ones(len(samp_df))
        if samp_wind == 1:
            self._wgt = CV(samp_df)

    def setweight(self,windex):
        if windex == 0:
            self._wgt = np.ones(len(self._df))
        if windex == 1:
            self._wgt = CV(self._df)

    def dist(self,s_var):
        return self._df[s_var].values

    def cosdist(self,s_var):
        return np.cos(self._df[s_var].values)

    def applycut(self,s_cut):
        newhist = copy.deepcopy(self)
        newhist._df = self._df.query(s_cut)
        newhist.setweight(self._wi)
        return newhist

class SimpleHisto:

    def __init__ (self,df_df,f_scale,i_wgt,s_color,s_label):
        temp_df = df_df.copy()
        temp_df['myscale'] = f_scale
        self.mydf = temp_df
        self.mycolor = s_color
        self.mylabel = s_label
        self.iwgt = i_wgt

        self.mycut = 'Enu_1m1p > 0'

    def AddCut(self,s_cut):
        self.mycut  = s_cut

    def ClearCut(self):
        self.mycut = 'Enu_1m1p > 0'

    def GetHist(self,s_varname):
        temp_df = self.mydf.query(self.mycut)
        if self.iwgt == 1:
            myweight = CV(temp_df)
        elif self.iwgt == 0:
            myweight = np.ones(len(temp_df))
        return temp_df[s_varname].values,self.mycolor,self.mylabel,myweight,temp_df['myscale'].values


class StackedHisto:

    def  __init__ (self,a_df_mc,a_df_scale):
        self.mystrata = []
        self.stratalabel = []
        self.stratacolor = []
        self.mylayer = []
        self.layerlabel = []
        self.layercolor = []
        self.layeriwgt = []

        temp_a_df = []
        for i in range(len(a_df_mc)):
            temp_df = a_df_mc[i].copy()
            temp_df['myscale'] = a_df_scale[i]
            temp_a_df.append(temp_df)
        self.mymc = pd.concat(temp_a_df)


    def ClearCut(self):
        self.mycut = 'Enu_1m1p > 0'

    def AddCut(self,s_cut):
        self.mycut  = s_cut

    def AddLayer(self,df_layer,df_scale,i_wgt,s_label,s_color):
        temp_df = df_layer.copy()
        temp_df['myscale'] = df_scale
        self.mylayer.append(temp_df)
        self.layerlabel.append(s_label)
        self.layercolor.append(s_color)
        self.layeriwgt.append(i_wgt)

    def AddStrata(self,s_strata,s_label,s_color):
        self.mystrata.append(s_strata)
        self.stratalabel.append(s_label)
        self.stratacolor.append(s_color)

    def GetHists(self,s_varname):
        a_vals = []
        a_cols = []
        a_labels = []
        a_wgts = []
        a_scale = []

        for i in range(len(self.mystrata)):
            temp_df = self.mymc.query(self.mystrata[i])
            a_vals.append(temp_df[s_varname].values)
            a_cols.append(self.stratacolor[i])
            a_labels.append(self.stratalabel[i])
            a_wgts.append(CV(temp_df))
            a_scale.append(temp_df['myscale'].values)

        for i in range(len(self.mylayer)):
            temp_df = self.mylayer[i]
            a_vals.append(temp_df[s_varname].values)
            a_cols.append(self.layercolor[i])
            a_labels.append(self.layerlabel[i])
            if self.layeriwgt[i] == 1:
                a_wgts.append(CV(temp_df))
            elif self.layeriwgt[i] == 0:
                a_wgts.append(np.ones(len(temp_df)))
            a_scale.append(temp_df['myscale'].values)

        return  a_vals,a_cols,a_labels,a_wgts,a_scale

    def ClearLayers(self):
        self.mylayer = []
        self.layerlabel = []
        self.layercolor = []
        self.layerwgt = []

    def ClearStrata(self):
        self.mystrata = []
        self.stratalabel = []
        self.stratacolor = []


class distVar:

    def __init__(self,s_name,n_range,s_label='',s_cov=''):
        self.myname = s_name
        self.myrange = n_range
        self.myscov = s_cov
        if s_label == '':
            self.mylabel = s_name
        else:
            self.mylabel = s_label


def distplot_wratio_davio(myvar,nbins,myrange,stackedhists,datahist,stxcoord,legpos=0,ymax=-1,normshift=1,fs=(16,11),s_cov=''):

    fig,ax = plt.subplots(figsize=fs)
    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)
    vals_mc_raw = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)
    yerr_mc_sys = np.zeros(nbins)
    ndof = 0

    a_labels_evts = []

    gh_vals,gh_cols,gh_labels,gh_wgts,gh_scale = stackedhists.GetHists(myvar)
    data_vals,_,data_label,data_wgt,data_scale = datahist.GetHist(myvar)

    # for poisson error bars:
    a_hstack = []
    a_hscale = []

    for i in range(len(gh_vals)):
        h1_raw,_ = np.histogram(gh_vals[i],nbins,range=myrange,weights=gh_wgts[i])     # hist of raw event weights
        h1,_ = np.histogram(gh_vals[i],nbins,range=myrange,weights=gh_wgts[i]*gh_scale[i])

        vals_mc_raw += h1_raw
        vals_mc += h1

        scalesort = gh_scale[i].argsort()
        sorted_vals = gh_vals[i][scalesort[::-1]]
        sorted_scale = gh_scale[i][scalesort[::-1]]
        sorted_wgt = gh_wgts[i][scalesort[::-1]]
        for sc in np.unique(sorted_scale):
            subvals = sorted_vals[sorted_scale==sc]
            subwgts = sorted_wgt[sorted_scale==sc]
            subh1,_ = np.histogram(subvals,nbins,range=myrange,weights=subwgts)
            yerrsq_mc += np.power(np.sqrt(subh1)*sc,2)
            a_hstack.append(subh1)
            a_hscale.append(sc)

        a_labels_evts.append(gh_labels[i]+' (%.2f)'%h1.sum())

    vals_data_raw,binedges = np.histogram(data_vals,nbins,range=myrange,weights=data_wgt)
    vals_data,_ = np.histogram(data_vals,nbins,range=myrange,weights=np.multiply(data_wgt,data_scale))
    bincenters = np.diff(binedges)/2 + binedges[:-1]

    if s_cov == '':
        m_cov = np.zeros((nbins,nbins))
        yerr_mc_total = np.sqrt(yerrsq_mc)
    else:
        m_cov = np.genfromtxt(s_cov,delimiter=',')
        for i in range(nbins):
            for j in range(nbins):
                if(vals_mc[i] > 0 and vals_mc[j]):
                    m_cov[i][j] *= vals_mc[j]*vals_mc[j]
                else:
                    m_cov[i][j] = 0
        yerr_mc_sys = np.sqrt(np.diag(m_cov))
        print(yerr_mc_sys)
        yerr_mc_total = np.sqrt(np.diag(m_cov) + yerrsq_mc)
        print(yerr_mc_total)

    yerrsq_data = np.zeros(nbins)
    for i in range(nbins):
        for j in range(nbins):
            if (m_cov[i][j]!=m_cov[i][j]):
                m_cov[i][j] = 999

    for i in range(nbins):
        if vals_mc_raw[i] > 0:
            ndof += 1
            if vals_data[i] > 0:
                yerrsq_data[i] += (3.0*vals_data_raw[i]*vals_mc[i]*normshift)/(vals_mc[i]*normshift+2.0*vals_data_raw[i])
            else:
                yerrsq_data[i] += vals_mc[i]*normshift/2.0
            print("sqdata",yerrsq_data[i],i)
            m_cov[i][i] += yerrsq_data[i]*np.power(data_scale[0],2)
        else:
            m_cov[i][i] += 999
    print("data_scale", data_scale)
#     yerr_data = np.sqrt(yerrsq_data)*data_scale[0]
    yerr_data=np.sqrt(vals_data/1.0)

#     print('test',data_scale)

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, .75])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    if ymax < 0: ymax = max(vals_data.max(),vals_mc.max())*1.35

    ax0.set_ylim(0,ymax)
    ax0.set_xlim(myrange)
    ax1.set_ylim(.5,1.5)
    ax1.set_xlim(myrange)
    ax1.set_xlabel(stxcoord,fontsize=20)
    ax0.set_ylabel('Events in 7e20 POT',fontsize=20)
    ax1.set_ylabel('Data/MC',fontsize=20)
    ax0.set_title('MCC9 Data/MC',fontsize=30)

    ax0.hist(gh_vals,nbins,range=myrange,weights=[gh_wgts[i]*gh_scale[i] for i in range(len(gh_wgts))],color=gh_cols,stacked=True,linewidth=0,label=a_labels_evts,edgecolor=None)
    ax0.errorbar(bincenters,vals_data,fmt='.',yerr=yerr_data,color='black',capsize=5,label=data_label+' (%i)'%vals_data.sum())

    if s_cov!='':
        errboxes = []
        for i in range(len(bincenters)):
            rect = Rectangle((binedges[i],(vals_mc[i]-yerr_mc_sys[i])),binedges[i+1]-binedges[i],yerr_mc_sys[i]*2)
            errboxes.append(rect)
        pc = PatchCollection(errboxes,facecolor=None,alpha=.1,hatch='/',zorder=12)
        ax0.add_collection(pc)
        ax0.hist(np.zeros(1),(1,2),facecolor=None,alpha=.1,hatch='//',zorder=0,label='Flux + Cross Section Systematics')
        print("errboxes",len(errboxes))

    ax0.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    er_rat_dots = np.true_divide(yerr_data,vals_mc,where=vals_mc!=0)
    er_rat_line = np.true_divide(yerr_mc_total,vals_mc,where=vals_mc!=0)
    er_rat_line_sys = np.true_divide(yerr_mc_sys,vals_mc,where=vals_mc!=0)

    chisq = 0.0
#     print("mcov",m_cov)
    invcov = np.linalg.inv(m_cov)
    # calc chi2
    for i in range(nbins):
        for j in range(nbins):
            if vals_mc_raw[i] > 0 and vals_mc_raw[j] > 0:
                chisq += (vals_mc[i]*normshift - vals_data[i]) * (vals_mc[j]*normshift - vals_data[j]) * invcov[i][j]
#                 print(invcov[i][j])
    pval = 1 - stats.chi2.cdf(chisq, ndof)

    errboxes_tot = []
    errboxes_sys = []
    for i in range(len(bincenters)):
        rect_tot = Rectangle((binedges[i],normshift-er_rat_line[i]),binedges[i+1]-binedges[i],er_rat_line[i]*2)
        errboxes_tot.append(rect_tot)
        if s_cov != '':
            rect_sys = Rectangle((binedges[i],normshift-er_rat_line_sys[i]),binedges[i+1]-binedges[i],er_rat_line_sys[i]*2)
            errboxes_sys.append(rect_sys)
    pc_tot = PatchCollection(errboxes_tot, facecolor='red', alpha=.3)
    pc_sys = PatchCollection(errboxes_sys, facecolor=None, alpha=.1,hatch='/',zorder=12)
    ax1.add_collection(pc_tot)
    if s_cov != '':
        ax1.add_collection(pc_sys)
        ax1.hist(np.zeros(1),(1,2),facecolor=None,alpha=.1,hatch='//',zorder=0,label='Sys Error')

    ax1.hist(np.zeros(1),(1,2),facecolor='red',alpha=.3,zorder=0,label='Total Error')
    ax1.errorbar(bincenters,np.true_divide(vals_data,vals_mc),yerr=er_rat_dots,fmt='o',color='maroon',capsize=0)


    legloc = ['upper left','lower left','lower right']
    ax1.legend(loc=legloc[legpos],fontsize=15,frameon=False)

    ax1.axhline(1,color='black',linestyle=':')
    if normshift != 1:
        ax1.axhline(normshift,color='maroon')
    ax0.annotate(r'$\sum$data/$\sum$MC = %.2f'%(vals_data.sum()/float(vals_mc.sum())),xy=(.01,.92),xycoords='axes fraction',fontsize=20,bbox=dict(boxstyle="square", fc="ghostwhite",alpha=.8))

    ax1.annotate('p-value: %.3f\n'%pval+r'$\chi^2_{CNP}/%i  (dof)$: %.3f'%(ndof,chisq/float(ndof)),xy=(.85,.7), xycoords='axes fraction',fontsize=15,bbox=dict(boxstyle="round4", fc="w",alpha=.9))

    plt.tight_layout()
    print('Events:',vals_data.sum())
    print('Min Bin ct (data):',vals_data.min())
    return fig,ax0,pval


def truthplot(myvar,nbins,myrange,predhists,stxcoord,legpos=0,ymax=-1,normshift=1):

    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)
    vals_troof = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)

    a_labels_evts = []

    for i in range(len(predhists)):
        h1_raw,binedges = np.histogram(predhists[i].dist(myvar),nbins,range=myrange,weights=predhists[i]._wgt)     # hist of raw event weights
        h1 = h1_raw * predhists[i]._scale
        a_scaledweight.append(predhists[i]._wgt*predhists[i]._scale)

        yerrsq_mc += np.power(np.sqrt(h1_raw)*predhists[i]._scale,2)

        vals_mc += h1
        a_labels_evts.append(predhists[i]._label+' (%.2f)'%h1.sum())

        htroof,_ = np.histogram(predhists[i].dist('MC_energyInit'),nbins,range=myrange,weights=predhists[i]._wgt)
        vals_troof += htroof* predhists[i]._scale


    yerr_mc = np.sqrt(yerrsq_mc)

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    if ymax < 0: ymax = vals_mc.max()*1.35

    ax.set_ylim(0,ymax)
    ax.set_xlim(myrange)
    ax.set_ylabel('Events in 7e20 POT',fontsize=20)
    ax.set_title('MCC9 MC/Truth',fontsize=30)

    ax.hist(gh_vals,nbins,range=myrange,weights=a_scaledweight,color=gh_cols,stacked=True,linewidth=0,label=a_labels_evts,edgecolor=None)
    plt.scatter(bincenters,vals_troof,color='black',label='Truth',zorder=10)

    ax.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    return ax

def mcplot(_vartest,_myrange,_mybins,_predhists,_obshist):
    _pvals,_,_,_pwgts,_ = _predhists.GetHists(_vartest)
    _ovals,_,_,_owgts,_ = _obshist.GetHist(_vartest)

    phist = np.zeros(_mybins)
    for i in range(len(_pvals)):
        _ph,_ = np.histogram(_pvals[i],bins=_mybins,range=_myrange,weights=_pwgts[i])
        phist += _ph

    ohist,binedges = np.histogram(_ovals,bins=_mybins,range=_myrange,weights=_owgts)
    bincenters = np.diff(binedges)/2 + binedges[:-1]

    fig,ax = plt.subplots(figsize=(16,11))
    plt.scatter(bincenters,np.true_divide(phist,ohist,where=ohist!=0))
    ax.set_ylabel('Raw MC / Data',fontsize=25)
    ax.set_xlabel(_vartest)
    print('MC raw',phist.sum())
    print('Data raw',ohist.sum())


def gaus(x,a,x0,sigma):
    g = np.zeros(len(x))
    for i in range(len(x)):
        g[i] = a*exp(-(x[i]-x0)**2/(2*sigma**2))
    return g

def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    return a, b
