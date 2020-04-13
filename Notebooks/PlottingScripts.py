import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
from numpy import mean
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from numpy import asarray
from root_pandas import read_root
from matplotlib import gridspec
from scipy import stats
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from textwrap import wrap
import copy

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
        self._df = samp_df
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

    def applycut(self,s_cut):
        newhist = copy.deepcopy(self)
        newhist._df = self._df.query(s_cut)
        newhist.setweight(self._wi)
        return newhist

class distVar:

    def __init__(self,s_name,n_range):
        self.myname = s_name
        self.myrange = n_range


def distplot_wratio(myvar,nbins,myrange,predhists,datahist,stxcoord,legpos=0,ymax=-1,normshift=1,fs=(16,11)):

    fig,ax = plt.subplots(figsize=fs)
    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)
    vals_mc_raw = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)

    a_scaledweight = []
    a_labels_evts = []
    a_colors = []
    a_hists = []

    for i in range(len(predhists)):
        h1_raw,binedges = np.histogram(predhists[i].dist(myvar),nbins,range=myrange,weights=predhists[i]._wgt)     # hist of raw event weights
        h1 = h1_raw * predhists[i]._scale
        vals_mc_raw += h1_raw
        vals_mc += h1
        yerrsq_mc += np.power(np.sqrt(h1_raw)*predhists[i]._scale,2)

        a_scaledweight.append(predhists[i]._wgt*predhists[i]._scale)
        a_labels_evts.append(predhists[i]._label+' (%.2f)'%h1.sum())
        a_colors.append(predhists[i]._color)
        a_hists.append(predhists[i].dist(myvar))

    yerr_mc = np.sqrt(yerrsq_mc)
    vals_data_raw,_ = np.histogram(datahist.dist(myvar),nbins,range=myrange,weights=datahist._wgt)
    vals_data = vals_data_raw * datahist._scale
    yerr_data = np.sqrt(vals_data_raw)*datahist._scale

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, .75])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    if ymax < 0: ymax = max(vals_data.max(),vals_mc.max())*1.35

    ax0.set_ylim(0,ymax)
    ax0.set_xlim(myrange)
    ax1.set_ylim(.5,1.5)
    ax1.set_xlim(myrange)
    ax1.set_xlabel(stxcoord,fontsize=20)
    ax0.set_ylabel('Events in 5e19 POT',fontsize=20)
    ax1.set_ylabel('Data/MC',fontsize=20)
    ax0.set_title('MCC9 Data/MC',fontsize=30)

    ax0.hist(a_hists,nbins,range=myrange,weights=a_scaledweight,color=a_colors,stacked=True,linewidth=0,label=a_labels_evts,edgecolor=None)
    ax0.errorbar(bincenters,vals_data,fmt='.',yerr=yerr_data,color='black',capsize=5,label=datahist._label+' (%i)'%vals_data.sum())

    ax0.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    # Ok. we're gonna do error differently!
    er_rat_dots = np.true_divide(yerr_data,vals_mc,where=vals_mc!=0)
    er_rat_line = np.true_divide(yerr_mc,vals_mc,where=vals_mc!=0)

    #chisq = np.true_divide(np.power(vals_data-vals_mc,2),np.power(yerr_mc,2)+np.power(yerr_data,2)).sum()     # poisson approx
    chisq = 0        # CNP approx
    ndof = 0
    for i in range(len(vals_data)):
        if vals_mc_raw[i] > 0:
            ndof += 1
            if vals_data[i] > 0:
                chisq += np.true_divide(np.power(vals_data[i]-(vals_mc[i]*normshift),2),(3.0*vals_data[i]*vals_mc[i]*normshift)/(vals_mc[i]*normshift+2.0*vals_data[i]))
            else:
                chisq += np.true_divide(np.power(vals_data[i]-(vals_mc[i]*normshift),2),(vals_mc[i]*normshift/2.0))

    #chisq = np.true_divide(np.power(vals_data-(vals_mc*normshift),2),(3*vals_data*vals_mc)/(vals_mc+2*vals_data)).sum() #Chi2 CNP

    pval = 1 - stats.chi2.cdf(chisq, ndof)

    errboxes = []
    for i in range(len(bincenters)):
        rect = Rectangle((binedges[i],normshift-er_rat_line[i]),binedges[i+1]-binedges[i],er_rat_line[i]*2)
        errboxes.append(rect)
    pc = PatchCollection(errboxes, facecolor='red', alpha=.3)
    ax1.add_collection(pc)
    ax1.errorbar(bincenters,np.true_divide(vals_data,vals_mc),yerr=er_rat_dots,fmt='o',color='maroon',capsize=0)

    ax1.hist(np.zeros(1),(1,2),facecolor='red',alpha=.3,zorder=0,label='MC Stat Error')
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
    return fig,ax0



def truthplot(myvar,nbins,myrange,predhists,stxcoord,legpos=0,ymax=-1,normshift=1):

    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)
    vals_troof = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)

    a_scaledweight = []
    a_labels_evts = []
    a_colors = []
    a_hists = []

    for i in range(len(predhists)):
        h1_raw,binedges = np.histogram(predhists[i].dist(myvar),nbins,range=myrange,weights=predhists[i]._wgt)     # hist of raw event weights
        h1 = h1_raw * predhists[i]._scale
        a_scaledweight.append(predhists[i]._wgt*predhists[i]._scale)

        yerrsq_mc += np.power(np.sqrt(h1_raw)*predhists[i]._scale,2)

        vals_mc += h1
        a_labels_evts.append(predhists[i]._label+' (%.2f)'%h1.sum())
        a_colors.append(predhists[i]._color)
        a_hists.append(predhists[i].dist(myvar))

        htroof,_ = np.histogram(predhists[i].dist('MC_energyInit'),nbins,range=myrange,weights=predhists[i]._wgt)
        vals_troof += htroof* predhists[i]._scale


    yerr_mc = np.sqrt(yerrsq_mc)

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    if ymax < 0: ymax = vals_mc.max()*1.35

    ax.set_ylim(0,ymax)
    ax.set_xlim(myrange)
    ax.set_ylabel('Events in 5e19 POT',fontsize=20)
    ax.set_title('MCC9 MC/Truth',fontsize=30)

    ax.hist(a_hists,nbins,range=myrange,weights=a_scaledweight,color=a_colors,stacked=True,linewidth=0,label=a_labels_evts,edgecolor=None)
    plt.scatter(bincenters,vals_troof,color='black',label='Truth',zorder=10)

    ax.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    return ax
