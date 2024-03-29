"""
Library: lib_plot.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
# Base
import os
import re
import time
import math
import numpy as np
import pandas as pd
import hydroeval as he
import datetime as dtt

# Analysis
from scipy import special as sp
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter as sfilter

# Graphics
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import rc_context

# -----------------------------------------------------------------------------

# Plotter for all columns of data

def plot_columns(df):
    fig, ax = plt.subplots(figsize=(14, 6))

    for column in df.columns:
        if df[column].dtype == 'object':
            continue  # Skip non-numeric columns

        ax.plot(df.index, df[column], label=column)

    ax.legend()
    plt.show()

#############################################################################
# Plotting options (rcParams)
#############################################################################

rc_dict = {
    'axes.titlesize' : 16,
    'axes.labelsize' : 16,
    # lines.linewidth : 3,
    # lines.markersize : 10,
    'xtick.labelsize' : 16,
    'ytick.labelsize' : 16,
    'legend.fontsize' : 'x-large',
}

#############################################################################

def bias(obs, sim):
    """distance between obs' and sim's mean values"""
    if len(obs)==len(sim):
        return np.mean(obs-sim)
    else: raise ValueError(
        f'obs and sim must have same first dimension, but have shapes {np.shape(obs)} and {np.shape(sim)}')

    
def timeseries(dates, data):
    """Returns a matrix (list of type(dates,data)) in the format [dates,data]"""
    
    if len(dates)==len(data):
        return [[dates[i],data[i]] for i in range(len(dates))]
    else: raise ValueError(
        f'dates and data must have same first dimension, but have shapes {np.shape(dates)} and {np.shape(data)}')


#############################################################################
# Plotting options (rcParams)
#############################################################################

rc_dict = {
    'axes.titlesize': 16,
    'axes.labelsize': 16,
    # lines.linewidth : 3,
    # lines.markersize : 10,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 'x-large',
}


#############################################################################
# Triple plot
#############################################################################

# @mplt.rc_context(rc_dict)
def plot_triple(fig, ax, times1: list, data1: list, data1_label: str,
                input1: list, input1_label: str,
                times2: list, data2: list, data2_label: str,
                input2: list, input2_label: str,
                times3: list, data3: list, data3_label: list,
                user_input: list,
                ):
    """3 subplots with 2 sim VS obs timeseries and eventual inputs

    Inputs
    ------
    - data1: list(obs, sim)
    - data1_label: label of the quantity
    - input1: list()
    - data2: list(obs, sim)
    - data2_label: label of the quantity
    - input2: list()
    - ...

    Note: figure, saving options to be defined outside
    """

    mplt.rcParams.update(rc_dict)

    # ----------------------------------------------------------------------
    # Plot 1 sim vs obs timeseries

    obs = data1[0]; obs_label = data1_label + '_obs'
    sim = data1[1]; sim_label = data1_label + '_sim'
    labely = data1_label
    times = [x.to_pydatetime() for x in times1]
    marker = 'o'; linestyle = '-'
    units = r' $[dB]$'

    # RMSE, R, bias, KGE calculation
    RMSE = np.mean((sim - obs) ** 2) ** 0.5;
    print('RMSE =', RMSE)
    R = np.corrcoef(sim, obs)[0][1];
    print('R=', R)
    BIAS = bias(sim, obs);
    print('bias =', BIAS)
    KGE = he.evaluator(he.kge, sim, obs)[0, :][0];
    print('KGE=', KGE)

    title = f'{sim_label} VS {obs_label} - RMSE={RMSE:.2f}, R={R:.2f}, bias={BIAS:.2f}, KGE={KGE:.2f}'

    # ax[0].set_xlim(xmin=times[0], xmax=times[-1]) # xmax=times[len(times)-1]
    ax[0].plot(times, sim, c='tab:red', label=sim_label,
               linestyle=linestyle, lw=1, marker=marker,
               ms=1)  # alpha=.4, zorder=-1)
    ax[0].plot(times, obs, c='tab:blue', label=obs_label,
               linestyle=linestyle, lw=1, marker=marker, ms=1,
               alpha=.4, zorder=-1)
    ax[0].legend(loc='upper left')
    ax[0].set_title(title)
    ax[0].set_ylabel(labely + units, size='xx-large')

    units = r' $[-]$'
    ax0 = ax[0].twinx()
    ax0.plot(times, input1, label=input1_label, color='tab:green', lw=1)
    ax0.legend(loc='upper right')
    ax0.set_ylabel(input1_label + units, size='xx-large')

    # -----------------------------------------------------------------------
    # Plot 2 sim vs obs timeseries

    obs = None
    if len(data2) > 1:
        obs = data2[0];
        obs_label = data2_label + '_obs'
        sim = data2[1];
        sim_label = data2_label + '_sim'
    else:
        sim = data2;
        sim_label = data2_label + '_sim'
    labely = data2_label
    times = [x.to_pydatetime() for x in times2]
    units = r' $[m^3/m^3]$'

    # RMSE, R, bias calculation
    simmatrix = np.array([[sim[i], obs[i]] for i in range(len(sim))
                          if not np.isnan(obs[i])])
    RMSE = np.nanmean((sim - obs) ** 2) ** 0.5;
    print('RMSE =', RMSE)
    R = np.corrcoef(simmatrix, rowvar=False)[0][1];
    print('R (sim vs obs) =', R)
    BIAS = bias(np.array([e[0] for e in simmatrix]),
                np.array([e[1] for e in simmatrix]))
    KGE_SM = he.evaluator(he.kge, sim, obs)[0, :][0];
    print('KGE=', KGE_SM)

    title = f'{sim_label} VS {obs_label} - RMSE={RMSE:.2f}, R={R:.2f}, bias={BIAS:.2f}, KGE={KGE_SM:.2f}'

    ax[1].set_xlim(xmin=times[0], xmax=times[-1])
    # if obs is not None: ax[1].plot(times, obs, c='tab:blue', label=obs_label)
    ax[1].plot(times, sim, c='tab:red', label=sim_label)
    ax[1].legend(loc='upper left')
    ax[1].set_title(title)
    ax[1].set_ylabel(data2_label + units, size='xx-large')

    if len(input2) is not 0:
        ax1 = ax[1].twinx()
        ax1.plot(times, input2, c='tab:blue', label=input2_label)
        ax1.set_ylabel(input2_label + units, size='xx-large')
        ax1.legend(loc='upper right')

    # -----------------------------------------------------------------------
    # Plot of inputs P, EPOT

    label1, label2 = data3_label
    times = [x.to_pydatetime() for x in times3]
    units = r' $[mm]$'

    ax[2].bar(times, data3[0], color='tab:gray', label=label1, )
    ax[2].legend(loc='upper left')
    ax[2].set_ylabel(label1 + ', ' + label2 + units, size='xx-large')

    ax2 = ax[2].twinx()
    ax2.plot(times, data3[1], label=label2, color='tab:green')
    ax2.legend(loc='upper right')
    ax2.set_ylabel(label2 + r' $[mm/day]$', size='xx-large')

    # -----------------------------------------------------------------------

    ax[0].set_ylim(-18, -7)
    ax0.set_ylim(0, 1)
    # ax[1].set_ylim(0.2, 1)
    ax[2].set_ylim(0, 100)
    ax2.set_ylim(0, 10)


#############################################################################
# Scatter plot
#############################################################################

def plot_sim_vs_obs(sim: list, obs: list, quantity: str, um: str):
    import matplotlib.gridspec as gridspec

    def linear(x, a, b):
        return a + b * x

    title = f'{quantity} obs VS simul - '  # y VS x
    xlabel = f'{quantity}_sim {um}'
    ylabel = f'{quantity}_obs {um}'

    data = pd.DataFrame({'sim': sim, 'obs': obs})
    data.dropna(inplace=True)
    x = data.sim.values
    y = data.obs.values

    fig = plt.figure(figsize=(6, 6), dpi=200)
    gs = gridspec.GridSpec(nrows=1, ncols=1, width_ratios=[1],
                           height_ratios=[1])
    ax = plt.subplot(gs[0])
    ax.plot(x, y, marker='o', linestyle='', color='tab:blue')
    min_common = np.min([x, y]) - 0.1 * abs(np.mean([x, y]))
    max_common = np.max([x, y]) + 0.1 * abs(np.mean([x, y]))
    ax.set_xlim(min_common, max_common)
    ax.set_ylim(min_common, max_common)
    lin_grid = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1],
                           100);  # ax.plot(lin_grid, lin_grid, color='k') # y = x

    # Fit
    popt, pcov = curve_fit(linear, x, y)
    ax.plot(lin_grid, linear(np.array(lin_grid), *popt), color='k')

    RMSE = np.nanmean((sim - obs) ** 2) ** 0.5;
    print('RMSE =', RMSE)
    R = np.corrcoef(x, y)[0][1];
    print('R=', R, 'R^2=', R ** 2)
    BIAS = bias(x, y);
    print('bias=', BIAS)

    ax.set_xlabel(xlabel);
    ax.set_ylabel(ylabel)
    xtext = 0.2 * (max_common - min_common) + min_common
    ytext = 0.9 * (max_common - min_common) + min_common
    ax.text(xtext, ytext,
            f'y={popt[0]:.2f}+{popt[1]:.2f}x',
            ha="center", va="center", size=15,
            bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=2,
                      alpha=.5))

    ax.set_title(
        title + f'RMSE={RMSE:.2f}, R={R:.2f},' + r' $R^2$' + f'={R ** 2:.2f}, bias={BIAS:.2f}')
    ax.set_aspect('equal', adjustable='box', share=True)
    
    

#############################################################################
# Double plot
#############################################################################

# @mplt.rc_context(rc_dict)
def plot_double(fig, ax,
                times1: list, data1: list, data1_label: str,
                input1: list, input1_label: str,
                data1_mask: np.ndarray,
                times3: list, data3: list, data3_label: list,
                user_input: list,
                ):
    """2 subplots with 2 sim VS obs timeseries and eventual inputs

    Inputs
    ------
    - data1: list(obs, sim)
    - data1_label: label of the quantity
    - input1: list()
    - data2: list(obs, sim)
    - data2_label: label of the quantity
    - input2: list()
    - ...

    Note: figure, saving options to be defined outside
    """

    mplt.rcParams.update(rc_dict)

    # ----------------------------------------------------------------------
    # Plot 1 sim vs obs timeseries

    obs = data1[0]; obs_label = data1_label + '_obs'
    sim = data1[1]; sim_label = data1_label + '_sim'
    labely = data1_label
    times = [x.to_pydatetime() for x in times1]
    marker = 'o'; linestyle = ''
    units = r''

    # RMSE, R, bias, KGE calculation
    sim_gap = sim[data1_mask]
    obs_gap = obs[data1_mask]
    RMSE = np.nanmean((sim_gap - obs_gap) ** 2) ** 0.5;
    print('RMSE =', RMSE)
    R = np.corrcoef(sim_gap, obs_gap)[0][1];
    print('R=', R)
    BIAS = bias(sim_gap, obs_gap);
    print('bias =', BIAS)
    KGE = he.evaluator(he.kge, sim_gap, obs_gap)[0, :][0];
    print('KGE=', KGE)

    title = f'{sim_label} VS {obs_label} - RMSE={RMSE:.3f}, R={R:.2f}, bias={BIAS:.3f}, KGE={KGE:.3f}'

    # ax[0].set_xlim(xmin=times[0], xmax=times[-1]) # xmax=times[len(times)-1]
    ax[0].plot(times, sim, c='tab:red', label=sim_label,
               linestyle=linestyle, lw=1, marker=marker,
               ms=1)  # alpha=.4, zorder=-1)
    ax[0].plot(times, obs, c='tab:blue', label=obs_label,
               linestyle=linestyle, lw=1, marker=marker, ms=1,
               alpha=.4, zorder=-1)
    ax[0].legend(loc='upper left')
    ax[0].set_title(title)
    ax[0].set_ylabel(labely + units, size='xx-large')
    
    units = r' $[-]$'
    ax0 = ax[0].twinx()
    ax0.plot(times, input1, label=input1_label, color='tab:green', lw=1, zorder=-20)
    ax0.legend(loc='upper right')
    ax0.set_ylabel(input1_label + units, size='xx-large')
    
    ax[0].set_zorder(ax0.get_zorder()+1)
    ax[0].set_frame_on(False)
    

    # -----------------------------------------------------------------------
    # Plot of inputs P, EPOT

    label1, label2 = data3_label
    times = [x.to_pydatetime() for x in times3]
    units = r' $[mm]$'

    ax[1].bar(times, data3[0], color='tab:gray', label=label1, )
    ax[1].legend(loc='upper left')
    ax[1].set_ylabel(label1 + ', ' + label2 + units, size='xx-large')

    ax2 = ax[1].twinx()
    ax2.plot(times, data3[1], label=label2, color='tab:green')
    ax2.legend(loc='upper right')
    ax2.set_ylabel(label2 + r' $[mm/day]$', size='xx-large')

    # -----------------------------------------------------------------------

    
    

# #############################################################################
# # Triple plot
# #############################################################################
#
# # @mplt.rc_context(rc_dict)
# def plot_triple(fig, ax, times1:list, data1:list, data1_label:str,
#                 input1:list, input1_label:str,
#                 times2:list, data2:list, data2_label:str,
#                 input2:list, input2_label:str,
#                 times3:list, data3:list, data3_label:list,
#                 user_input:list,
#                ):
#     """3 subplots with 2 sim VS obs timeseries and eventual inputs
#
#     Inputs
#     ------
#     - data1: list(obs, sim)
#     - data1_label: label of the quantity
#     - input1: list()
#     - data2: list(obs, sim)
#     - data2_label: label of the quantity
#     - input2: list()
#     - ...
#
#     Note: figure, saving options to be defined outside
#     """
#
#     irri = user_input
#     #----------------------------------------------------------------------
#     # Plot 1 sim vs obs timeseries
#
#     obs = data1[0]; obs_label=data1_label+'_obs'
#     sim = data1[1]; sim_label=data1_label+'_sim'
#     labely = data1_label
#     times = times1
#     marker='o'; linestyle='-'
#     units=r' $[dB]$'
#
#     # RMSE, R, bias, KGE calculation
#     RMSE=np.mean((sim-obs)**2)**0.5; print('RMSE =', RMSE)
#     R=np.corrcoef(sim,obs)[0][1]; print('R=', R)
#     BIAS=bias(sim,obs); print('bias =', BIAS)
#     KGE=he.evaluator(he.kge, sim, obs)[0,:][0]; print('KGE=', KGE)
#
#     title=f'{sim_label} VS {obs_label} - RMSE={RMSE:.2f}, R={R:.2f}, bias={BIAS:.2f}, KGE={KGE:.2f}'
#
#     ax[0].set_xlim(xmin=times[0], xmax=times[len(times)-1])
#     ax[0].plot(times, sim, c='tab:red', label=sim_label,
#                linestyle=linestyle, marker=marker, )#alpha=.4, zorder=-1)
#     ax[0].plot(times, obs, c='tab:blue', label=obs_label,
#                linestyle=linestyle, marker=marker, alpha=.4, zorder=-1)
#     ax[0].legend(loc='upper left')
#     ax[0].set_title(title)
#     ax[0].set_ylabel(labely+units, fontsize=16)
#
#     units = r' $[-]$'
#     ax1 = ax[0].twinx()
#     ax1.plot(times, input1, label=input1_label, color='tab:green')
#     ax1.legend(loc='upper right')
#     ax1.set_ylabel(input1_label+units, fontsize=16)
#
#
#     #-----------------------------------------------------------------------
#     # Plot 2 sim vs obs timeseries
#
#     obs = data2[0]; obs_label=data2_label+'_obs'
#     sim = data2[1]; sim_label=data2_label+'_sim'
#     labely = data2_label
#     times = times2
#     units = r' $[m^3/m^3]$'
#
#     # RMSE, R, bias calculation
#     simmatrix = np.array( [ [sim[i], obs[i]] for i in range(len(sim))
#                            if not np.isnan(obs[i]) ] )
#     RMSE=np.nanmean((sim-obs)**2)**0.5; print('RMSE =', RMSE)
#     R=np.corrcoef(simmatrix,rowvar=False)[0][1]; print('R (sim vs obs) =', R)
#     BIAS=bias(np.array([e[0] for e in simmatrix]),
#               np.array([e[1] for e in simmatrix]))
#     KGE_SM=he.evaluator(he.kge, sim, obs)[0,:][0]; print('KGE=', KGE_SM)
#
#     if irri==True:
#         IRRmatrix = np.array( [ [IRR[i], IRR_obs[i]] for i in range(len(IRR))
#                                if not np.isnan(IRR_obs[i]) ] )
#         R_IRR=np.corrcoef(IRRmatrix,rowvar=False)[0][1]
#         print('R_IRR (IRR vs IRR_obs)=', R_IRR)
#         B_IRR=bias(np.array([e[0] for e in IRRmatrix]),
#                    np.array([e[1] for e in IRRmatrix]))
#         irri_title = f'sumIRR_obs={np.sum(IRR_obs):.2f}, '+\
#                      f'sumIRR_sim={np.sum(IRR):.2f}, '+\
#                      f'R_IRR={R_IRR:.2f}, '+\
#                      f'bias_IRR={B_IRR:.2f}, '
#     else: irri_title=''
#
#     title=f'{sim_label} VS {obs_label} - RMSE={RMSE:.2f}, R={R:.2f}, bias={BIAS:.2f}, KGE={KGE_SM:.2f}'+' '+f'{irri_title}'
#
#     ax[1].set_xlim(xmin=times[0], xmax=times[-1])
#     ax[1].plot(times, sim, c='tab:red', label=sim_label)
#     ax[1].plot(times, savgol_filter(obs, window_length=36, polyorder=1, mode='mirror'),
#                c='tab:blue', label=obs_label,
#                linestyle='-', alpha=.7, zorder=-1)
#     ax[1].legend(loc='upper left')
#     ax[1].set_title(title)
#     ax[1].set_ylabel(data2_label+units, fontsize=16)
#
#     #-----------------------------------------------------------------------
#     # Plot of inputs P, IRR, veg
#
#     label1, label2, label3 = data3_label
#     times = times3
#     units = r' $[mm]$'
#
#     ax[2].bar(times, data3[0], color='tab:gray', label=label1)
#     ax[2].bar(times, data3[1], color='tab:blue', label=label2, zorder=2)
#     ax[2].legend(loc='upper left')
#     ax[2].set_ylabel(label1+', '+label2+units, fontsize=16)
#
#     ax2 = ax[2].twinx()
#     ax2.plot(times, data3[2], label=label3, color='tab:green')
#     ax2.legend(loc='upper right')
#     ax2.set_ylabel(label3+r' $[mm/h]$', fontsize=16)
#
#
#
#     ax[0].set_ylim(-17, -3)
#     ax[1].set_ylim(0.06, .4)
#     ax1.set_ylim(0,1)
#     ax[2].set_ylim(0, 25)
#     ax2.set_ylim(0,.3)
#
#
# #############################################################################
# # Scatter plot
# #############################################################################
#
# def plot_sim_vs_obs(sim:list, obs:list, quantity:str, um:str):
#
#     import matplotlib.gridspec as gridspec
#
#     def linear(x,a,b):
#         return a+b*x
#
#     title = f'{quantity} obs VS simul - ' # y VS x
#     xlabel = f'{quantity}_sim {um}'
#     ylabel = f'{quantity}_obs {um}'
#
#     data = pd.DataFrame({'sim': sim,'obs': obs})
#     data.dropna(inplace=True)
#     x = data.sim.values
#     y = data.obs.values
#
#     fig = plt.figure(figsize=(6, 6), dpi=200)
#     gs = gridspec.GridSpec(nrows=1, ncols=1, width_ratios=[1], height_ratios=[1])
#     ax = plt.subplot(gs[0])
#     ax.plot(x, y, marker='o', linestyle='', color='tab:blue')
#     min_common = np.min([x,y])-0.1*abs(np.mean([x,y]))
#     max_common = np.max([x,y])+0.1*abs(np.mean([x,y]))
#     ax.set_xlim(min_common, max_common)
#     ax.set_ylim(min_common, max_common)
#     lin_grid = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100); # ax.plot(lin_grid, lin_grid, color='k') # y = x
#
#     # Fit
#     popt, pcov = curve_fit(linear, x, y)
#     ax.plot(lin_grid, linear(np.array(lin_grid),*popt), color='k')
#
#     RMSE=np.nanmean((sim-obs)**2)**0.5; print('RMSE =', RMSE)
#     R=np.corrcoef(x,y)[0][1]; print('R=', R, 'R^2=', R**2)
#     BIAS=bias(x,y); print('bias=', BIAS)
#
#     ax.set_xlabel(xlabel); ax.set_ylabel(ylabel, fontsize=16)
#     xtext=0.2*(max_common-min_common)+min_common
#     ytext=0.9*(max_common-min_common)+min_common
#     ax.text(xtext, ytext,
#             f'y={popt[0]:.2f}+{popt[1]:.2f}x',
#             ha="center", va="center", size=15,
#             bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=2, alpha=.5))
#
#     ax.set_title(title+f'RMSE={RMSE:.2f}, R={R:.2f},'+r' $R^2$'+f'={R**2:.2f}, bias={BIAS:.2f}')
#     ax.set_aspect('equal', adjustable='box', share=True)
    

#############################################################################
# Data analysis, fit
# Fitting functions
#############################################################################


def linear(x,a,b):
    return a+b*x
    

def gauss(x, A, mean, dev):
    """Not-normalized, shifted gaussian distribution."""
    
    import math
    
    pdf = (1/(dev*np.sqrt(2*math.pi)))*np.exp(-(x-mean)**2/(2*dev**2))
    return A*pdf


def skew_gauss(x, A, mean, dev, alpha,):
    """Skew, not-normalized and shifted gaussian distribution.

    References:
    - https://www.wolframalpha.com/input?i=skew+gaussian+distribution
    - https://stackoverflow.com/questions/15400850/scipy-optimize-curve-fit-unable-to-fit-shifted-skewed-gaussian-curve
    - https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html
    """
    
    import math
    import scipy.special as sp
    
    pdf = (1/(dev*np.sqrt(2*np.pi)))*np.exp(-pow((x-mean),2)/(2*pow(dev,2)))
    cdf = sp.erfc((-alpha*(x-mean))/(dev*np.sqrt(2)))
    return A*pdf*cdf


#############################################################################
# Snippets, old code
#############################################################################
    
# NS=1-np.nansum((sim-obs)**2)/np.nansum((obs-np.nanmean(obs))**2); print('NS =', NS)
#     NS_radQ=1-np.nansum((np.sqrt(sim+0.00001)-np.sqrt(obs+0.00001))**2)/np.nansum((np.sqrt(obs+0.00001)-np.nanmean(np.sqrt(obs+0.00001)))**2)
#     NS_lnQ=1-np.nansum((np.log(sim+0.0001)-np.log(obs+0.0001))**2)/np.nansum((np.log(obs+0.0001)-np.nanmean(np.log(obs+0.0001)))**2)
#     NS_lnQ=NS_lnQ.real; # print(NS_lnQ) 
#     NS_radQ=NS_radQ.real; # print(NS_radQ)

#----------------------------------------------------------------------------
# Etichette su hist_gauss_fit

# xtext = 0.5*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] # 0.8*(max(x)-min(x))+min(x);
# ytext = 0.5*max(counts) #0.5*(max(counts)-min(counts))+min(counts)
# t = plt.text(xtext, ytext,
#              f'tot #={len(data)}\nmean={popt[1]:.2f}\ndev={popt[2]:.2f}\
#              ({popt[2]/abs(popt[1])*100:.1f}%)',
#              ha="center", va="center", size=15,
#              bbox=dict(boxstyle="round,pad=0.3", fc="tab:orange", ec="k", lw=2, alpha=.5))