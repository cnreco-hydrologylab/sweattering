"""
Library: lib_generic.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
#!/usr/bin/env python
# coding: utf-8

# Ausiliary functions (your @staticmethods!)
# To import all the functions in this module, run:
#
# from dir.file import *
#
# being dir the directory of file.

# Base
import os
import re
import time
import math
import numpy as np
import pandas as pd
import datetime as dtt

# Analysis
import pyswarms as ps
from scipy import special as sp
from scipy.stats import norm, gamma, f, chi2
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter as sfilter

# Geospatial
import hydroeval as he

# Graphics
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import IPython.display as disp
get_ipython().run_line_magic('matplotlib', 'inline')

#----------------------------------------------------------------------------# Statistics and data cleaning / normalization

def lin_db(x):
    """linear to dB"""
    return 10*np.log10(x)


def db_lin(x):
    """dB to linear"""
    return 10**(x/10)


def norm_simple(x):
    """min-max normalization of x data"""
    return (x-np.min(x))/(np.max(x)-np.min(x))


def norm_fit(x, a, b):
    """norm in range a,b
    
    x = data, a<b
    """
    return a+(x-np.nanmin(x))*(b-a)/(np.nanmax(x)-np.nanmin(x))


def bias(obs, sim):
    """distance between obs' and sim's mean values"""
    if len(obs)==len(sim):
        return np.mean(obs-sim)
    else: raise ValueError(
        f'obs and sim must have same first dimension, but have shapes {np.shape(obs)} and {np.shape(sim)}')

    
def Rvalue(x:list,y:list)->float:
    """compute Pearson's R between x,y data"""
    if len(x)!=len(y):
        raise ValueError(
            f'x and y must have same first dimension,'+
            'but have shapes{np.shape(x)} and {np.shape(y)}')
        
    mask = ~(np.isnan(x) | np.isnan(y))
    x_filtered = np.array(x)[mask]
    y_filtered = np.array(y)[mask]
    xy_filtered = np.transpose(np.column_stack((x_filtered, y_filtered)))
    return np.corrcoef(xy_filtered)[0][1]


def timeseries(dates, data):
    """Returns a matrix (list of type(dates,data)) in the format [dates,data]"""
    
    if len(dates)==len(data):
        return [[dates[i],data[i]] for i in range(len(dates))]
    else: raise ValueError(
        f'dates and data must have same first dimension, but have shapes {np.shape(dates)} and {np.shape(data)}')
    

def significant_figures_str(master, slave):
    master_rounded = float("%.1g" % master)
    master_rounded_str = str(master_rounded)
    digits = master_rounded_str[::-1].find('.')
    slave_rounded_str = str(f'%.{digits}f' % slave)
    
    return [master_rounded_str, slave_rounded_str]


def substitute_keywords(template, **kwargs):
    for key, value in kwargs.items():
        template = template.replace('{' + key + '}', str(value))
    return template
    

#----------------------------------------------------------------------------
# Data analysis, fit
# Fitting functions


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

    
def HIST_norm(ref_mean, ref_std, obs:list):
    """HIST normalization
    Ref. Mladenova, 2013, https://ieeexplore.ieee.org/document/6264094
    
    obs = [value, mean, std]
    """
    value, mean, std = obs
    return ref_mean+ref_std/std*(value-mean)
