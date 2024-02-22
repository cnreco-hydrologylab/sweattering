"""
Library: lib_models_components.py
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

import hydroeval as he

# Graphics
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import IPython.display as disp
get_ipython().run_line_magic('matplotlib', 'inline')

# -----------------------------------------------------------------------------

def calc_depth(freq:float, soil:list, water:list, angle:float)->list:
    depth = [doi(freq=freq, sand=soil[0], clay=soil[1], water=w, angle=angle)*1000 for w in water]
    return depth


def calc_rho(rho_st:float, Kc:list, Kc0:float, EPOT:list)->list:
    rho = [rho_st+0.04*(5-Kc[i]*Kc0*EPOT[i]) for i in range(len(Kc))]
    return rho


def calc_Ks(WW:list, rho:list, WW_fc, WW_w)->list:
    t = len(WW)
    Ks = []
    for i in range(t):
        if WW[i-1]>=((1-rho)*WW_fc+rho*WW_w):
            Ks.append(1)
        elif (WW[i-1]<((1-rho)*WW_fc+rho*WW_w))and(WW[i-1]>WW_w):
            Ks.append((WW[i-1]-WW_w)/((1-rho)*(WW_fc-WW_w)))
        else: Ks.append(0)
    return Ks


def calc_Ks_single(WW:float, rho:float, WW_fc:float, WW_w:float)->float:
    Ks=None
    if WW>=((1-rho)*WW_fc+rho*WW_w):
        Ks=1.
    elif (WW<((1-rho)*WW_fc+rho*WW_w))and(WW>WW_w):
        Ks=((WW-WW_w)/((1-rho)*(WW_fc-WW_w)))
    else: Ks=0.
    return Ks