"""
Library: lib_EPOT.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
# EPOT

import datetime
import sys

# sys.path.append('../')
from pyeto import *

import numpy as np

def hamon(tavg, jdate, lat, par=1.2):
    # ' @title Hamon Potential Evapotranspiration Equation
    # ' @description The Hamon method is also considered as one of the simplest estimates
    # '     of potential Evapotranspiration.
    # ' @param par proportionality coefficient (unitless),
    # PAR=1.2 as in https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1752-1688.2005.tb03759.x
    # ' @param tavg vector of mean daily temperature (deg C)
    # ' @param lat latitude ()
    # ' @param jdate a day number of the year (julian day of the year)  timetj = df.index.day
    # ' @return outputs potential evapotranspiration (mm day-1)
    # ' @details For details see Haith and Shoemaker (1987) 

    var_theta = 0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.0086 * (jdate-186)))
    var_pi = np.arcsin(0.39795 * np.cos(var_theta))
    daylighthr = 24 - 24 / np.pi * np.arccos((np.sin(0.8333 * np.pi / 180) + np.sin(lat * np.pi / 180) * np.sin(var_pi)) / (np.cos(lat * np.pi / 180) * np.cos(var_pi)))
    esat = 0.611 * np.exp(17.27 * tavg / (237.3 + tavg))
    
    return par * 29.8 * daylighthr * (esat / (tavg + 273.2))


# Ref. https://pyeto.readthedocs.io/en/latest/hargreaves.html


def hargre(lat_deg, dates, temp_min, temp_max, temp_mean):
    """Hargreaves-Samani model for ET0 estimation from temperature input.

    Params
    ------
    - lat_deg: float
    - dates: timestamp
    - temp_*: float


    """
    lat = deg2rad(lat_deg)  # Convert latitude in degrees to radians
    day_of_year = dates.dayofyear
    sol_decli = sol_dec(day_of_year)  # Solar declination
    sha = sunset_hour_angle(lat, sol_decli)
    ird = inv_rel_dist_earth_sun(day_of_year)
    et_radia = et_rad(lat, sol_decli, sha, ird)  # Extraterrestrial radiation
    return hargreaves(temp_min, temp_max, temp_mean, et_radia)

## Example use with pandas.DataFrame df
## Final dataframe is made by custom function timeseries

# lat_deg = 44.570842547510622 # latitude of Budrio (deg)
# temp_min = df.min()['Temperatura[°C]'].values
# temp_max = df.max()['Temperatura[°C]'].values
# temp_mean = df.mean()['Temperatura[°C]'].values
# dates = df.asfreq().index
# eto = timeseries( dates,
#                  [ hargre(lat_deg, dates[i] , temp_min[i], temp_max[i], temp_mean[i])
#                   for i in range(len(dates)) ] )
# eto_df = pd.DataFrame(eto).rename(columns={0:'Date',1:'EPOT'}).set_index('Date')