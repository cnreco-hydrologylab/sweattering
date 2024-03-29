"""
Library: lib_SOIL.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
import numpy as np


def hallikainen(freq: int, sand: float, clay: float, water: np.array):
    """
    coeff: dict('real':[],'img':[])
    """

    # Part,Frequency(GHz),a0,a1,a2,b0,b1,b2,c0,c1,c2
    coeff = {
        1.4:
            {
                'real': [2.862, -0.012, 0.001, 3.803, 0.462,
                         -0.341, 119.006, -0.500, 0.633],
                'img': [0.356, -0.003, -0.008, 5.507, 0.044,
                        -0.002, 17.753, -0.313, 0.206]
            },
        4:
            {
                'real': [2.927, -0.012, -0.001, 5.505, 0.371,
                         0.062, 114.826, -0.389, -0.547],
                'img': [0.004, 0.001, 0.002, 0.951, 0.005,
                        -0.01, 16.759, 0.192, 0.29]
            },
        6:
            {
                'real': [1.993, 0.002, 0.015, 38.086, -0.176,
                         -0.633, 10.72, 1.256, 1.522],
                'img': [-0.123, 0.002, 0.003, 7.502, -0.058,
                        -0.116, 2.942, 0.452, 0.543]
            }
    }

    coeff_r = np.array(coeff[freq]['real'])
    real = coeff_r[0] + coeff_r[1] * sand + coeff_r[2] * clay + \
           (coeff_r[3] + coeff_r[4] * sand + coeff_r[5] * clay) * water + \
           (coeff_r[6] + coeff_r[7] * sand + coeff_r[8] * clay) * (water ** 2)

    coeff_i = np.array(coeff[freq]['img'])
    img = coeff_i[0] + coeff_i[1] * sand + coeff_i[2] * clay + \
          (coeff_i[3] + coeff_i[4] * sand + coeff_i[5] * clay) * water + \
          (coeff_i[6] + coeff_i[7] * sand + coeff_i[8] * clay) * (water ** 2)

    return real, img


# ----------------------------------------------------------------------------

def doi(freq: float, sand: float, clay: float, water: np.array, angle: float):
    """
    freq [GHz]
    angle [°]
    
    return depth [m]
    """
    c = 299792458  # m/s
    theta = angle * np.pi / 180.  # angle [rad]

    coeff = {
        1.4:
            {
                'real': [2.862, -0.012, 0.001, 3.803, 0.462,
                         -0.341, 119.006, -0.500, 0.633],
                'img': [0.356, -0.003, -0.008, 5.507, 0.044,
                        -0.002, 17.753, -0.313, 0.206]
            },
        4:
            {
                'real': [2.927, -0.012, -0.001, 5.505, 0.371,
                         0.062, 114.826, -0.389, -0.547],
                'img': [0.004, 0.001, 0.002, 0.951, 0.005,
                        -0.01, 16.759, 0.192, 0.29]
            },
        6:
            {
                'real': [1.993, 0.002, 0.015, 38.086, -0.176,
                         -0.633, 10.72, 1.256, 1.522],
                'img': [-0.123, 0.002, 0.003, 7.502, -0.058,
                        -0.116, 2.942, 0.452, 0.543]
            }
    }

    coeff_r = np.array(coeff[freq]['real'])
    real = coeff_r[0] + coeff_r[1] * sand + coeff_r[2] * clay + \
           (coeff_r[3] + coeff_r[4] * sand + coeff_r[5] * clay) * water + \
           (coeff_r[6] + coeff_r[7] * sand + coeff_r[8] * clay) * (water ** 2)

    coeff_i = np.array(coeff[freq]['img'])
    img = coeff_i[0] + coeff_i[1] * sand + coeff_i[2] * clay + \
          (coeff_i[3] + coeff_i[4] * sand + coeff_i[5] * clay) * water + \
          (coeff_i[6] + coeff_i[7] * sand + coeff_i[8] * clay) * (water ** 2)

    depth = c / (2 * np.pi * freq * 1e9) * (np.sqrt(real) / img) * np.cos(
        theta)

    return depth


# ----------------------------------------------------------------------------
def SaxtonRawls(Sand, Clay, OM):
    Sand = Sand / 100
    Clay = Clay / 100

    SM33t = -0.251 * Sand + 0.195 * Clay + 0.011 * OM + 0.006 * Sand * OM - 0.027 * Clay * OM + 0.452 * Sand * Clay + 0.299
    SM33 = SM33t + 1.283 * SM33t ** 2 - 0.374 * SM33t - 0.015

    SM1500t = -0.024 * Sand + 0.487 * Clay + 0.006 * OM + 0.005 * Sand * OM - 0.013 * Clay * OM + 0.068 * Sand * Clay + 0.031
    SM1500 = SM1500t + 0.14 * SM1500t - 0.02

    SMsat33t = 0.278 * Sand + 0.034 * Clay + 0.022 * OM - 0.018 * Sand * OM - 0.027 * Clay * OM - 0.584 * Sand * Clay + 0.078
    SMsat33 = SMsat33t + 0.636 * SMsat33t - 0.107

    lam = (np.log(SM33) - np.log(SM1500)) / (np.log(1500) - np.log(33))

    SMsat = SM33 + SMsat33 - 0.097 * Sand + 0.043

    Ksat = 1930 * (SMsat - SM33) ** (3 - lam)  # [mm/h]
    Ksat = Ksat * 24  # [mm/day]

    return lam, SMsat, Ksat
