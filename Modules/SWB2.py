"""
Library: SWB2.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
#!/usr/bin/env python
# coding: utf-8


# Graphics

from Modules.lib_models_components import *


#############################################################################
# Soil Water Balance
# No irrigation
# Percolation, infiltration
#############################################################################

class SWB2:

    def __init__(self, inputs, user_in, PAR):
        self.inputs = inputs
        self.user_in = user_in
        self.PAR = PAR


    def model(self):
        """Soil Water Balance

        The soil water balance model (SWB) produces an estimate of the soil water
        content WW [%] that is used to simulate $\sigma^0$ by a water cloud
        model (WCM).

        Inputs
        ----------
        - PAR: initial guess values for parameters to calibrate
        - inputs: input quantities for calibration

        Return
        -------
        KGE from hydroeval between sigma0 observed and simulated.

        """

        # User_in options
        # opt_static = user_in # deprecated, here for compatibility

        # Inputs
        a, n1, n2, Ksat1, Ksat2 = self.PAR
        t, t_sat, P, EPOT, Kc, veg, angle, s0, freq, rho_st, soil, WW_sat, WW_fc, WW_w, sm = self.inputs

        # Depth and soil parameters
        fixed_depth = 100
        sand, clay, om = soil  # soil texture

        Kc0 = 1

        # Saturated hydraulic conductivity and saturation
        # lam1, WWsat1, Ksat1 = SaxtonRawls(sand, clay, om)
        # lam2, WWsat2, Ksat2 = SaxtonRawls(sand, clay, om)
        # n1 = 3+2/lam1
        # n2 = 3+2/lam2

        WWsat1 = WW_sat;
        WWsat2 = WW_sat

        angle_m = np.mean(angle)
        Ks = [0.]  # water stress coefficient
        rho = [0.]  # depletion fraction
        IE = [0.]
        WW1 = [.2]  # water content I layer [m3/m3]
        WW2 = [.4]  # water content II layer [m3/m3]
        PERC1 = [0.]; PERC2 = [0.]
        EPOT1 = [0.]; EPOT2 = [0.]
        ET1 = [0.]; ET2 = [0.]

        COST = .0  # additional cost to KGE
        LAMBDA = 1000  # Lagrange multiplier

        for i in [i + 1 for i in range(len(t) - 1)]:

            # Compute depth of layers
            WW1_before = WW1[i - 1]
            depth1 = doi(freq=freq,
                         sand=sand, clay=clay,
                         water=WW1_before,
                         angle=angle_m) * 1000 \
                     + 50  # *1000 accounts for going from [m] to [mm]
            # 2 cm are added to improve stability
            depth2 = 100  # fixed_depth-depth1+20 # fixed max depth - depth [mm]

            # Infiltration excess (only for I layer)
            IE.append(P[i] * (WW1[i - 1] / WWsat1) ** a)  # [mm]

            # Evapotranspiration
            rho = rho_st
            WW_mean = (WW1[i - 1] * depth1 + WW2[i - 1] * depth2) / fixed_depth
            Kstress = calc_Ks_single(WW_mean, rho, WW_fc, WW_w)
            ET = EPOT[i] * Kc[i] * Kc0 * Kstress
            ET1.append(ET * depth1 / fixed_depth)
            ET2.append(ET * depth2 / fixed_depth)

            # percolation [mm/day]
            PERC1_base = Ksat1 * (WW1[i - 1] / WWsat1) ** n1 if WW2[
                                                                    i - 1] < WWsat2 else 0  # [mm/day]
            PERC1_amm = WW1[i - 1] * depth1 + P[i] - ET1[i] - IE[i]
            if PERC1_amm < 0: PERC1_amm = 0
            PERC1_result = min(
                [PERC1_base, PERC1_amm]) if PERC1_base > 0 and PERC1_amm > 0 else 0
            PERC1.append(PERC1_result)

            PERC2_base = Ksat2 * (WW2[i - 1] / WWsat2) ** n2  # [mm/day]
            PERC2_amm = WW2[i - 1] * depth2 + PERC1[i] - ET2[i]
            if PERC2_amm < 0: PERC2_amm = 0
            PERC2_result = min(
                [PERC2_base, PERC2_amm]) if PERC2_base > 0 and PERC2_amm > 0 else 0
            PERC2.append(PERC2_result)

            # Water balance I layer [mm]
            WW1.append(WW1[i - 1] + (P[i] - IE[i] - PERC1[i] - ET1[i]) / depth1)
            WW2.append(WW2[i - 1] + (PERC1[i] - PERC2[i] - ET2[i]) / depth2)

            # Computation of deep percolation (water above field capacity)
            # if WW1[i]>WW_fc: WW1[i]=WW_fc
            # if WW2[i]>WW_fc: WW2[i]=WW_fc
            if WW1[i] < WW_w: WW1[i] = WW_w
            if WW2[i] < WW_w: WW2[i] = WW_w
            if WW1[i] > WWsat1: WW1[i] = WWsat1
            if WW2[i] > WWsat2: WW2[i] = WWsat2

            # Regularization for non-physical data below wilting point
            # if WW[i]<WW_w: COST += (WW[i]-WW_w)**2

        KGE = he.evaluator(he.kge, WW2, sm)[0][0]

        return [WW1, WW2, PERC1, PERC2, IE, ET1, ET2, KGE]


    def pso_calib_irri(self):
        """Ausiliary function for PSO optimization"""
        n_particles = self.PAR.shape[0]
        err = np.zeros(n_particles)
        for i in range(n_particles):
            KGE = self.model(self.PAR[i], self.inputs, self.user_in)[-1]
            err[i] = 1 - KGE
        return err