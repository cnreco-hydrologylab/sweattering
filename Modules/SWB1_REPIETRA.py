"""
Library: SWB1_REPIETRA.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
#!/usr/bin/env python
# coding: utf-8

from Modules.lib_models_components import *

# ----------------------------------------------------------------------------

class SWB1:
    
    def __init__(self, PAR_str, inputs, user_in):
        self.inputs = inputs
        self.user_in = user_in
        self.PAR_str = PAR_str
    
    
    def model(self, PAR, PAR_str, inputs, user_in):
        """Soil Water Balance 1 layer (bucket model)

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

        # Inputs
        t, mask_gap_sm, P, EPOT, sm, rho_st = inputs
        a, Kc0, n1, Ksat1, depth, WW_fc, WW_w = PAR

        WW1   = [sm[0]]  # water content I layer [m3/m3]
        IE    = [0]
        PERC1 = [0]
        ET1   = [0]

        for i in [i + 1 for i in range(len(t) - 1)]:

            # Infiltration excess (only for I layer)
            IE.append(P[i] * (WW1[i - 1] / WW_fc) ** a)  # [mm]

            # Evapotranspiration
            rho = rho_st
            Kstress = calc_Ks_single(WW1[i-1], rho, WW_fc, WW_w)
            ET1.append(EPOT[i]*(WW1[ i - 1 ] / WW_fc)*Kc0*Kstress) # assuming static Kc=1

            # Percolation
            PERC1_base = Ksat1 * (WW1[i - 1] / WW_fc) ** n1
            PERC1_amm = WW1[i - 1] * depth + P[i] - ET1[i] - IE[i]
            if PERC1_amm < 0: PERC1_amm = 0
            PERC1_result = min(
                [PERC1_base, PERC1_amm]) if PERC1_base > 0 and PERC1_amm > 0 else 0
            PERC1.append(PERC1_result)

            # Water balance I layer [mm]
            WW1.append(WW1[i - 1] + (P[i] - IE[i] - PERC1[i] - ET1[i]) / depth)

            # Computation of deep percolation (water above field capacity)
            if WW1[i] < WW_w: WW1[i]  = WW_w
            if WW1[i] > WW_fc: WW1[i] = WW_fc
            

        WW1=np.array(WW1); PERC1=np.array(PERC1); IE=np.array(IE); ET1=np.array(ET1)
        WW1_sm = WW1[mask_gap_sm]
        KGE = he.evaluator(he.kge, WW1_sm, sm)[0][0]
        if np.isnan(KGE): raise ValueError

        return [WW1, PERC1, IE, ET1, KGE]


    def pso_calib_irri(self, PAR):
        """Ausiliary function for PSO optimization"""
        n_particles = PAR.shape[0]
        err = np.zeros(n_particles)
        for i in range(n_particles):
            KGE = self.model(PAR[i], self.PAR_str, self.inputs, self.user_in)[-1]
            err[i] = 1 - KGE
        return err