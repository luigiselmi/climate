#! /usr/bin/python
###############################################
# atmolib v.0.1 (29.06.2026)
###############################################
import numpy as np

def eq_clapeyron(T):
    '''
    This function calculates the saturation water vapour pressure 
    from the temperature given in Kelvin using the 
    Clausius-Clapeyron equation. The return values are in kPa.
    '''
    R = 8.314  # J/(Kmol) gas constant
    M = 18.015 * 10 ** (-3) # Kg/mol vapour molar mass 
    p_0 = 0.6113 # kPa vapour pressure at T = 273.15
    T_0 = 273.15 # K
    L = 2.5 * 10 ** 6 # J/Kg latent heat of vaporization
    R_v = R / M 
    return p_0 * np.exp((L/R_v) * (1/T_0 - 1/T)) # kPa