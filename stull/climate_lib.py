#! /usr/bin/python
###############################################
# atmolib v.0.1 (29.06.2026)
###############################################
import numpy as np
from metpy import constants as c
from metpy.units import units

#------------------------------------------------
# Clausius Clapeyron
#------------------------------------------------
def eq_clapeyron(T):
    '''
    This function calculates the saturation water vapour pressure 
    from the temperature given in Kelvin using the 
    Clausius-Clapeyron equation. The return values are in Pa.
    '''
    R = 8.314 * units.J / (units.K * units.mol)          # J/(Kmol) gas constant
    M = 18.015 * 10 ** (-3) * units.kilogram / units.mol # Kg/mol vapour molar mass 
    e_0 = 0.6113 * 1000 * units.pascal                          # Pa vapour pressure at T = 273.15
    T_0 = 273.15 * units.K                               # K
    L = 2.5 * 10 ** 6 * units.J / units.kilogram         # J/Kg latent heat of vaporization
    R_v = R / M 
    return e_0 * np.exp((L/R_v) * (1/T_0 - 1/T))         # Pa

#------------------------------------------------
# Dry adiabats
#------------------------------------------------
def dry_adiabats():
    '''
    This function computes a set of dry adiabats. The returned values
    are in Kelvin.
    '''
    temperatures = []
    p_0 = 100 * 1000 * units.pascal # Pa 
    p = np.linspace(10, 100) * 1000 * units.pascal # pressure in Pa
    t_0 = np.array([-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60]) * units.degK + c.T0
    for i in range(len(t_0)):
        T = t_0[i] * (p / p_0) ** (c.Rd / c.Cp_d)
        temperatures.append(T)
    return temperatures

#------------------------------------------------
# Isohumes
#------------------------------------------------
def isohumes():
    '''
    This function computes a set of isohumes to be used
    as baselines for a thermodiagram. The returned values
    are in Kelvin.
    '''
    temperatures = []
    p = np.linspace(10, 100) * 1000 * units.pascal # Pa
    e_0 = 0.6113 * 1000 * units.pascal # water vapor pressure in Pa at T_0
    r_s_values = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]) / 1000  # in gg^(-1)
    for r_s in r_s_values: 
        d = 1 / c.T0 - (c.Rv / c.Lv) * np.log(r_s * p / (e_0 * c.epsilon))
        T = 1 / d 
        temperatures.append(T)
    return temperatures

#------------------------------------------------
# Moist adiabatic lapse rate
#------------------------------------------------
def moist_adiabat_lapse_rate(p, T):
    '''
    This function calculates the moist adiabatic lapse rate
    from pressure and temperature. The input temperature is in K
    and the input pressure is in Pa. The return value is in K/Pa
    '''
    e = eq_clapeyron(T)
    r_s = c.epsilon * e / (p - e)
    n = (c.Rd / c.Cp_d) * T + (c.Lv / c.Cp_d) * r_s
    d = (1 + c.Lv ** 2 * r_s * c.epsilon / (c.Cp_d * c.Rd * T ** 2)) * p
    G_s = n / d
    return G_s

#------------------------------------------------
# Moist adiabat
#------------------------------------------------
def moist_adiabat_deprecated(T_0):
    '''
    This function computes the moist adiabat starting from
    the input temperature (Kelvin). The returned values are in Kelvin.
    '''
    num_steps = 50
    ground_p = 100 * 1000 * units.pascal
    top_p = 10 * 1000 * units.pascal
    step_p = (ground_p - top_p) /  num_steps # Pa
    m_adiabat = np.zeros(num_steps) * units.degK
    pressure_levels = np.linspace(ground_p, top_p, num_steps)
    i = 0
    T = T_0
    m_adiabat[i] = T_0
    p = pressure_levels[i]
    for i in range(1, num_steps):
        T = T + moist_adiabat_lapse_rate(p, T) * (- step_p)
        m_adiabat[i] = T
        p = pressure_levels[i]
        
    return m_adiabat

#------------------------------------------------
# Moist adiabats (a set of moist adiabat)
#------------------------------------------------
def moist_adiabats_deprecated():
    '''
    This function computes a set of moist adiabats
    '''
    adiabats = [] 
    t_0 = np.array([-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60]) * units.degK + c.T0
    for i in range(len(t_0)):
        T = moist_adiabat_deprecated(t_0[i])
        adiabats.append(T)
    return adiabats