#! /usr/bin/python
###############################################
# atmolib v.0.1 (29.06.2026)
###############################################
import numpy as np
from metpy import constants as c
from metpy.units import units

def eq_clapeyron(T):
    '''
    This function calculates the saturation water vapour pressure 
    from the temperature given in Kelvin using the 
    Clausius-Clapeyron equation. The return values are in Pa.
    '''
    R = 8.314 * units.J / (units.K * units.mol)          # J/(Kmol) gas constant
    M = 18.015 * 10 ** (-3) * units.kilogram / units.mol # Kg/mol vapour molar mass 
    p_0 = 0.6113 * 1000 * units.pascal                          # Pa vapour pressure at T = 273.15
    T_0 = 273.15 * units.K                               # K
    L = 2.5 * 10 ** 6 * units.J / units.kilogram         # J/Kg latent heat of vaporization
    R_v = R / M 
    return p_0 * np.exp((L/R_v) * (1/T_0 - 1/T))         # Pa

def p_dry_adiabats():
    '''
    This function computes a set of dry adiabats
    '''
    temperatures = []
    p_0 = 100 # kPa 
    p = np.linspace(10, 100) # pressure in kPa
    t_0 = np.array([-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60]) + 273.15
    for i in range(len(t_0)):
        T = t_0[i] * (p / p_0) ** (c.Rd / c.Cp_d)
        temperatures.append(T - 273.15)
    return temperatures

def p_isohumes():
    '''
    This function computes a set of isohumes to be used
    as baselines for a thermodiagram.
    '''
    temperatures = []
    p = np.linspace(10, 100)
    e_0 = 0.6113 # water vapor pressure at T_0
    r_s_values = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]) / 1000  # in gg^(-1)
    for r_s in r_s_values: 
        d = 1 / c.T0 - (c.Rv / c.Lv) * np.log(r_s * p / (e_0 * c.epsilon))
        T = 1 / d - c.T0
        temperatures.append(T)
    return temperatures


def moist_adiabat_lapse_rate(p, T):
    '''
    This function calculates the moist adiabatic lapse rate
    from pressure and temperature.
    '''
    e = eq_clapeyron(T)
    r_s = c.epsilon * e / (p - e)
    n = (c.Rd / c.Cp_d) * T + (c.Lv / c.Cp_d) * r_s
    d = (1 + c.Lv ** 2 * r_s * c.epsilon / (c.Cp_d * c.Rd * T ** 2)) * p
    G_s = n / d
    return G_s
    
def moist_adiabat(T_0):
    '''
    This function calculates the moist adiabat starting 
    from the temperature at the reference level of 100 kPa
    '''
    p_ground = 100 # kPa
    p_top = 10 # kPa
    num_steps = 450
    step_p = 0.2 # kPa
    m_adiabat = np.zeros(450)
    pressure_levels = np.flip(np.linspace(p_top, p_ground, 450))
    i = 0
    T = T_0 + c.T0
    m_adiabat[i] = T_0
    p = pressure_levels[i]
    for i in range(1, num_steps):
        T = T + moist_adiabat_lapse_rate(p, T) * (- step_p)
        m_adiabat[i] = T - 273.15
        p = pressure_levels[i]
        
    return np.flip(m_adiabat)

def p_moist_adiabats():
    '''
    This function computes a set of moist adiabats
    '''
    adiabats = [] 
    t_0 = np.array([-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60])
    for i in range(len(t_0)):
        T = moist_adiabat(t_0[i])
        adiabats.append(T)
    return adiabats