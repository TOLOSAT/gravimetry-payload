# -*- coding: utf-8 -*-
"""
    Global variables for the whole project, derived from WGS-84.
    Figures are in SI units.

    WIP 
"""

global omegaEarth
global radiusEarth 
global muEarth

omegaEarth = 7292115E-11
radiusEarth = 6378136
muEarth = 3.986004415E05
G = 6.673E-11 # m^3/s^2/kg : Gravitational constant
g = 9.80665 # m/s^2 : average surface acceleration
ro = 5515 # kg/m^3 : earth average density
W_0 = 62636856.0 # m^2/s^2 : average geopotential (US value, must find better)
W_GH = 62601662.83663869 # m^2/s^2 : average geopotential at topographic surface calculated using GH3 (in GH_harmonics)


# WGS84 reference ellipsoid model
ref_e = "WGS 84" # https://en.wikipedia.org/wiki/World_Geodetic_System
GM_e = 3986004.418E8 # m^3/s^2 : standard gravitational parameter
wo = 7292115E-11 # rad/s : angular velocity of Earth
a_e = 6378137.00 # m : equatorial radius or ellipsoid model
f = 1/298.257223563 # flat parameter

b_e = a_e * (1-f) # m : polar radius
E = np.sqrt(a_e**2 - b_e**2) # linear eccentricity
e_1 = E/a_e
e_2 = E/b_e
m = wo**2 * a_e*2 * b_e / GM_e # just to simplify the code
g_a = GM_e/(a_e*b_e) * (1 - 3/2*m - 3/14*e_2*m) # m/s^2 : gravity acc. at equator
g_b = GM_e/(a_e**2) * (1 - m - 3/7*e_2*m) # m/s^2 : gravity acc. at poles

# EGS2008 potential model
ref_g = "EGS2008"
a_g = 6378136.3 # m : Reference radius for the potential model
GM_g = 3986004.415E8 # m^3/s^2 : standard gravitational parameter in the potential model
