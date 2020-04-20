"""

@authors:

# =============================================================================
 Information:

    The functions in this script are used to calculate values about 
    Earth, Geophysics, and mathematics

# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from numpy import pi, sin, cos
from scipy.special import lpmn
import math
import matplotlib.pyplot as plt
from time import gmtime, strftime

#import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


from GH_convert import cart2sphA

# =============================================================================
# GLOBAL VARIABLES
# =============================================================================
class Constants:
    """
    A variable with all the constants inside it
    Many equations come from the geoid cook book pdf found in documentation, 
    also at:http://mitgcm.org/~mlosch/geoidcookbook.pdf
    """   
    G = 6.673E-11 # m^3/s^2/kg : Gravitational constant
    GM_g = 3986004.415E8 # m^3/s^2 : standard gravitational parameter for the EGM2008 model
    a_g = 6378136.3 # m : Reference radius
    g = 9.80665 # m/s^2 : average surface acceleration
    ro = 5515 # kg/m^3 : earth average density 
    w = 0 # rad/s : angular velocity of Earth
    a = 6378137 # m : equatorial radius
    f = 1/298.257223563 # flat parameter
    b = a * (1-f) # m : polar radius
    E = np.sqrt(a**2-b**2) # linear _eccentricity
    e_1 = E/a
    e_2 = E/b
    m = w**2 * a*2 * b / GM_g # just to simplify the code
    g_a = GM_g/(a*b) * (1 - 3/2*m - 3/14*e_2*m) # m/s^2 : gravity acc. at equator
    g_b = GM_g/(a**2) * (1 - m - 3/7*e_2*m) # m/s^2 : gravity acc. at poles
    W_0 = 62636856.0 # m^2/s^2 : average geopotential (US value, must find better)
    
cst = Constants() 



# =============================================================================
# FUNCTIONS - GEODESY
# =============================================================================
def Get_Ellipsoid_Radius (Lat):
    """
    Returns the radius of the reference elipsoid in meters
    https://gis.stackexchange.com/questions/20200/how-do-you-compute-the-earths-radius-at-a-given-geodetic-latitude
    
    Input:
        Lat: latitude, inclination from the z axis in radians
    Output:
        R: Radius in meters
    """
    a = 6378137 # m : equatorial radius
    f = 1/298.257223563 # flat parameter
    b = a * (1-f) # m : polar radius
    
    Lat = pi/2 - Lat # THIS IS TEMPORARY AND MUST CHANGE
    
    deno = np.sqrt(a**2*sin(Lat)**2 + b**2*cos(Lat)**2)
    R = a*b/deno
    
#    numer = (a**2*cos(Lat))**2 + (b**2*sin(Lat))**2
#    denom = (a   *cos(Lat))**2 + (b   *sin(Lat))**2
#    R = np.sqrt(numer/denom)
    
#    R = a*(1-f*sin(Lat)**2)
    
    return R


def Get_Normal_Gravity (Lat):
    """
    Returns the normal gravity, on the ellipsoid, at the given latitude
    This equaion is made up, until the Get_Normal_Gravity2 function starts
    making sense
    """
    c = Constants()
    g_a=c.g_a; g_b=c.g_b 
    
    f = (g_a - g_b)/g_a    
    g_0_lat = g_a*(1-f*cos(Lat)**2)
    
    return g_0_lat

def Get_Normal_Gravity2 (Lat):
    """
    Returns the normal gravity, on the ellipsoid, at the given latitude
    This equaion is BUUUUUSTED
    """
    c = Constants()
    a=c.a; b=c.b; g_a=c.g_a; g_b=c.g_b; e=c.e_1
    
    k = (b*g_b - a*g_a) / a*g_a
    g_0_lat = g_a * (1 + k**2 * sin(Lat)**2) / np.sqrt(1 - e**2 * sin(Lat)**2)
    
    return g_0_lat

# =============================================================================
# FUNCTIONS - MATHEMATICAL VALUES
# =============================================================================
def Pol_Legendre (l, m, x):
    """
    returns an array[m+1,n+1] of the values of the associated Legendre function
    of all integer degrees l and order m, at point x
    """
    Plm_z, Plm_dz = lpmn(m, l, x)
    return Plm_z, Plm_dz # use Plm_z[m, l]


def Normalize (l, m):
    """
    Returns the normalization coefficient of degree l and order m
    Equation obtained fron the dicumentation that came with the EGM2008 coeffs
    """
    d_om = 0
    if (m == 0) :
        d_om = 1

    P1 = math.factorial(l - m)
    P2 = (2*l + 1)
    P3 = (2 - d_om)
    P4 = math.factorial(l + m)

    N = np.sqrt(P1*P2*P3/P4)
    return N


def dichotomy_grad (f, arg_before, z_e, arg_after, w_0, de, grad):
    """
    f : function for f(*in_first,z,*in_after) = w
    arg_before, arg_after : f function arguments that come before and after z_e
    z_e : firstguess for the value to be tested
    w_0: target value for w
    de : delta error  
    grad : gradient
    """
    w_i = f(*arg_before, z_e, *arg_after)
    di = w_0 - w_i
    z_i = z_e    
    c = 0
#    print(f"dicho start\nw_0={w_0:.2f}; z_i={z_i:.2f}; w_i={w_i:.2f}; di={di:.2f}; add={(di/grad):.2f}; "); 
    while (abs(di) >= de):
        c+=1        
        z_i += di/grad
        w_i = f(*arg_before, z_i, *arg_after)
        di = w_0 - w_i        
#        sleep(1); 
#        print(f"w_0={w_0:.2f}; z_i={z_i:.2f}; w_i={w_i:.2f}; di={di:.2f}; add={(di/grad):.2f}; "); 
#    print(f"dichotomy_grad: {c} steps")
    return z_i



# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_Radius ():
    Lats = np.arange(-90, 91, 10)
    Rads = Get_Ellipsoid_Radius(Lats*pi/180)
    plt.plot(Lats, Rads)


def TEST_Constants():
    cts = Constants()
    print(f"g = {cts.g} m/s^2")
    f = Get_Normal_Gravity(pi/180 * 50)
    print(f"g_0 at Lat = 50 is {f}")
    
    
def TEST_gravity ():
    Lats = np.arange(-90, 91, 10)
    Grav1 = Get_Normal_Gravity(Lats*pi/180)    
    Grav2 = Get_Normal_Gravity2(Lats*pi/180)
    plt.figure()
    plt.clf()
    plt.title("Gravity acceleration (m/s^2) vs lattitute")
    plt.plot(Lats, Grav1)
    plt.plot(Lats, Grav2)



# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    
#    TEST_Radius()
    
#    TEST_Constants()
    
    TEST_gravity()
    
    print("\nGH_import done")

