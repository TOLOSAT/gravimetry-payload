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
import GH_convert      as conv
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
    also at: http://mitgcm.org/~mlosch/geoidcookbook.pdf
    More about WGS84 at: https://earth-info.nga.mil/GandG///wgs84/gravitymod/egm2008/egm08_wgs84.html
    """
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
def ALF_norm_gcb (N, M, phi):
    """
    returns an array[m+1,n+1] of the values of the Associated Legendre Function
    of all integer degrees l and order m, at point x
    Array is normalized, equations from the geoid cook book
    This method is called the "standard forward colums method" explained in:
    https://link.springer.com/article/10.1007/s00190-002-0216-2

    todo:
        compute the derivatives as well
        figure out if it is cos(phi) or sin(phi)
        the thing maxes out at some point
    """
    t = sin(phi)
    u = cos(phi)

    POL = np.zeros((N+1, M+1))
    POL[0,0] = 1
    POL[1,0] = t
    POL[1,1] = u * np.sqrt(3)

    a_nm = lambda n, m : np.sqrt( (2*n+1)*(2*n-1) / ((n-m)*(n+m)) )
    b_nm = lambda n, m : np.sqrt( (2*n+1)*(n+m-1)*(n-m-1) / ((n-m)*(n+m)*(2*n-3)) )

    for n in range(2, M+1):
        POL[n,n] = u*np.sqrt((2*n+1)/(2*n))*POL[n-1,n-1]
#    for m in range(1, M+1):
#        prod_i = 1
#        for i in range (1, m): prod_i = prod_i * np.sqrt( (2*i + 1) / (2*i) )
#        POL[m,m] = u**m * np.sqrt(3) * prod_i

    for n in range (1, N+1): # m
        POL[n, 0] = a_nm(n,0) * t * POL[n-1,0] - b_nm(n,0)*POL[n-2,0]

    for m in range (1, M+1): # m
        for n in range (m+1, N+1): # n
            POL[n, m] = a_nm(n,m)*t*POL[n-1,m] - b_nm(n,m)*POL[n-2,m]

    return POL


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
    Equation obtained from the GCB
    """
    k = 2
    if (m == 0) : k = 1
    P1 = math.factorial(l - m)
    P2 = k*(2*l + 1)
    P3 = math.factorial(l + m)
    N = np.sqrt(P1*P2/P3)
    return N


def Normalize1 (l, m):
    """
    Returns the normalization coefficient of degree l and order m
    Equation obtained from the documentation that came with the EGM2008 coeffs
    thisfunction zeros out below 2e-162
    """
    d_om = 0
    if (m == 0) : d_om = 1
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


def TEST_Normalize():
    """
    This test showed that the normalization coefficient calculated becomes 0.0
    when its (positive) value goes below 2.27e-162. The function must change
    This value is reched for orders above 55
    """
    lmax = 100; Norm_lm = np.zeros((lmax+1,lmax+1))
    for l in range (0,lmax+1):
        for m in range(0,l+1):
            Norm_lm[l,m] = Normalize(l, m)
    return Norm_lm


def TEST_Normalize2():
    """ yes, all zeros are in the output. """
    lmax = 100; Norm_lm = np.zeros((lmax+1,lmax+1))
    for l in range (0,lmax+1):
        for m in range(0,l+1):
            Norm_lm[l,m] = Normalize(l, m) -Normalize1(l, m)
    return Norm_lm


def TEST_APF():
    lmax = 155
    phi = 180/pi
    return ALF_norm_gcb(lmax, lmax, phi*pi/180)





# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

#    TEST_Radius()

#    TEST_Constants()

#    TEST_gravity()

    Nn_lm2 = TEST_Normalize()
#    zeros_ = TEST_Normalize2()

#    a = Normalize(590, 60)
#    b = Normalize(591, 60)
#    c = Normalize(590, 61)


#    alf_mat = TEST_APF()
    aa = ALF_norm_gcb(200, 200, 1)

    Plm_z, Plm_dz = lpmn (10, 10,pi)



    print("\nGH_import done")