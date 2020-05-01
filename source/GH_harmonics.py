"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to calculate the sums from spherical 
    harmonic coefficients, and generate grids of data mapped out over the 
    surface of the Earth
    
    Generally used variables:
        R, Lat, Long = coordinates in the geocentric CRS
        lmax = maximum degree to calculate to for the geopotantial
        HC, HS = Geopotential stokes coefficients
        lmax_topo, HC_topo, HS_topo = same, but for topography. 
        
        limits = [Western_long, Eastern_long, Southern_lat, Northern_lat]
        tens = how large the grid should be - look at init_grid() to understand
        
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import matplotlib.pyplot as plt 
import numpy as np
from numpy import pi, sin, cos
from time import sleep

import GH_import       as imp
import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
import GH_terminal     as term
#import GH_harmonics    as harm
import GH_geoMath      as gmath
#import GH_earthMap     as emap



# =============================================================================
# FUNCTIONS TO GENERATE DATA ARRAYs
# =============================================================================
def init_grid (tens, limits):
    """
    Initiates the grid variables based on the number of points wanted
    within the given limits
    """
    size_long = 1 + 36*tens
    size_lat  = 1 + 18*tens
    dim = limits * pi/180
    
    Line_theta = np.linspace(dim[0], dim[1], size_long)
    Line_phi  = np.linspace(dim[2], dim[3], size_lat)
    
    G_theta, G_phi = np.meshgrid(Line_theta, Line_phi)
    G_Grid = np.zeros((size_lat, size_long))
    
    return G_Grid, G_theta, G_phi


def Gen_Grid (tens, Get_FUNCTION, in_args, limits):
    """
    This function generates a grid of the desired spherical harmonic model
    at Lat/Long coordinates
    Input: 
        tens: how large the array should be
        Get_FUNCTION: the callable function that must be used
        *in_args: the arguments to the callable function besides R, Lat, Long
    Output:
        G_Grid: grid of Get_FUNCTION(R,Lat,Long,*in_args)
        G_Long: grid of longitudes
        G_Lat:  grid of latitudes
    """
    G_Grid, G_Long, G_Lat = init_grid(tens, limits)   
    print(f"Making a grid with \"{Get_FUNCTION.__name__}()\", with {G_Grid.size} points\n",end="\r")

#    if (measure == "geopot"):
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    
    for j in range(0, G_Lat.shape[0]):
        term.printProgressBar(j+1, G_Lat.shape[0])
        Lat = pi/2 - G_Lat[j][0]
        R = gmath.Get_Ellipsoid_Radius(Lat)
        
        for i in range(0, G_Long.shape[1]):
            Long = G_Long[0][i] - pi
            print(f"\rLong =  {Long} ;Lat {Lat}",end="\r")
            G_Grid[j,i] = Get_FUNCTION(R, Lat, Long, *in_args)
          
    return G_Grid, G_Long*180/pi, G_Lat*180/pi # in degrees now



# =============================================================================
# FUNCTIONS TO CALCULATE SPHERICAL HARMONIC SUMS
# =============================================================================
def Get_Topo_Height (R, Lat, Long,    lmax_topo, HC_topo, HS_topo):
    """
    This function returns the height of Earth's estimated topography at Lat/Long 
    coordinates
    The solution is calculated up to degree lmax in the [HC, HS] model
    """
    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax_topo, lmax_topo, cos(Lat)) # I am allowed to write that.
    for l in range (0, lmax_topo+1):
        Sum2 = 0
        for m in range (0, l+1):
#            print(f"\rl={l} ; m={m}",end="\r")
            Sum2 += (HC_topo[l,m]*cos(m*Long) + HS_topo[l,m]*sin(m*Long)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum1 += Sum2

    return Sum1


def Get_Geo_Pot (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    """
    cst = gmath.Constants()
    
    R+= Get_Topo_Height (R, Lat, Long,    lmax_topo, HC_topo, HS_topo)
    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(Lat))
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum1 += (cst.a_g/R)**l * Sum2

    geopot = cst.GM_g/R*(1 + Sum1)

    return geopot


def Get_Geoid_Height (R, Lat, Long,    lmax, HC, HS):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The cosine coefficients for even l and m=0 are corrected to remove the 
    reference ellipsoid from the results
    
    Equations from the geoid cook book
    """
    cst = gmath.Constants()    
    g_0 = gmath.Get_Normal_Gravity(Lat)
    Lat_gc = conv.geodes2geocen(Lat)
    
#    R = cst.a_g 
    g_0 = cst.g 
    Lat_gc = Lat
    
    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, sin(Lat_gc) )
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, 1):#═l+1):            
            HC_lm = CorrCos_lm(l, m, HC[l,m])            
            Sum2 += (HC_lm*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * gmath.Normalize(l, m)
            
        Sum1 +=  (cst.a_g/R)**l * Sum2 

    Geo_H = cst.GM_g * Sum1 / (R*g_0)

    return Geo_H


def Get_Geoid_Height2 (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    
    Equations from the GFZ handbook, eq.116
    """
    c = gmath.Constants() 
    a_g=c.a_g; GM_g=c.GM_g; # g_e=c.g 
    G=c.G; ro=c.ro;  
    
    R_e = gmath.Get_Ellipsoid_Radius(Lat)
    g_e = gmath.Get_Normal_Gravity(Lat)
#    Lat_gc = conv.geodes2geocen(Lat)
    Lat_gc = Lat
    
    Sum_geo = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, sin(Lat_gc) )
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long)+HS[l,m]*sin(m*Long)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum_geo +=  (a_g/R_e)**l * Sum2 

    Sum_topo = 0
    P_lm, _ = gmath.Pol_Legendre(lmax_topo, lmax_topo, cos(Lat))
    for l in range (0, lmax_topo+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC_topo[l,m]*cos(m*Long) + HS_topo[l,m]*sin(m*Long)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum_topo += Sum2

    Geo_H = GM_g/(R_e*g_e) * Sum_geo  -  2*pi*G*ro/g_e * (R_e*Sum_topo)**2

    return Geo_H


def Get_acceleration (R, Lat, Long,    lmax, HC, HS):
    """
    This function returns the acceleration at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The acceleration is the gratient of the geopotential, and is calculated 
    over a distance "d" 
    """
    c = gmath.Constants()
    a_g=c.a_g; GM_g=c.GM_g; 
    
    d = 1 # m
    R_1 = R - d/2
    R_2 = R + d/2
    
    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(Lat))
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum1 += (a_g/R)**l * Sum2

    PG_1 = GM_g/R_1*(1 + Sum1)
    PG_2 = GM_g/R_2*(1 + Sum1)
    
    local_acc = (PG_1 - PG_2) / d
    
    return local_acc


def Get_isopot (R, Lat, Long,    W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo): 
    """
    This function returns the Height at given Lat/Long coordinates at which the 
    geopotential is equal to the given W_0
    The solution is calculated up to degree lmax in the HC HS model
    The approach is a dichotomy within a width of "e" m error
    """    
    de = 10
    grad = -10 # 9.81
    
    arg_after = [Lat, Long, lmax, HC, HS, lmax_topo, HC_topo, HS_topo]
    R_iso = gmath.dichotomy_grad (Get_Geo_Pot, [], R, arg_after, W_0, de, grad)
    Height = R_iso - R
    return Height # , R_iso 



# =============================================================================
# SUB FUNCTIONS BUT STILL HARMONICS
# =============================================================================
def CorrCos_lm (l, m, HC_lm):
    """ correct the HC_lm coef if the conditions are correct - GCB """
    lmax_corr=10
    if ( (m==0) and (l<=lmax_corr) and (l%2 == 0) ):
        HC_lm += Cosine_Correction(l)
    return HC_lm

def Cosine_Correction (N): 
    """ 
    Returns the corrected Cosine coefficient 
    The maths in this book cmes from the Geoid cook book and the fortran code 
    I dissasembled
    """    
    n = N/2
    
    c=gmath.Constants()    
    GM = c.GM_e
    GM_g = c.GM_e
    a = c.a_g
    a_g = c.a_g
    E = c.E 
    m = c.m
    e_2 = c.e_2
    
    q0 = 1/2*((1+3/e_2**2) * np.arctan(e_2) - 3/e_2)    
    
    CA_ME2 = 1/3*(1 - 2/15*m*e_2/q0)
    
    J_2n = (-1)**(n+1) * 3*(E/a)**(2*n) / ((2*n+1)*(2*n+3))
    J_2n = J_2n * (1- n + 5*n*CA_ME2)
    
    Corr = J_2n / np.sqrt(4*n + 1) * GM/GM_g * (a/a_g)**n

    return Corr 


def Cosine_Correction2 (N):
    """ raw data from the hsynth fortran code """    
    C = np.zeros(21)
    C[ 2] =  -0.484166774985E-03                            
    C[ 4] =   0.790303733524E-06                            
    C[ 6] =  -0.168724961164E-08                            
    C[ 8] =   0.346052468485E-11                            
    C[10] =  -0.265002226377E-14                            
    C[12] =  -0.410790140988E-16                            
    C[14] =   0.447177356743E-18                            
    C[16] =  -0.346362564558E-20                            
    C[18] =   0.241145603096E-22                            
    C[20] =  -0.160243292770E-24 
    return C[N]
    
    
    

# =============================================================================
# TEST FUNCTIONS
# =============================================================================    
def Math_calc_geopot_basic(z):
    """ some function needed in TEST_plot_radius """
    G = 6.673E-11
    M = 5.975E24
    a = 6.378E6    
    P = G*M*(1/a + 1/(a+float(z)))
    return P

def TEST_plot_radius(fignum, lmax, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    Plots the geopotential at given coordinates from Earth's center to the 
    surface, and beyond into space.
    """    
    R_earth = gmath.Get_Ellipsoid_Radius(Lat)
    Rs = np.linspace(95, 105, 100)
    G_Pot = np.zeros(Rs.shape)
    G_Pot_Basic = np.zeros(Rs.shape)
    
    for i in range (len(Rs)):
        G_Pot[i] = Get_Geo_Pot (Rs[i]*R_earth, Lat*pi/180, Long*pi/180, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)
        G_Pot_Basic[i] = Get_Geo_Pot (Rs[i]*R_earth, Lat*pi/180, Long*pi/180, 2, HC, HS, lmax_topo, HC_topo, HS_topo)
#        G_Pot_Basic[i] = Math_calc_geopot_basic(Rs[i]*R_earth)
       
    plt.figure(fignum)
    plt.plot(Rs*R_earth, G_Pot, label=f"{Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot - G_Pot_Basic, label=f"{Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot_Basic, label=f"basic {Lat}-{Long}; {lmax}")
    plt.title("geopotential against the radius of the Earth, loop lmax")
    plt.xlabel("Distance from the center to the surface of the Earth (in %)")
    plt.ylabel("local value of the geopotential (m^2/s^2)")
    plt.legend(fontsize = 8)
    plt.show(block=False)


def TEST_plotGeoPot_radius():
    HC, HS = imp.Fetch_Coef() 
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    plt.figure()
    plt.clf()
    for i in range (2, 15):
        TEST_plot_radius(1, i*2, 50, 50, HC, HS, 10, HC_topo, HS_topo)
    

def TEST_lmax_loop_long_line():
    """ plots the geoid height at the equator, around the world """
    plt.figure()
    plt.clf()
    plt.grid(True)
    
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()    
    lmax_topo = 10

    lmaxs = np.arange(3, 25, 2)
    for lmax in lmaxs:
        Lat = pi/180 * 40
        R = gmath.Get_Ellipsoid_Radius(Lat)
        Longs = np.linspace(0, 2*pi, 91)
        
        Geo_H = np.zeros(len(Longs))
        
        for i in range(len(Longs)):
            Long = Longs[i]
#            Geo_H[i] = Get_acceleration   (R, Lat, Long,    lmax, HC, HS); title_spec="Acceleration"
#            Geo_H[i] = Get_Topo_Height   (R, Lat, Long,    lmax_topo, HC_topo, HS_topo); title_spec="Topography height"
#            Geo_H[i] = Get_Geo_Pot       (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="GeoPot"
            Geo_H[i] = Get_Geoid_Height  (R, Lat, Long,    lmax, HC, HS); title_spec="Geoid height"
#            Geo_H[i] = Get_Geoid_Height2 (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="Geoid height"
        
        Longs = (Longs-pi) * 180/pi
        plt.plot(Longs, Geo_H, label=f"lx={lmax}")
    
    plt.suptitle(f"{title_spec} at equator (m) vs Longitude; loop lmax")
    plt.legend()

    
def TEST_lmax_loop_lat_line():
    """ plots the geoid height at the equator, around the world """
    plt.figure()
    plt.clf()
    plt.grid(True)
    
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()    
    lmax_topo = 10

    lmaxs = np.arange(3, 25, 2)
    for lmax in lmaxs:
        Long = pi/180 * 80
        Lats = np.linspace(0, pi, 91)
        
        Geo_H = np.zeros(len(Lats))
        
        for i in range(len(Lats)):
            Lat = Lats[i]
            R = gmath.Get_Ellipsoid_Radius(Lat)
            
#            Geo_H[i] = Get_acceleration   (R, Lat, Long,    lmax, HC, HS); title_spec="Acceleration"
#            Geo_H[i] = Get_Topo_Height   (R, Lat, Long,    lmax_topo, HC_topo, HS_topo); title_spec="Topography height"
#            Geo_H[i] = Get_Geo_Pot       (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="GeoPot"
            Geo_H[i] = Get_Geoid_Height  (R, pi/2 - Lat, Long,    lmax, HC, HS); title_spec="Geoid height"
#            Geo_H[i] = Get_Geoid_Height2 (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="Geoid height"
        
        Lats = (pi/2-Lats)
        plt.plot(Lats*180/pi, Geo_H, label=f"lx={lmax}")
    
    plt.suptitle(f"{title_spec} at equator (m) vs Latitude; loop lmax")
    plt.legend()
    
 
def Get_isopot_average ():
    """ Returns the geopotential average at the surface of the ellipsoid """
#    HC, HS = imp.Fetch_Coef()
#    HC_topo, HS_topo = imp.Fetch_Topo_Coef() 
#    
#    lmax_av = 29; lmax_topo_av = 48; tens_av = 10
#    Grid, _, _ = Gen_Grid (tens_av, Get_Geo_Pot, [lmax_av, HC, HS, lmax_topo_av, HC_topo, HS_topo])
#    mm = np.amin(Grid)
#    MM = np.amax(Grid)
#    W_0 = np.average(Grid)   
#    print(f"lmax_av = {lmax_av}; lmax_topo_av={lmax_topo_av}; tens={tens_av}; \n\tmm={mm}; \n\tMM={MM}; \n\tW_0={W_0}")   
#  
#    mm=62499132.77072437
#    MM=62680791.364166744 
    W_0=62601662.83663869
    return W_0


def TEST_Get_isopot ():
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef() 
    Lat =  11
    Long =  10
    lmax = 5
    lmax_topo = 10
    
    W_0 = Get_isopot_average()
    R_e = gmath.Get_Ellipsoid_Radius(Lat)
    R_iso, height = Get_isopot(R_e, pi/180*Lat, pi/180*Long, W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)
    
    print(f"Average potential at {Lat} {Long} is at R={R_iso}, H={height}")


def TEST_Cosine_corr():
    HC, HS = imp.Fetch_Coef()
    lmax = 10
    CCos_2n0 = np.zeros(lmax)
    
    print("Ref. Gravity Potential Even Zonal Terms (C-form)")
    print("Subtracted From the Input Coefficients:")
    
    for l in range(1,lmax+1):
        CCos_2n0[l-1] = Cosine_Correction(2*l, HC[2*l, 0])
        print(f"C({2*l},0) =\t{CCos_2n0[l-1]}")
    
    
    
# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    
#    TEST_plotGeoPot_radius()
    TEST_lmax_loop_lat_line()
#    TEST_lmax_loop_long_line()
#    TEST_Get_isopot () 
#    TEST_Cosine_corr()
    
    print("\nGH_generate done")

