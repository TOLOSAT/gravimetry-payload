"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to calculate the sums from spherical 
    harmonic coefficients

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
#import GH_basemap      as bmp
#import GH_harmonics    as harm


# =============================================================================
# FUNCTIONS TO CALCULATE SPHERICAL HARMONIC SUMS
# =============================================================================
def Get_Topo_Height (lmax_topo, Lat, Long, HC_topo, HS_topo):
    """
    This function returns the height of Earth's estimated topography at Lat/Long 
    coordinates
    The solution is calculated up to degree lmax in the [HC, HS] model
    """
    Sum1 = 0
    P_lm, _ = imp.Pol_Legendre(lmax_topo, lmax_topo, cos(Lat)) # I am allowed to write that.
    for l in range (0, lmax_topo+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC_topo[l,m]*cos(m*Long) + HS_topo[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum1 += Sum2

    return Sum1



def Get_Geo_Pot (lmax, R, Lat, Long, HC, HS):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    """
    cst = imp.Constants()

    Sum1 = 0
    P_lm, _ = imp.Pol_Legendre(lmax, lmax, cos(Lat))
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum1 += (cst.a_g/R)**l * Sum2

    geopot = cst.GM_g/R*(1 + Sum1)

    return geopot



def Get_Geoid_Height (lmax, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    
    Equations from geoid cook book
    """
    cst = imp.Constants()    
    R = imp.Get_Ellipsoid_Radius(Lat)
    g_0 = imp.Get_Normal_Gravity(Lat)
    Lat_gc = conv.geodes2geocen(Lat)
    
    R = cst.a_g 
    g_0 = cst.g 
    Lat_gc = Lat
    
    
    Sum1 = 0
    P_lm, _ = imp.Pol_Legendre(lmax, lmax, sin(Lat_gc) )
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum1 +=  (cst.a_g/R)**l * Sum2 

    Geo_H = cst.GM_g * Sum1 / (R*g_0)

    return Geo_H


def Get_Geoid_Height2 (lmax, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    
    Equations from GHZ handbook
    """
    c = imp.Constants() 
    a_g=c.a_g; GM_g=c.GM_g; # g_e=c.g 
    G=c.G; ro=c.ro;  
    
    R_e = imp.Get_Ellipsoid_Radius(Lat)
    g_e = imp.Get_Normal_Gravity(Lat)
#    Lat_gc = conv.geodes2geocen(Lat)
    Lat_gc = Lat
    
    
    Sum_geo = 0
    P_lm, _ = imp.Pol_Legendre(lmax, lmax, sin(Lat_gc) )
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long)+HS[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum_geo +=  (a_g/R_e)**l * Sum2 


    Sum_topo = 0
    P_lm, _ = imp.Pol_Legendre(lmax_topo, lmax_topo, cos(Lat))
    for l in range (0, lmax_topo+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC_topo[l,m]*cos(m*Long) + HS_topo[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum_topo += Sum2

    Geo_H = GM_g/(R_e*g_e) * Sum_geo  -  2*pi*G*ro/g_e * (R_e*Sum_topo)**2


    return Geo_H


def Get_delta_g (lmax, R, Lat, Long, HC, HS):
    """
    This function returns delta_g, whatever that is, from the GFZ handbook
    """
    GM, a, g = imp.Get_Grav_constants()
        
    Sum1 = 0
    P_lm, _ = imp.Pol_Legendre(lmax, lmax, cos(Lat))
    for l in range (0, lmax+1):
        Sum2 = 0
        for m in range(0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum1 +=  (l-1) * (a/R)**l * Sum2
    D_g = GM * Sum1 / R**2
    return D_g


def Get_acceleration (lmax, R, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the acceleration at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The acceleration is the gratient of the geopotential, and is calculated 
    over a "d" distance 
    """
    c = imp.Constants()
    a_g=c.a_g; GM_g=c.GM_g; 
    
    d = 1 # m
    R_1 = R - d/2
    R_2 = R + d/2
    
    Sum1 = 0
    P_lm, _ = imp.Pol_Legendre(lmax, lmax, cos(Lat))
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * P_lm[m, l] * imp.Normalize(l, m)
        Sum1 += (a_g/R)**l * Sum2

    PG_1 = GM_g/R_1*(1 + Sum1)
    PG_2 = GM_g/R_2*(1 + Sum1)
    
    local_acc = (PG_1 - PG_2) / d
    
    return local_acc


# =============================================================================
# DICOTOMY POINT BY POINT
# =============================================================================
def dichotomy_grad (f, in_first, z_e, in_after, w_0, de, grad):
    """
    f : function for f(*in_first,z,*in_after) = w
    in_first, in_after : f function variables that come before and after z_e
    z_e : firstguess for the value to be tested
    w_0: target value for w
    de : delta error  
    grad : gradient
    """
    w_i = f(*in_first, z_e, *in_after)
    di = w_0 - w_i
    z_i = z_e    
    c = 0
    print(f"w_0={w_0:.2f}; z_i={z_i:.2f}; w_i={w_i:.2f}; di={di:.2f}; add={(di/grad):.2f}; "); 
    
    while (abs(di) >= de):
        c+=1        
        z_i += di/grad
        w_i = f(*in_first, z_i, *in_after)
        di = w_0 - w_i        
#        sleep(1); 
        print(f"w_0={w_0:.2f}; z_i={z_i:.2f}; w_i={w_i:.2f}; di={di:.2f}; add={(di/grad):.2f}; "); 
    print(f"dichotomy_grad: {c} steps")
    return z_i


def Get_isopot (lmax, R_, W, Lat, Long, HC, HS):
    """
    This function returns the Height at given Lat/Long coordinates at which the 
    geopotential is equal to the given W
    The solution is calculated up to degree lmax in the HC HS model
    The approach is a dichotomy within a width of "e" m error
    """
#    c = imp.Constants()
#    a_g=c.a_g; GM_g=c.GM_g; W = c.W_0
    
    de = 10
    grad = -10 # 9.81
    R_e = imp.Get_Ellipsoid_Radius(Lat) #+ Get_Topo_Height(lmax_topo, Lat, Long, HC_topo, HS_topo) # add Earth's radius !!
    
    R_iso = dichotomy_grad (Get_Geo_Pot, [lmax], R_e, [Lat, Long, HC, HS], W, de, grad)
    Height = R_iso - R_e
    return R_iso, Height




# =============================================================================
# FUNCTIONS TO GENERATE DATA ARRAYs
# =============================================================================
def init_grid (tens):
    """
    Initiates the grid variables based on the number of points wanted
    """
    size_long = 1 + 36*tens
    size_lat  = 1 + 18*tens

    Line_long = np.linspace(0, 2*pi, size_long) # 0 to 360 ; must subtract 180
    Line_lat  = np.linspace(0, pi, size_lat) # 0 to 180 ; must do 90 - theta
    G_Long, G_Lat = np.meshgrid((Line_long - pi), (pi/2 - Line_lat))

    G_Grid = np.zeros((size_lat, size_long))
    
    return G_Grid, G_Long, G_Lat


def Gen_Grid (measure, lmax, HC, HS, tens, lmax_topo=0, HC_topo=[], HS_topo=[]):
    """
    This function generates a grid of the desired spherical harmonic model
    at Lat/Long coordinates
    Input:
        measure: "geopot" or "geoid", the measurement to be mapped on Earth
        lmax: degree to which the topography should be calculated
        HC: Harmonic cosine coefficients 
        HS: Harmonic sine coefficients 
        tens: how large the array should be
    Output:
        G_Height: array of grid height
        G_long: grid of longitudes
        G_lat: grid of latitudes

    """
    GM, a, g = imp.Get_Grav_constants()
    G_Grid, G_Long, G_Lat = init_grid(tens)   
    print(f"Generating {measure} grid for lmax = {lmax}, {G_Grid.size} points")

#    if (measure == "geopot"):
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    
    for j in range(0, G_Lat.shape[0]):
        term.printProgressBar(j+1, G_Lat.shape[0])
        Lat = pi/2 - G_Lat[j][0]
    
        for i in range(0, G_Long.shape[1]):
            Long = G_Long[0][i] - pi

            if (measure == "geopot"):
                R = imp.Get_Ellipsoid_Radius(Lat) + Get_Topo_Height(lmax_topo, Lat, Long, HC_topo, HS_topo) # add Earth's radius !!
                G_Grid[j,i] = Get_Geo_Pot(lmax, R, Lat, Long, HC, HS)
                
            elif (measure == "geoid"):
#                a = imp.Get_Radius(Lat)
                G_Grid[j,i] = Get_Geoid_Height(lmax, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo)
                
            elif (measure == "geoid2"):
                G_Grid[j,i] = Get_Geoid_Height2(lmax, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo)
                
            elif (measure == "acceleration"):
                R = imp.Get_Ellipsoid_Radius(Lat)
                G_Grid[j,i] = Get_acceleration(lmax, R, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo)
                
            elif (measure == "delta g"):
                a = 6378136.3 # m
                a = imp.Get_Radius(Lat)
                G_Grid[j,i] = Get_delta_g(lmax, a, Lat, Long, HC, HS)   

                
                
    return G_Grid, G_Long*180/pi, G_Lat*180/pi # in degrees now


def Gen_Topo (lmax_topo, HC_topo, HS_topo, tens):
    """
    This function generates an array containing Earth's topography at Lat/Long
    coordinates 
    Input:
        lmax: degree to which the topography should be calculated
        HC_topo: Harmonic cosine coefficients to earth's topography
        HS_topo: Harmonic sine coefficients to earth's topography
        tens: token of the volume of data desired
    Output:
        G_Height: array of grid height
        G_long: grid of longitudes
        G_lat: grid of latitudes

    """    
    G_Grid, G_Long, G_Lat = init_grid(tens)   
    print(f"Generating Topography grid for lmax = {lmax_topo}, {G_Grid.size} points")

    for j in range(0, G_Lat.shape[0]):
        term.printProgressBar(j+1, G_Lat.shape[0])
        Lat = pi/2 - G_Lat[j][0]
    
        for i in range(0, G_Long.shape[1]):
            Long = G_Long[0][i] - pi
            
            G_Grid[j,i] = Get_Topo_Height(lmax_topo, Lat, Long, HC_topo, HS_topo)
            
    return G_Grid, G_Long*180/pi, G_Lat*180/pi # in degrees now


def Gen_isopot(lmax, tens, HC, HS, lmax_av=5, tens_av=1):
    """ Generates the isopotential surface arround W_0 """
    
    Grid, _, _ = Gen_Grid ("geopot", lmax_av, HC, HS, tens_av)
    mm = np.amin(Grid)
    MM = np.amax(Grid)
    W_0 = np.average(Grid) # ------------------ the important one
    print(f"lmax_av = {lmax_av}; \n\tmin={mm}; \n\tMax={MM}; \n\tAvg={W_0}")
#    mm=62551117.58580653; MM=62660433.92621861; Av=62609103.69878714
    
    
    G_Grid, G_Long, G_Lat = init_grid(tens)   
    print(f"Generating isopotential grid for lmax = {lmax}, {G_Grid.size} points")
    
    for j in range(0, G_Lat.shape[0]):
        term.printProgressBar(j+1, G_Lat.shape[0])
        Lat = pi/2 - G_Lat[j][0]
        R_e = imp.Get_Ellipsoid_Radius(Lat)
    
        for i in range(0, G_Long.shape[1]):
            Long = G_Long[0][i] - pi
            
            R_iso, Height = Get_isopot(lmax, R_e, W_0, pi/180*Lat, pi/180*Long, HC, HS)
            G_Grid[j,i] = Height
                
    return G_Grid, G_Long*180/pi, G_Lat*180/pi # in degrees now


    
    
# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def Plot_GeoPot_height(fignum, lmax, Lat, Long, HC, HS):
    """
    Plots the geopotential at given coordinates from Earth's center to the 
    surface, and beyond into space.
    """
    
    R_earth = imp.Get_Ellipsoid_Radius(Lat)
    Rs = np.linspace(95, 105, 100)
    G_Pot = np.zeros(Rs.shape)
    G_Pot_Basic = np.zeros(Rs.shape)
    
    for i in range (len(Rs)):
        term.printProgressBar(i+1, len(Rs))
        G_Pot[i] = Get_Geo_Pot (lmax, Rs[i]*R_earth, Lat*pi/180, Long*pi/180, HC, HS)
        G_Pot_Basic[i] = Get_Geo_Pot (2, Rs[i]*R_earth, Lat*pi/180, Long*pi/180, HC, HS)
#        G_Pot_Basic[i] = Math_calc_geopot_basic(Rs[i]*R_earth)
        
        
    plt.figure(fignum)
    plt.plot(Rs*R_earth, G_Pot, label=f"{Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot - G_Pot_Basic, label=f"{Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot_Basic, label=f"basic {Lat}-{Long}; {lmax}")
    plt.title("geopotential against the radius of the Earth")
    plt.xlabel("Distance from the center to the surface of the Earth (in %)")
    plt.ylabel("local value of the geopotential (m^2/s^2)")
    plt.legend(fontsize = 8)
    plt.show(block=False)
    
    
 
def TEST_plotGeoPot_Height():
    HC, HS = imp.Fetch_Coef() 
    plt.figure(1)
    plt.clf()
    for i in range (2, 15):
        Plot_GeoPot_height(1, i*2, 50, 50, HC, HS)
    
    
def Math_calc_geopot_basic(z):
    G = 6.673E-11
    M = 5.975E24
    a = 6.378E6
    
    P = G*M*(1/a + 1/(a+float(z)))
    return P

def TEST_Geoid_Line():
    """ plots the geoid height at the equator, around the world """
    plt.figure(1)
    plt.clf()
    plt.grid(True)
    plt.suptitle("Geoid height at equator (m) vs Longitude")
    
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    
    lmax_topo = 10

    lmaxs = np.arange(3, 25, 2)

    for lmax in lmaxs:
        Lat = pi/180 * 40
        R = imp.Get_Ellipsoid_Radius(Lat)
        Longs = np.linspace(0, 2*pi, 91)
        
        Geo_H = np.zeros(len(Longs))
        
        for i in range(len(Longs)):
            Long = Longs[i]
            Geo_H[i] = Get_acceleration(lmax, R, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo)
        
        Longs = (Longs-pi) * 180/pi
        plt.plot(Longs, Geo_H, label=f"{lmax}")
    
    plt.legend()
    
    
def TEST_Get_isopot ():
    HC, HS = imp.Fetch_Coef()
    
    lmax_av = 5; 
    tens_av = 1;    
    Grid, _, _ = Gen_Grid ("geopot", lmax_av, HC, HS, tens_av)
    mm = np.amin(Grid)
    MM = np.amax(Grid)
    Av = np.average(Grid) # ------------------ the important one
    print(f"lmax_av = {lmax_av}; \n\tmin={mm}; \n\tMax={MM}; \n\tAvg={Av}")
#    mm=62551117.58580653; MM=62660433.92621861; Av=62609103.69878714
    
    Lat =  11
    Long =  10
    lmax = 5
    
    R_e = imp.Get_Ellipsoid_Radius(Lat)
    R_iso, height = Get_isopot(lmax, R_e, Av, pi/180*Lat, pi/180*Long, HC, HS)
    
    print(f"Average potential at {Lat} {Long} is at R={R_iso}, H={height}")
    
    
    
    
    

    
# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    
#    TEST_plotGeoPot_Height()
#    TEST_Geoid_Line()
    TEST_Get_isopot ()    
    
    
    
    """
    HC, HS = imp.Fetch_Coef() 
    plt.figure(1)
    
    Lat = 0
    Long = 0
    R_earth = imp.Get_Radius(Lat)

    
    
#    Rs = np.linspace(95, 105, 100)
#    G_Pot = np.zeros(Rs.shape)
#    G_Pot_Basic = np.zeros(Rs.shape)    
#    for i in range (len(Rs)):
#        term.printProgressBar(i+1, len(Rs))
#        G_Pot[i] = harm.Get_Geo_Pot (lmax, Rs[i]*R_earth, Lat, Long, HC, HS)
#        G_Pot_Basic[i] = Math_calc_geopot_basic(Rs[i]*R_earth)
#    plt.plot(Rs, G_Pot,       label=f"      {Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot_Basic, label=f"basic {Lat}-{Long}; {lmax}")
#    plt.title("geopotential against the radius of the Earth")
#    plt.xlabel("Distance from the center to the surface of the Earth (in %)")
#    plt.ylabel("local value of the geopotential (m^2/s^2)")
    
    
    lmax = np.arange(0, 20, 1)
    G_Pot = np.zeros(lmax.shape)   
    for i in range (len(lmax)):
        term.printProgressBar(i+1, len(lmax))
        G_Pot[i] = harm.Get_Geo_Pot (lmax[i], R_earth, Lat, Long, HC, HS)
    plt.plot(lmax, G_Pot,       label=f"{Lat}-{Long}")
    plt.title("geopotential against lmax")
    plt.xlabel("lmax")
    plt.ylabel("local value of the geopotential (m^2/s^2)")   
    
    
    plt.legend(fontsize = 8)
    plt.show(block=False)
    """
    
    
    print("\nGH_generate done")

