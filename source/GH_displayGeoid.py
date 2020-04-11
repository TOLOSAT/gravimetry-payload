"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to display various graphs and maps about
    Geoid coefficients

# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
# You might need to comment these two lines out
import os
os.environ['PROJ_LIB'] = r'C:\Users\Xavier\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

import GH_import       as imp
#import GH_convert      as conv
import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
import GH_terminal     as term
import GH_basemap      as bmp
#import GH_harmonics    as harm


"""
# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================
def Plot_Array_Diff(HS_nm_slv, HC_nm_slv, fig_num = 6):
    print("plotting coeff difference")


    #resize the official coef
    HC, HS = imp.Fetch_Coef()
    HS_nm_sz = HS[:len(HS_nm_slv), :len(HS_nm_slv)]
    HC_nm_sz = HC[:len(HC_nm_slv), :len(HC_nm_slv)]

    #subtract calculated coeffs
    HS_nm_sz -= HS_nm_slv
    HC_nm_sz -= HC_nm_slv

    fig_HC = plt.figure(fig_num)
    plt.clf()
    plt.suptitle("Harmonic coeff difference between official and solved; degree: "+str(len(HS_nm_sz)-1))

    for n in range (0, len(HC_nm_sz)):
        Ms_n = np.arange(0, n+1)

        HC_ni = HC_nm_sz[n, :n+1]
        HS_ni = HS_nm_sz[n, :n+1]

        plt.subplot(211)
        plt.plot(Ms_n, HC_ni,'-*', label='n='+str(n))

        plt.subplot(212)
        plt.plot(Ms_n, HS_ni,'-*', label='n='+str(n))

    plt.subplot(211)
    plt.ylabel("COSINE coeff diff")
    plt.grid(True)
#    plt.xlabel("order m of derivation (log)")
#    plt.ylabel("value of HC_nm")
    plt.legend(loc = 'upper right', title = 'Degree n', fontsize = 5)

    plt.subplot(212)
    plt.ylabel("SINE coeff diff")
    plt.grid(True)
    plt.xlabel("order m of derivation (log)")
#    plt.ylabel("value of HS_nm")
#    plt.legend(loc = 'lower right', title = 'Degree n', fontsize = 5)

    plt.show()
"""


def Map_Geoid (fignum, lmax, HC, HS, tens, levels, title, lmax_topo):
    """
    Makes a map of given geoid coefficients
    """
    # Get geoid grid
    G_Geoid, G_Long, G_Lat = gen.Gen_Grid ("delta g", lmax, HC, HS, tens, lmax_topo)

    print("Plotting Geoid map\n")
    FIG = plt.figure(fignum)
    plt.clf()
    AX = FIG.add_subplot(111)

    """plot parameters"""
    alpha = 1
    map_colors = "jet"

    # Make map
    MAP = bmp.Gen_Basemap(FIG.number)
    MAP.drawcoastlines(linewidth = 0.4)

    # Display of Gm_Height, with coordinates G_phi and G_theta
    MAP.contourf(G_Long, G_Lat, G_Geoid, latlon = True,
                levels = levels, alpha = alpha,
                cmap=plt.get_cmap(map_colors))

    """plot apperance"""
    plt.suptitle(title)
    plot_specs = f"{1+36*tens}x{1+18*tens} points; lmax_topo = {lmax_topo} degrees; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)

    # add a colorbar
    CBAR = MAP.colorbar(location='bottom',pad="5%")
    CBAR.set_label("Gravitational potential in m^2/s^2")

    plt.axis('off')
    plt.show(block=False)
    return FIG

def Plot_GeoPot_height(fignum, lmax, Lat, Long, HC, HS):
    """
    Plots the geopotential at given coordinates from Earth's center to the 
    surface, and beyond into space.
    """
    
    R_earth = imp.Get_Radius(Lat)
    Rs = np.linspace(95, 105, 100)
    G_Pot = np.zeros(Rs.shape)
    G_Pot_Basic = np.zeros(Rs.shape)
    
    for i in range (len(Rs)):
        term.printProgressBar(i+1, len(Rs))
        G_Pot[i] = gen.Get_Geo_Pot (lmax, Rs[i]*R_earth, Lat, Long, HC, HS)
        G_Pot_Basic[i] = Math_calc_geopot_basic(Rs[i]*R_earth)
        
        
    plt.figure(fignum)
    plt.plot(Rs, G_Pot,       label=f"      {Lat}-{Long}; {lmax}")
    plt.plot(Rs, G_Pot_Basic, label=f"basic {Lat}-{Long}; {lmax}")
    plt.title("geopotential against the radius of the Earth")
    plt.xlabel("Distance from the center to the surface of the Earth (in %)")
    plt.ylabel("local value of the geopotential (m^2/s^2)")
    plt.legend(fontsize = 8)
    plt.show(block=False)
    
    
    
    
# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_plotGeoid():
    HC, HS = imp.Fetch_Coef()
    lmax = 10
    lmax_topo = 1
    tens = 3
    levels = 35
    title = f"TEST map of geoid"
    fig1 = Map_Geoid(2, lmax, HC, HS, tens, levels, title, lmax_topo)
    


def TEST_plotGeoPot_Height():
    HC, HS = imp.Fetch_Coef()    
    for i in range (1, 4):
        Plot_GeoPot_height(1, i*5, 0, 0, HC, HS)
    
    
def Math_calc_geopot_basic(z):
    G = 6.673E-11
    M = 5.975E24
    a = 6.378E6
    
    P = G*M*(1/a + 1/(a+float(z)))
    return P
    
# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_plotGeoid()
#    TEST_plotGeoPot_Height()
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
#        G_Pot[i] = gen.Get_Geo_Pot (lmax, Rs[i]*R_earth, Lat, Long, HC, HS)
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
        G_Pot[i] = gen.Get_Geo_Pot (lmax[i], R_earth, Lat, Long, HC, HS)
    plt.plot(lmax, G_Pot,       label=f"{Lat}-{Long}")
    plt.title("geopotential against lmax")
    plt.xlabel("lmax")
    plt.ylabel("local value of the geopotential (m^2/s^2)")   
    
    
    
    
    
    
    
    
    
    plt.legend(fontsize = 8)
    plt.show(block=False)
    
    print("\nGH_displayGeoid done")
"""
