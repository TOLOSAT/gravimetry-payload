"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to display various graphs and maps about
    Topology coefficients

# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
# You might need to comment these two lines out
import os
os.environ['PROJ_LIB'] = r'C:\Users\Xavier\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from numpi import pi, sin, cos


import GH_import       as imp
#import GH_convert      as conv
import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
import GH_export       as exp
#import GH_displayTopo  as dtopo
import GH_terminal     as term
import GH_basemap      as bmp
#import GH_harmonics    as harm

# =============================================================================
# Generate functions
# =============================================================================
def Gen_Topo (lmax, HC, HS, tens, lmax_topo=10, HC_topo=[], HS_topo=[]):
    """
    This function generates an array containing Earth's topology at Lat/Long
    coordinates.
    Input:
        measure: "topo", "geopot" or "geoid", the measurement to be mapped on Earth
        lmax: degree to which the topology should be calculated
        HC_topo: Harmonic cosine coefficients to earth's topology
        HS_topo: Harmonic sine coefficients to earth's topology
        tens: how large the array should be
    Output:
        G_Height: array of grid height
        G_long: grid of longitudes
        G_lat: grid of latitudes

    """
    
    G_Grid, G_Long, G_Lat, Line_long, Line_lat, size_long, size_lat, points = bmp.init_grid(tens)   
    print(f"Generating Topology grid for lmax = {lmax}, {points} points")

    
    for i in range(0, size_long):
        term.printProgressBar(i+1, size_long)
        Long = Line_long[i]

        for j in range(0, size_lat):
            Lat = Line_lat[j]
            G_Grid[j,i] = harm.Get_Topo_Height(lmax, Lat, Long, HC, HS)
              
    return G_Grid, G_Long*180/pi, G_Lat*180/pi # in degrees now

















# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================
def Map_Topo (fignum, lmax, HC_topo, HS_topo, tens, levels, title):
    """
    makes a map of the topological coefficients
    """
    print("Plotting The Topology")

    # Get height grid
    G_Height, G_Long, G_Lat = gen.Gen_Grid ("topo", lmax, HC_topo, HS_topo, tens)

    print("Plotting Topology map")
    FIG = plt.figure(fignum)
    plt.clf()
    AX = FIG.add_subplot(111)

    """plot parameters"""
    alpha = 1
    map_colors = "terrain"
#    map_colors = "gist_earth"

    # Make map
    MAP = bmp.Gen_Basemap(FIG.number)
    MAP.drawcoastlines(linewidth = 0.4)

    # Display of Gm_Height, with coordinates G_phi and G_theta
    MAP.contourf(G_Long, G_Lat, G_Height, latlon = True,
                levels = levels, alpha = alpha,
                cmap=plt.get_cmap(map_colors))

    """plot apperance"""
#    FIG.set_size_inches(8, 6) #(36, 24)
    font_s = 10
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{1+36*tens}x{1+18*tens} points; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize = font_s)

    # add a colorbar
    CBAR = MAP.colorbar(location='bottom',pad="5%")
    CBAR.set_label("Height from sea level in meters")

    plt.axis('off')
    plt.show(block=False)
    return FIG


# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def TEST_Plots():
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax = 10
    tens = 2
    levels = 35
    title = f"TEST map of topology"
    fig1 = Map_Topo(1, lmax, HC_topo, HS_topo, tens, levels, title)
    exp.Store_Figure(fig1.number, "test")

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_Plots();
    print("\nGH_displayCoef done")

