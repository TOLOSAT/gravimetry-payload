"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to display various graphs and maps about
    Topography coefficients

# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
# You might need to comment these two lines out
#import os
#os.environ['PROJ_LIB'] = r'C:\Users\Xavier\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
from numpy import pi, sin, cos

import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
import GH_basemap      as bmp
import GH_harmonics    as harm


# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================
def Map_Topo (fignum, lmax, HC_topo, HS_topo, tens, levels, title):    
    """
    Makes a Matplotlib figure with the map and labels for the topography
    """
    G_Grid, G_Long, G_Lat = harm.Gen_Topo (lmax, HC_topo, HS_topo, tens)
    map_colors = "terrain"
#    map_colors = "gist_earth"
    FIG, AX, MAP, CBAR = bmp.Make_Map (fignum, G_Grid, G_Long, G_Lat, levels, map_colors)
    
    # Adapt labels
    plt.figure(FIG.number)
    font_s = 10
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{G_Grid.size} points; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize = font_s)
    CBAR.set_label("Height from sea level in meters") # geopot
    return FIG


# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def TEST_Map_topo():
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax_topo = 10
    tens = 1
    levels = 50
    title = f"TEST map of topography"
    fig1 = Map_Topo(1, lmax_topo, HC_topo, HS_topo, tens, levels, title)
#    exp.Store_Figure(fig1.number, "test")

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_Map_topo() 
    print("\nGH_displayCoef done")

