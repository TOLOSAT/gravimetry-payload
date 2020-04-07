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
import os   
os.environ['PROJ_LIB'] = r'C:\Users\Xavier\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

from mpl_toolkits.basemap import Basemap
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
#import GH_terminal     as term
import GH_basemap      as bmp


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
#    font_s = 10 # 50  
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{1+36*tens}x{1+18*tens} points; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs) #, fontsize = font_s)

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
    lmax = 8
    tens = 1
    levels = 35
    title = f"TEST map of topology"
    fig1 = Map_Topo(1, lmax, HC_topo, HS_topo, tens, levels, title) 
    
    
# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    TEST_Plots();
    print("\nGH_displayCoef done")

