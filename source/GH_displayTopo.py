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
import cartopy.crs as ccrs
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
import GH_harmonics    as harm
#import GH_geoMath      as gmath
import GH_earthMap     as emap


# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================
def Map_Topo (lmax_topo, HC_topo, HS_topo, tens, levels, title, style="map"):    
    """ 
    Makes a Matplotlib figure with the map, topography and labels 
    """
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (tens, harm.Get_Topo_Height, [lmax_topo, HC_topo, HS_topo])
    map_color = "terrain"
#    map_colors = "gist_earth"
    
    if (style == "ball"):
        FIG, AX = emap.Make_Map_3D()
        CBAR = emap.Plot_surface_3D(G_Grid, G_Long, G_Lat, AX, map_color=map_color) 
    elif (style == "relief"):
        FIG, AX = emap.Make_Map_3D()
        CBAR = emap.Plot_surface(G_Grid, G_Long, G_Lat, AX, map_color=map_color)  
        AX.set_zlabel("Height (m)",rotation=90) 
    else:        
        FIG, AX = emap.Make_Map()# proj = ccrs.Mollweide)
        CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels, map_color=map_color) 
    
    
    # Adapt labels
    font_s = 10
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{G_Grid.size} points; lmax = {lmax_topo} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize = font_s)
    CBAR.set_label("Height from sea level in meters") 
    return FIG


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def Demo_Map_Topo(lmax_topo, HC_topo, HS_topo, tens, levels, title):    
    """ 
    Makes a Matplotlib figure with the map, topography and labels 
    """
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (tens, harm.Get_Topo_Height, [lmax_topo, HC_topo, HS_topo])
    map_color = "terrain"
#    map_colors = "gist_earth"

    FIG, AX = emap.Make_Map_3D()
    CBAR = emap.Plot_surface_3D(G_Grid, G_Long, G_Lat, AX, map_color=map_color) 
    font_s = 10
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{G_Grid.size} points; lmax = {lmax_topo} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize = font_s)
    CBAR.set_label("Height from sea level in meters")

    FIG, AX = emap.Make_Map_3D()
    CBAR = emap.Plot_surface(G_Grid, G_Long, G_Lat, AX, map_color=map_color)  
    AX.set_zlabel("Height (m)",rotation=90) 
    font_s = 10
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{G_Grid.size} points; lmax = {lmax_topo} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize = font_s)
    CBAR.set_label("Height from sea level in meters")
    
    FIG, AX = emap.Make_Map()# proj = ccrs.Mollweide)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels, map_color=map_color) 
    font_s = 10
    plt.suptitle(title) #, fontsize = font_s)
    plot_specs = f"{G_Grid.size} points; lmax = {lmax_topo} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize = font_s)
    CBAR.set_label("Height from sea level in meters") 
    
    return FIG

def TEST_Map_Topo():
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax_topo = 10
    tens = 1
    levels = 50
    title = f"TEST map of topography"    
    fig = Demo_Map_Topo(lmax_topo, HC_topo, HS_topo, tens, levels, title)
#    exp.Store_Figure(fig.number, "test")

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_Map_Topo() 
    print("\nGH_displayCoef done")

