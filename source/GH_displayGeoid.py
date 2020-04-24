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
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import numpy as np

import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
import GH_harmonics    as harm
#import GH_geoMath      as gmath
import GH_earthMap     as emap


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

# =============================================================================
# MAPPING FUNCTIONS
# =============================================================================

def Map_Geoid (tens, levels, title,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """ Makes a Matplotlib figure with the map, geoid and labels """
    # Get the data
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (tens, harm.Get_Geoid_Height, [lmax, HC, HS])
    # Make a map    
    FIG, AX = emap.Make_Map()#proj = ccrs.Mollweide)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels)     
    # Adapt labels
    plt.figure(FIG.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; lmax_topo = {lmax_topo} degrees; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Geoid height in m")    
    return FIG


def Map_GeoPot (tens, levels, title,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """ Makes a Matplotlib figure with the map, geopotential and labels """
    # Get the data
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (tens, harm.Get_Geo_Pot, [lmax, HC, HS, lmax_topo, HC_topo, HS_topo])
    # Make a map
    FIG, AX = emap.Make_Map()#proj = ccrs.Mollweide)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels)    
    # Adapt labels        
    plt.figure(FIG.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; lmax_topo = {lmax_topo} degrees; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Gravitational potential in m^2/s^2")    
    return FIG


def Map_isoPot (tens, levels, title,     W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """ Makes a Matplotlib figure with the map, isopotential and labels """
    # Get the data
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (tens, harm.Get_isopot, [W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo])
    # Make a map    
    FIG, AX = emap.Make_Map()
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels)   
    # Adapt labels
    plt.figure(FIG.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; lmax_topo = {lmax_topo} degrees; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Height above reference ellipsoid where W(R)=W_0 (m)")    
    return FIG



# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def TEST_Map_Geoid():
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax = 5; lmax_topo = 5; tens = 1; levels = 50; 
    title = f"TEST map of Geoid"
    _ = Map_Geoid(tens, levels, title, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)


def TEST_Map_GeoPot():
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax = 5; lmax_topo = 5; tens = 1; levels = 50 
    title = f"TEST map of GeoPotential"
    _ = Map_GeoPot(tens, levels, title, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)


def TEST_Map_isoPot():
    W_0 = harm.Get_isopot_average() 
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax = 5; lmax_topo = 5; tens = 1; levels = 50
    title = f"TEST map of isopotential W_0={W_0:.2f} m^2/s^2"
    _ = Map_isoPot(tens, levels, title, W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)
   
    
    
# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
#    TEST_Map_GeoPot()
    
#    TEST_Map_isoPot()
    
    TEST_Map_Geoid()
    
    print("\nGH_displayGeoid done")


