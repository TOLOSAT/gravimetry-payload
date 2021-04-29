"""
@authors:
# =============================================================================
 Information:
    The purpose of this script is to display various graphs and maps about
    Geoid coefficients
    Generally used variables:
        lmax   = maximum degree to calculate to for the geopotantial
        HC, HS = Geopotential stokes coefficients
        lmax_topo, HC_topo, HS_topo = same, but for topography.
        limits = [Western_long, Eastern_long, Southern_lat, Northern_lat]
        mins = The grid resolution in arc minutes
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
import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
import GH_harmonics    as harm
#import GH_geoMath      as gmath
import GH_earthMap     as emap


'''
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
'''

# =============================================================================
# MAPPING FUNCTIONS
# =============================================================================

def Map_Geoid (mins, levels, title,    lmax, HC, HS, limits=np.array([-180,180,-90,90])):
    """ Makes a Matplotlib figure with the map, geoid and labels
    """
    # Get the data
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (mins, harm.Get_Geoid_Height,
                                           [lmax, HC, HS],
                                           limits)
    # save grid
    detail = f"grid geoid l{lmax}"
    exp.Store_temp_GLl(G_Grid, G_Long, G_Lat, detail)


    # Make a map
    FIG, AX = emap.Make_Map(limits=limits)#proj = ccrs.Mollweide)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels)
#    FIG, AX = emap.Make_Map_3D()
#    CBAR = emap.Plot_surface(G_Grid, G_Long, G_Lat, AX)
#    AX.set_zlabel("Geoid Height (m)",rotation=90)

    # Adapt labels
    plt.figure(FIG.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Geoid height in m")
    return FIG, [G_Grid, G_Long, G_Lat]


def Map_GeoPot (mins, levels, title,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo, limits=np.array([-180,180,-90,90])):
    """ Makes a Matplotlib figure with the map, geopotential and labels
    """
    # Get the data
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (mins, harm.Get_Geo_Pot,
                                           [lmax, HC, HS, lmax_topo, HC_topo, HS_topo],
                                           limits)
    # Make a map
    FIG, AX = emap.Make_Map(limits=limits)#proj = ccrs.Mollweide)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels)
#    FIG, AX = emap.Make_Map_3D()
#    CBAR = emap.Plot_surface_3D(G_Grid, G_Long, G_Lat, AX)
#    AX.set_zlabel("Geopotential (m)",rotation=90)

    # Adapt labels
    plt.figure(FIG.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; lmax_topo = {lmax_topo} degrees; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Gravitational potential in m^2/s^2")
    return FIG, [G_Grid, G_Long, G_Lat]


def Map_isoPot (mins, levels, title,     W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo, limits=np.array([-180,180,-90,90])):
    """ Makes a Matplotlib figure with the map, isopotential and labels
    """
    # Get the data
    G_Grid, G_Long, G_Lat = harm.Gen_Grid (mins, harm.Get_isopot,
                                           [W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo],
                                           limits)
    # Make a map
    FIG, AX = emap.Make_Map(limits=limits)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels)
    # Adapt labels
    plt.figure(FIG.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; lmax_topo = {lmax_topo} degrees; lmax = {lmax} degrees; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Height above reference ellipsoid where W(R)=W_0 (m)")
    return FIG, [G_Grid, G_Long, G_Lat]


def Map_Geoid_grid(detail="grid geoid l100", title="Geoid undulation High resolution"):
    """ Makes a Matplotlib figure with the map, geoid and labels
    """
    G_Grid, G_Long, G_Lat = imp.Load_GLl(detail)
    levels = 40
    limits = emap.get_limits(G_Long, G_Lat)

    FIG1, AX1 = emap.Make_Map(limits=limits) #, proj = ccrs.Mollweide)
    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX1, levels)
    plt.figure(FIG1.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; {detail}; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Geoid height in m")
    '''
    FIG2, AX2 = emap.Make_Map_3D()
    CBAR = emap.Plot_surface_3D(G_Grid, G_Long, G_Lat, AX2, ratio = 0.15)
    plt.figure(FIG2.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; {detail}; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Geoid height in m")
    FIG3, AX3 = emap.Make_Map_3D()
    CBAR = emap.Plot_surface(G_Grid, G_Long, G_Lat, AX3)
    plt.figure(FIG2.number)
    plt.suptitle(title)
    plot_specs = f"{G_Grid.size} points; {detail}; {levels} color levels"
    plt.title(plot_specs, fontsize=10)
    CBAR.set_label("Geoid height in m")
    '''
    return [G_Grid, G_Long, G_Lat]



# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_Map_Geoid():
    """  plots a quick geoid undulation map """
    HC, HS = imp.Fetch_Coef("full4")
    lmax = 10; mins = 600; levels = 70;
    title = f"Map of Geoid undulation"
    fig = Map_Geoid(mins, levels, title, lmax, HC, HS)
#    exp.Store_Figure(fig.number, f"test geoid", dpi=1000)

def TEST_Map_GeoPot():
    """ plots a quick geopotential map """
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax = 10; lmax_topo = 10; mins = 600; levels = 50
    title = f"TEST map of GeoPotential"
    _ = Map_GeoPot(mins, levels, title, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)


def TEST_Map_isoPot():
    """ plots a quick isopotential map """
    W_0 = harm.Get_isopot_average()
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    lmax = 5; lmax_topo = 5; mins = 600; levels = 50
    title = f"TEST map of isopotential W_0={W_0:.2f} m^2/s^2"
    _ = Map_isoPot(mins, levels, title, W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)



# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
#    TEST_Map_GeoPot()

#    TEST_Map_isoPot()

#    TEST_Map_Geoid()


    gx, LLx, llx = Map_Geoid_grid("Sulawesi", "Geoid undulation (Sulawesi)")


#    gx, LLx, llx = Map_Geoid_grid("grid geoid l100", "GH3 geoid")
#    gf, LLf, llf = Map_Geoid_grid("EGM2008 1h",         "F77 geoid")
#
#    '''
#    FIG, AX= emap.Make_Map()
#    emap.Plot_contourf(gx - gf, LLx, llx, AX, map_color="seismic")
#    plt.title("Difference between GH3 and F77 geoid interpretations")
#    '''


    print("\nGH_displayGeoid done")