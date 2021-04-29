"""
@authors:
# =============================================================================
 Information:
    The functions in this script all regard matters related to matplotlib
    figures and cartopy maps. This script replaces previous works done with the
    discontinued Basemap library.
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors

import numpy as np
from numpy import pi, sin, cos

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


# =============================================================================
# SETUP
# =============================================================================
colors.TwoSlopeNorm(vmin=-10000, vcenter=0., vmax=10000)


# =============================================================================
# FIGURE FUNCTIONS
# =============================================================================
def Make_Map_Fig (proj, fignum, ax_pos, shape, limits):
    """ Generates a mpl figure with the wanted coordiates system projection """
    FIG = plt.figure(*fignum, figsize=shape)
    AX = FIG.add_subplot(ax_pos, projection=proj(central_longitude=0))
    AX.set_extent(limits, crs=ccrs.PlateCarree())
    plt.show(block=False)
    return FIG, AX


def Make_Map (proj=ccrs.PlateCarree, fignum=[], ax_pos=111, shape=(7,5), limits=np.array([-180,180,-90,90])):
    """ Adds gridlines, credits and coastlines to a mpl figure """
    FIG, AX = Make_Map_Fig(proj, fignum, ax_pos, shape, limits)
    Add_Gridlines(AX, proj)
    Add_Credits(AX)
    AX.coastlines(linewidth = 0.6)
    return FIG, AX


def Make_Map_3D (fignum=[], ax_pos=111, shape=(7,5) ):
    """ Generates a 3d mpl figure """
    FIG = plt.figure(*fignum, figsize=shape)
    AX = FIG.add_subplot(ax_pos, projection='3d')
    return FIG, AX



# =============================================================================
# PLOT FUNCTIONS
# =============================================================================
def Plot_contourf(G_Grid, G_Long, G_Lat, AX=0, levels=75, proj=ccrs.PlateCarree, map_color="jet"):
    """
    Display of G_Grid, with coordinates G_Long and G_Lat
    map_colors = ["jet", "terrain", "gist_earth"]
    """
    if (AX==0): AX = plt.gca()
    alpha = 1

    data = AX.contourf(G_Long, G_Lat, G_Grid,
                       levels = levels, alpha = alpha,
                       transform = proj(), cmap=plt.get_cmap(map_color))
    CBAR = plt.colorbar(mappable=data, ax=AX, cmap=plt.get_cmap(map_color),
                        orientation='horizontal', pad=0.10)
    return CBAR


def Plot_surface (G_Grid, G_Long, G_Lat, AX=0, map_color="jet"):
    """
    3D Display of G_Grid surface, with coordinates G_Long and G_Lat
    map_colors = ["jet", "terrain", "gist_earth"]
    """
    if (AX==0): AX = plt.gca()
    alpha = 1

    data = AX.plot_surface(G_Long, G_Lat, G_Grid,
                           alpha = alpha, antialiased=False,
                           cmap=plt.get_cmap(map_color))
    AX.set_xlabel("Longitude",rotation=90)
    AX.set_ylabel("Latitude",rotation=90)
    CBAR = plt.colorbar(mappable=data, ax=AX, cmap=plt.get_cmap(map_color),
                        orientation='horizontal', pad=0.10)
    return CBAR


def Plot_surface_3D (G_Grid, G_Long, G_Lat, AX=0, ratio=0.15, map_color="jet"):
    """
    Ball representation of G_Grid + radius, with coordinates G_Long and G_Lat
    map_colors = ["jet", "terrain", "gist_earth"]
    ratio=0: sphere earth. ratio = 1: chaos earth
    """
    if (AX==0): AX = plt.gca()
    AX.figure.set_size_inches((6,6))

    ranges = abs(G_Grid.max() - G_Grid.min())
    R = (G_Grid - G_Grid.min())*ratio + ranges*(1-ratio)
    X, Y, Z = conv.sph2cart_Grid(R, G_Long, G_Lat) #G_Grid +
    norm = colors.Normalize()
    cmap = cm.get_cmap(map_color)
    m = cm.ScalarMappable(cmap=cmap)
    AX.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cmap(norm(G_Grid)) )
    m.set_array(G_Grid)

    Add_white_box(AX, X, Y, Z)
    AX.set_xticklabels('')
    AX.set_yticklabels('')
    AX.set_zticklabels('')
#    AX.axis('off')
    CBAR = plt.colorbar(mappable=m, ax=AX, cmap=cmap,
                        orientation='horizontal', pad=-0.10)
    return CBAR


def Rotating_map_gif():
    pass # https://tex.stackexchange.com/questions/268830/drawing-spherical-harmonic-density-plots-on-the-surface-of-a-sphere-in-tikz-pgfp



# =============================================================================
# FIGURE EDITING FUNCTIONS
# =============================================================================
def Add_Credits(AX, loc=4):
    """ Adds a textbox with Grav Harm 3 credits """
    TOOL = "Grav Harm 3"
    CREDIT = "TOLOSAT"
    TEXT_BOX = AnchoredText(f"{TOOL} by {CREDIT}", loc=loc, prop={'size': 8}, frameon=True)
    AX.add_artist(TEXT_BOX)


def Add_Gridlines(AX, proj=ccrs.PlateCarree):
    """ writes the gridlines and labens on the cartopy map """
#        http://balbuceosastropy.blogspot.com/2015/06/spherical-harmonics-in-python.html
    if (proj != ccrs.PlateCarree):
        AX.gridlines(color="gray", alpha=0.4)
    else:
        GL = AX.gridlines(crs=proj(), draw_labels=True, linewidth=1,
                          color='gray', alpha=0.4, linestyle='--')
        GL.xlabels_top = False
        GL.ylabels_left = False
#        GL.xlocator = mticker.FixedLocator([-5, -1, 0, 3])
        GL.xformatter = LONGITUDE_FORMATTER
        GL.yformatter = LATITUDE_FORMATTER
        GL.xlabel_style = {'size': 8}
        GL.ylabel_style = {'size': 8}
#        GL.xlabel_style = {'color': 'red', 'weight': 'bold'}


def Add_white_box(AX, X, Y, Z):
    """ function that might help get equal axes in a 3D plot (stackoverflow) """
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
       AX.plot([xb], [yb], [zb], 'w')



# =============================================================================
# DATA FUNCTIONS
# =============================================================================
def get_limits(G_Long, G_Lat):
    limits = [G_Long[0][0], G_Long[0][-1], G_Lat[0][0], G_Lat[-1][0]]
    return limits



# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def Map_Earth ():  #proj_crs=ccrs.Mollweide ):
    """
    Creates a Matplotlib figure with a map of the Earth, colored continents
    and oceans, showing parallels and meridians
    """
    FIG = plt.figure(figsize=(11,4))
    plt.clf()
    limits=np.array([0,180,0,90])
    FIG, AX1 = Make_Map_Fig(ccrs.Mollweide,   [FIG.number], 121, (11,4), limits)
    FIG, AX2 = Make_Map_Fig(ccrs.PlateCarree, [FIG.number], 122, (11,4), limits)

    plt.suptitle("Maps of the Earth - The endless possibilities of Cartopy")
    plt.show(block=False)

    # =========================================================================
    plt.axes(AX1)
    plt.title("Mollweide projection, stock_img", fontsize=10)
    AX1.set_global()
    AX1.gridlines()
    AX1.coastlines(linewidth = 1.5)
    # AX1.stock_img()

    # =========================================================================
    plt.axes(AX2)
    plt.title("PlateCarree projection, LAND & OCEAN(110m) COASTLINE(50m) features", fontsize=10)
    AX2.set_extent([-7, 4, 47, 54], crs=ccrs.PlateCarree())
    Add_Gridlines(AX2)
    Add_Credits(AX1)

    water_color = "lightcyan"
    land_color = "peachpuff"
    AX2.add_feature(cfeature.LAND, facecolor = land_color)
    AX2.add_feature(cfeature.OCEAN, facecolor = water_color)

    # high_res_coastline = cfeature.NaturalEarthFeature(
    #         category = "physical", name = "coastline", scale='50m')
    # AX2.add_feature(high_res_coastline, edgecolor = "black", facecolor = "None", linewidth = 0.6)



# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    Map_Earth()


    # limits= np.array([-180, 180, -90, 90])          # WORLD
    # limits= np.array([-7, 15, 40, 54])              # CHANNEL
    # limits= np.array([-25, 30, 15, 65])             # WC EUROPE AND NW AFRICA
    # limits= np.array([100, 170, -50, 10])           # AUSTRALIA
    # limits= np.array([-180, 180, -90, -40])         # ANTARCTICA
    # limits= np.array([115, 140, -15, 5])         # SULAWESI
    # Make_Map(limits=limits)

    # Make_Map_3D()

    print("\nGH_displayCoef done")