"""

@authors:

# =============================================================================
 Information:

    The functions in this script all regard matters related to matplotlib 
    figures and cartopy maps. This script replaces previous works done with the 
    discontinued Basemap library.

todo: replace basemap with cartopy
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
import numpy as np
from numpy import pi, sin, cos

#import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_basemap      as bmp
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


'''
# =============================================================================
# FUNCTIONS - BASEMAP PARAMETERS
# =============================================================================
def init_basemap (style = "crude mill"):

    if (style == "crude mill"):
        proj = "mill" # projection
        LatS = -90 # llcrnrlat
        LatN = 90 # urcrnrlat
        LongW = -180 # llcrnrlon
        LongE = 180 # urcrnrlon
        TS = 20 # lat_ts -- I don't know what this is but everyone online uses it so yeah
        Res = "c" # resolution, Crude, Low, [Intermediate, High, Full] > download extensions

    elif (style == "low mill"):
        proj = "mill" # projection
        LatS = -90 # llcrnrlat
        LatN = 90 # urcrnrlat
        LongW = -180 # llcrnrlon
        LongE = 180 # urcrnrlon
        TS = 20 # lat_ts -- I don't know what this is but everyone online uses it so yeah
        Res = "l" # resolution, Crude, Low, [Intermediate, High, Full] > download extensions

    else:
        proj = "mill" # projection
        LatS = -90 # llcrnrlat
        LatN = 90 # urcrnrlat
        LongW = -180 # llcrnrlon
        LongE = 180 # urcrnrlon
        TS = 20 # lat_ts -- I don't know what this is but everyone online uses it so yeah
        Res = "c" # resolution, Crude, Low, [Intermediate, High, Full] > download extensions

    # Bm_Param = [proj, LatS, LatN, LongW, LongE, TS, Res]
    return proj, LatS, LatN, LongW, LongE, TS, Res


# =============================================================================
# FUNCTIONS TO MAKE MAPS
# =============================================================================
def Gen_Basemap (fignum, style = "crude mill"):
    """
    Generates a Basemap map *projection* in the figure numbered fignum
    """
    plt.figure(fignum)

    proj, LatS, LatN, LongW, LongE, TS, Res = init_basemap(style)

    MAP = Basemap(projection = proj,
                llcrnrlat = LatS,
                urcrnrlat = LatN,
                llcrnrlon = LongW,
                urcrnrlon = LongE,
                lat_ts = TS,
                resolution = Res)
    return MAP
'''

# =============================================================================
# FIGURE EDITING FUNCTIONS
# =============================================================================
def Add_Credits(AX):
    TOOL = "Grav Harm 3"
    CREDIT = "TOLOSAT"
    TEXT_BOX = AnchoredText(f"{TOOL} by {CREDIT}",
                            loc=4, prop={'size': 8}, 
                            frameon=True)
    AX.add_artist(TEXT_BOX)


def Add_Gridlines(AX, cr_sys=ccrs.PlateCarree):
    GL = AX.gridlines(crs=cr_sys(), draw_labels=True, linewidth=1, 
                      color='gray', alpha=0.5, linestyle='--')
    GL.xlabels_top = False
    GL.ylabels_left = False
#    GL.xlocator = mticker.FixedLocator([-5, -1, 0, 3])
    GL.xformatter = LONGITUDE_FORMATTER
    GL.yformatter = LATITUDE_FORMATTER
    GL.xlabel_style = {'size': 8}
    GL.ylabel_style = {'size': 8}
#    GL.xlabel_style = {'color': 'red', 'weight': 'bold'}



# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================
def Plot_contourf(G_Grid, G_Long, G_Lat, AX, levels=35, proj=ccrs.PlateCarree, map_color="jet"):
    """
    Display of G_Grid, with coordinates G_Long and G_Lat 
    map_colors = ["jet", "terrain", "gist_earth"]
    """
    alpha = 1
    plt.axes(AX)
    data = AX.contourf(G_Long, G_Lat, G_Grid,
                       levels = levels, alpha = alpha,
                       transform = proj(), cmap=plt.get_cmap(map_color))
    CBAR = plt.colorbar(mappable=data, ax=AX, cmap=plt.get_cmap(map_color), 
                        orientation='horizontal', pad=0.10)    
    return CBAR


def Make_Map_Fig (proj=ccrs.PlateCarree, fignum=[], ax_pos=111, shape=(7,5) ):
    """ Generates a mpl figure with the wanted coordiates system projection """
    FIG = plt.figure(*fignum, figsize=shape)
    AX = FIG.add_subplot(ax_pos, projection=proj() )
    plt.show(block=False)
    return FIG, AX


def Make_Map (proj=ccrs.PlateCarree, fignum=[], ax_pos=111, shape=(7,5) ):
    """ Adds gridlines, credits and coastlines to a mpl figure """
    FIG, AX = Make_Map_Fig(proj, fignum, ax_pos, shape)
    Add_Gridlines(AX)    
    Add_Credits(AX)
    AX.coastlines(linewidth = 0.6)
    return FIG, AX

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
    FIG, AX1 = Make_Map_Fig(ccrs.Mollweide,   [FIG.number], 121, (11,4) )
    FIG, AX2 = Make_Map_Fig(ccrs.PlateCarree, [FIG.number], 122, (11,4) )
 
    plt.suptitle("Maps of the Earth")        
    plt.show(block=False)
    
    # =========================================================================
    plt.axes(AX1)
    plt.title("Mollweide projection, stock_img", fontsize=10)
    AX1.set_global() 
    AX1.gridlines()
    AX1.coastlines(linewidth = 1.5)
    AX1.stock_img()

    # =========================================================================
    plt.axes(AX2)
    plt.title("PlateCarree projection, LAND & OCEAN(110m) COASTLINE(50m) features", fontsize=10)
    AX2.set_extent([-7, 4, 47, 54])
    Add_Gridlines(AX2)    
    Add_Credits(AX1)
    
    water_color = "lightcyan"
    land_color = "peachpuff"
    AX2.add_feature(cfeature.LAND, facecolor = land_color)
    AX2.add_feature(cfeature.OCEAN, facecolor = water_color)
    
    high_res_coastline = cfeature.NaturalEarthFeature(
            category = "physical", name = "coastline", scale='50m')
    AX2.add_feature(high_res_coastline, edgecolor = "black", facecolor = "None", linewidth = 0.6)

    

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
#    Map_Earth()
#    Make_Map()
    print("\nGH_displayCoef done")

