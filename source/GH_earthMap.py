"""

@authors:

# =============================================================================
 Information:

    The functions in this script all regard matters related to mpl Basemap

todo: replace basemap with cartopy
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
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
# DISPLAY FUNCTIONS
# =============================================================================
def Make_Map (fignum, G_Grid, G_Long, G_Lat, levels=35, map_colors="jet"):
    """ 
    Generates a matplotlib figure, adds the G_grid as a contourf 
    """
    FIG = plt.figure(fignum)
    plt.clf()
    AX = FIG.add_subplot(111)

    """plot parameters"""
    alpha = 1
#    map_colors = "jet"
#    map_colors = "terrain"
#    map_colors = "gist_earth"
    
    # Make map
    MAP = Gen_Basemap(FIG.number)
    MAP.drawcoastlines(linewidth = 0.4)

    # Display of Gm_Height, with coordinates G_phi and G_theta
    MAP.contourf(G_Long, G_Lat, G_Grid, latlon = True,
                levels = levels, alpha = alpha,
                cmap=plt.get_cmap(map_colors))

    # add a colorbar
    CBAR = MAP.colorbar(location='bottom',pad="5%")

    plt.axis('off')
    plt.show(block=False)
    
    return FIG, AX, MAP, CBAR




def Map_Earth ():  #proj_crs=ccrs.Mollweide ):
    """
    Creates a Matplotlib figure with a map of the Earth, colored continents 
    and oceans, showing parallels and meridians
    """
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    FIG = plt.figure(figsize=(11,4))
    plt.clf()
    
    # =========================================================================
    AX1 = FIG.add_subplot(121, projection = ccrs.Mollweide())
    AX1.set_global()    
    AX1.coastlines(linewidth = 1.5)
    AX1.stock_img()
    plt.title("Mollweide projection, stock_img", fontsize=10)

    # =========================================================================
    AX2 = FIG.add_subplot(122, projection = ccrs.PlateCarree())
    AX2.set_extent([-7, 4, 47, 54])
    AX2.gridlines()
    
    high_res_coastline = cfeature.NaturalEarthFeature(
            category = "physical", 
            name = "coastline",
            scale='50m')
    
    water_color = "lightcyan"
    land_color = "peachpuff"
    AX2.add_feature(cfeature.LAND, facecolor = land_color)
    AX2.add_feature(high_res_coastline, edgecolor = "black", linewidth = 0.6)
    AX2.add_feature(cfeature.OCEAN, facecolor = water_color)
    gl = AX2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-5, -1, 0, 3])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    plt.title("PlateCarree projection, LAND COASTLINE OCEAN features", fontsize=10)


    plt.suptitle("Maps of the Earth")
    plt.show(block=False)




# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_MAP():
    Map_Earth()

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_MAP()

    print("\nGH_displayCoef done")

