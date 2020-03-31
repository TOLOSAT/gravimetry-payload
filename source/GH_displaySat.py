"""

@authors:

# =============================================================================
 Information: 
    
    The purpose of this script is to display various graphs, plots to show:
        anything that has to do with satellite trajectories

        
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import os   
os.environ['PROJ_LIB'] = r'C:\Users\Xavier\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi

import GH_convert     as conv
import GH_import      as imp
#import GH_generate    as gen
#import GH_solve       as solv
#import GH_displayCoef as dcoef



# =============================================================================
# DISPLAY FUNCTIONS  
# =============================================================================
def Plot2D_PosEarthfixed (fignum, Pos, Title="No given title"):
    """
    This function plots spherical coordinates on a basemp plot in a mpl figure
    Input:
        fignum: index of the matplotlib figure to be created
        Pos: Coordinates in spherical reherantial, radian degrees
        Title: Title of the plot to appear on the figure
    Output:
        FIG: matplotlib figure object created in this function
        (MAP: basemap plot created in this function)    
    """

#    Rs = Pos[:,0] # in km
    Lat = Pos[:,1] * 180/pi
    Long = Pos[:,2] * 180/pi
    
    FIG = plt.figure(fignum)    
    

    """Basemap parameters"""
    proj = "mill" # projection
    LatS = -90 # llcrnrlat
    LatN = 90 # urcrnrlat 
    LongW = -180 # llcrnrlon
    LongE = 180 # urcrnrlon
    TS = 20 # lat_ts -- I don't know what this is but everyone online uses it so yeah
    Res = "c" # resolution
    water_color = 'lightcyan'
    land_color = 'peachpuff'
#    parallels = np.arange(-60.,61,30.)
#    meridians = np.arange(0.,351.,30.)

    MAP = Basemap(projection = proj, 
                 llcrnrlat = LatS, 
                 urcrnrlat = LatN, 
                 llcrnrlon = LongW, 
                 urcrnrlon = LongE, 
                 lat_ts = TS, 
                 resolution = Res)
    
    """MAP details"""
    MAP.drawcoastlines()
#    MAP.fillcontinents(color=land_color,lake_color=water_color)
#    MAP.drawmapboundary(fill_color=water_color)
#    MAP.drawmeridians()
#    MAP.drawparallels()    
#    MAP.drawparallels(parallels)
#    MAP.drawmeridians(meridians)


    MAP.scatter(Long, Lat, s=0.1, c='r', latlon=True, alpha=1)
    plt.suptitle("Position projected on Earth")
    plt.title(Title)    
    
    plt.show(block=False)
    return FIG




def Plot3D_Pos (fignum, Pos, Title):    
    """
    This functions creates a matplotlib figure and plots Pos in 3D. 
    
    Input:
        Pos: array in spherical coordinates of the satellite position
        Title: The title to be displayed on the plot
    Output:
        fig: the index of the created matplotlib figure
        
# =============================================================================
# The spyder kernel keeps crashing when I run this. 
# I used to think it was a basemap issue but that isn't it.   
# =============================================================================
    """   
    X_1 = []
    Y_1 = []
    Z_1 = []
    zero = [] #for the center of the Referential
    zero.append(0)
    
    for i in range (0, len(Pos)):
        r, theta, phi = Pos[i]
        xi, yi, zi = conv.sph2cart(r, theta, phi)
        X_1.append(xi)
        Y_1.append(yi)
        Z_1.append(zi)
        
    FIG = plt.figure(fignum) 
    ax1 = FIG.add_subplot(111, projection = '3d') 
    
    ax1.plot3D( X_1, Y_1, Z_1, 
               'ro', markersize=0.5, alpha=0.5)
    
    ax1.plot3D(zero, zero, zero, 'b*')
    ax1.set_xlabel('x') 
    ax1.set_ylabel('y') 
    ax1.set_zlabel('z') 
    plt.suptitle("Position in the terrestrial referential")
    plt.title(Title)
    plt.show(block=False) 
    
    return FIG


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_Plots ():
    file_name = "Polar_400km_EarthFixed_1jour_1sec.e"
    days = 0.9
    Pos, Time = imp.Fetch_Pos(file_name, days)
    Title = f"Test_Plots: Earthfixed polar 400km orbit for {days} days"
    fig1 = Plot3D_Pos(1, Pos, Title)
    fig2 = Plot2D_PosEarthfixed(2, Pos, Title)
    return


# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    TEST_Plots()
    
    print("\nGH_displaySat done")

