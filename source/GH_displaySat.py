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
#import os   
#os.environ['PROJ_LIB'] = r'C:\Users\Xavier\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

import GH_convert     as conv
import GH_import      as imp
import GH_generate    as gen
import GH_solve       as solv
import GH_displayCoef as dcoef
import GH_displaySat  as dsat


# =============================================================================
# DISPLAY FUNCTIONS  
# =============================================================================
def Plot2D_PosEarthfixed (Pos, Title):
    """
# =============================================================================
#   Clean this one up
    Where are the dots
    Clear out why creating a plot in a function pauses code
# =============================================================================
    """
    # Basemap parameters
    proj = "mill" # projection
    LatS = -90 # llcrnrlat
    LatN = 90 # urcrnrlat 
    LongW = -180 # llcrnrlon
    LongE = 180 # urcrnrlon
    TS = 20 # lat_ts -- I don't know what this is but everyone online uses it so yeah
    Res = "c" # resolution
    water_color = 'lightcyan'
    land_color = 'peachpuff'

    Xm = []
    Ym = []
    Rs = []
    
    for i in range(0, len(Pos)):  
        r, theta, phi = Pos[i]
        Rs.append(r)
        Xm.append(phi )#* 180/np.pi)
        Ym.append(theta)# * 180/np.pi)
        
    
    fig = plt.figure(1)
    m1 = Basemap(projection = proj, 
                 llcrnrlat = LatS, 
                 urcrnrlat = LatN, 
                 llcrnrlon = LongW, 
                 urcrnrlon = LongE, 
                 lat_ts = TS, 
                 resolution = Res)
    m1.drawcoastlines()
    m1.fillcontinents(color=land_color,lake_color=water_color)
    m1.drawmapboundary(fill_color=water_color)
    #m1.drawmeridians()
    #m1.drawparallels()


    
    m1.plot(Xm, Ym, 'ro', markersize=0.7, alpha=1)
    plt.title("Position projected on Earth")
    plt.suptitle(Title)    
    
    plt.show(1)
    return fig




def Plot3D_Pos (Pos, Title):    
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
        
    fig= plt.figure(2) 
    ax1 = fig.add_subplot(111, projection = '3d') 
    
    ax1.plot3D( X_1, Y_1, Z_1, 
               'ro', markersize=0.5, alpha=0.5)
    
    ax1.plot3D(zero, zero, zero, 'b*')
    ax1.set_xlabel('x') 
    ax1.set_ylabel('y') 
    ax1.set_zlabel('z') 
    plt.title("Position {r theta phi}, in terrestrial referential")
    plt.suptitle(Title)
    plt.show(2) 
    return fig


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_Plots ():
    file_name = "Polar_400km_EarthFixed_1jour_1sec.e"
    days = 0.9
    Pos, Time = imp.Fetch_Pos(file_name, days)
    Title = "Test Earthfixed plot for " + str(days)+ " days"
    fig1 = Plot3D_Pos(Pos, Title)
#    fig2 = Plot2D_PosEarthfixed(Pos, Title)



# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    TEST_Plots()
    
    print("\nGH_displaySat done")

