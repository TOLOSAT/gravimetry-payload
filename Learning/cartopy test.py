# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:55:39 2020

@author: Xavier
"""

'''
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

from cartopy import config
import cartopy.crs as ccrs


# =============================================================================
tens = 5
size_long = 1 + 36*tens
size_lat  = 1 + 18*tens
Line_long = np.linspace(0, 2*pi, size_long) # 0 to 360 ; must subtract 180
Line_lat  = np.linspace(0, pi, size_lat) # 0 to 180 ; must do 90 - theta
G_Long, G_Lat = np.meshgrid((Line_long - pi), (pi/2 - Line_lat))

G_Grid = np.zeros((size_lat, size_long))

f = lambda x, y : x**2 + y**2

for i in range(0, size_long):
    long = G_Long[0][i]
    for j in range (0, size_lat):
        lat = G_Lat[j][0]
        G_Grid[j,i] = f(long,lat)


PROJ = ccrs.PlateCarree
PROJ = ccrs.Miller
PROJ = ccrs.Mollweide
#PROJ = ccrs.Miller


ax = plt.axes(projection=PROJ())


plt.contourf(G_Long*180/pi, G_Lat*180/pi, G_Grid, 60, alpha = 0.5,
             transform=ccrs.PlateCarree(), cmap="jet")


#ax.stock_img()
#ax.coastlines(resolution="110m")

plt.show()
'''

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
import numpy as np
from numpy import pi

def main():
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([80, 170, -45, 30])

    # Put a background image on for nice sea rendering.
#    ax.stock_img()



    ax.add_feature(cfeature.LAND, facecolor = np.array([ 0.859375, 0.859375, 0.859375]))
    ax.add_feature(cfeature.COASTLINE, edgecolor = 'black')
    ax.add_feature(cfeature.OCEAN, facecolor = np.array([ 0.59375 , 0.71484375, 0.8828125 ]))


    plt.show()


if __name__ == '__main__':
    main()














