# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:55:39 2020

@author: Xavier
"""
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


ax = plt.axes(projection=ccrs.PlateCarree())

plt.contourf(G_Long*180/pi, G_Lat*180/pi, G_Grid, 60, alpha = 0.5,
             transform=ccrs.PlateCarree(), cmap="jet")

ax.coastlines()

plt.show()

















