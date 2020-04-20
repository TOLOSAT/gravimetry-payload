# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:55:39 2020

@author: Xavier
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
import numpy as np

# =============================================================================
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([80, 170, -45, 30])
ax.stock_img()

ax.add_feature(cfeature.LAND, facecolor = cfeature.COLORS["land"])
ax.add_feature(cfeature.COASTLINE, edgecolor = 'black')
ax.add_feature(cfeature.OCEAN, facecolor = cfeature.COLORS["water"])

SOURCE = 'Natural Earth'
LICENSE = 'public domain'
text = AnchoredText(r'$\mathcircled{{c}}$ {}; license: {}'.format(SOURCE, LICENSE),
                    loc=4, prop={'size': 12}, frameon=True)
ax.add_artist(text)
plt.show()


# =============================================================================
lon = np.linspace(-80, 170, 25)
lat = np.linspace(-20, 70, 25)
lon2d, lat2d = np.meshgrid(lon, lat)
data = np.cos(np.deg2rad(lat2d) * 4) + np.sin(np.deg2rad(lon2d) * 4)

proj_crs = ccrs.Mollweide
data_crs = ccrs.PlateCarree

fig = plt.figure()
ax = plt.subplot(111, projection=proj_crs())
ax.set_global()
ax.coastlines()

ax.contourf(lon, lat, data, levels = 50, transform=data_crs(), cmap = "jet")
#ax.pcolormesh(lon, lat, data, transform=data_crs(), cmap = "jet")
plt.show()



