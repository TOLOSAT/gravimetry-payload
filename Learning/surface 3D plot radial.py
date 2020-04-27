#%% -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 02:47:57 2020


plot sherical hermonics in 3D with radial colormap
http://balbuceosastropy.blogspot.com/2015/06/spherical-harmonics-in-python.html

"""

from __future__ import division
import scipy as sci
import scipy.special as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors


#%% ===========================================================================
l = 4    #degree
m = 2    # order
PHI, THETA = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j] #arrays of angular variables
R = np.abs(sp.sph_harm(m, l, PHI, THETA)) #Array with the absolute values of Ylm
"""
THETA = pi/2 - G_Lat*pi/180
PHI = G_Long*pi/180 + pi
R = G_Grid + 50000
"""
#Now we convert to cartesian coordinates
# for the 3D representation
X = R * np.sin(THETA) * np.cos(PHI)
Y = R * np.sin(THETA) * np.sin(PHI)
Z = R * np.cos(THETA)

N = R/R.max()    # Normalize R for the plot colors to cover the entire range of colormap.
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(7,5))
im = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(N))
ax.set_title(r'$|Y^2_ 4|$', fontsize=20)
m = cm.ScalarMappable(cmap=cm.jet)
m.set_array(R)    # Assign the unnormalized data array to the mappable
                  #so that the scale corresponds to the values of R
fig.colorbar(m, shrink=0.8);





#%% ===========================================================================
l = 4    # degree
m = 2    # order
PHI, THETA = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
R = sp.sph_harm(m, l, PHI, THETA).real

X = R * np.sin(THETA) * np.cos(PHI)
Y = R * np.sin(THETA) * np.sin(PHI)
Z = R * np.cos(THETA)

#As R has negative values, we'll use an instance of Normalize
#see http://stackoverflow.com/questions/25023075/normalizing-colormap-used-by-facecolors-in-matplotlib
norm = colors.Normalize()
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(7,5))
m = cm.ScalarMappable(cmap=cm.jet)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm(R)))
ax.set_title('real$(Y^2_ 4)$', fontsize=20)
m.set_array(R)
fig.colorbar(m, shrink=0.8);




#%% ===========================================================================
l = 4    # degree
m = 2    # order
PHI, THETA = np.mgrid[0:2*np.pi:300j, 0:np.pi:150j]
R = sp.sph_harm(m, l, PHI, THETA).real

s = 1
X = (s*R+1) * np.sin(THETA) * np.cos(PHI)
Y = (s*R+1) * np.sin(THETA) * np.sin(PHI)
Z = (s*R+1) * np.cos(THETA)

norm = colors.Normalize()
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(7,5))
m = cm.ScalarMappable(cmap=cm.jet)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.terrain(norm(R)))
ax.set_title('1 + real$(Y^2_ 4)$', fontsize=20)
m.set_array(R)
fig.colorbar(m, shrink=0.8);


#%%



