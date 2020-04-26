'''
========================
3D surface (solid color)
========================

Demonstrates a very basic plot of a 3D surface using a solid color.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin, cos

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
R = 10
size = 100
# Make data
u = np.linspace(0, 1.5* pi, size)
v = np.linspace(0, pi, size)
x = R * np.outer(cos(u), sin(v))
y = R * np.outer(sin(u), sin(v))
z = R * np.outer(np.ones(u.size), cos(v))

# Plot the surface
cbar = np.zeros((3,size))
cbar[0,:] = np.linspace(0, 1, size)
cbar[1,:] = np.linspace(0, 1, size)
cbar[2,:] = np.linspace(0, 1, size)


ax.plot_surface(x, y, z, color=cbar)

plt.show()
