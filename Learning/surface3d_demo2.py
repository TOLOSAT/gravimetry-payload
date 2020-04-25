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
# Make data
u = np.linspace(0, 1.5* pi, 100)
v = np.linspace(0, pi, 100)
x = R * np.outer(cos(u), sin(v))
y = R * np.outer(sin(u), sin(v))
z = R * np.outer(np.ones(u.size), cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')

plt.show()
