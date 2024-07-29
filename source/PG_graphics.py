import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

import PG_maths
import PG_globalVars

h_order = 5             #harmonic order
num = 180               #number of points for the grid

#########################################################################################
#                                   3D Grids functions                                  #
#########################################################################################


def GeoPot_Grid(coeffs):

    """ parameters :
            - coeffs : harmonics coefficients

        output :
            - x,y,z grids for earth representation
            - geopontential grid """

    # definition and construction of spherical coords. grids
    theta = np.linspace(0, np.pi, num)
    phi = np.linspace(-np.pi, np.pi, num)
    theta, phi = np.meshgrid(theta, phi)

    # definition of the radius and geopot. grids
    pot_grid = np.zeros((num, num))
    r = np.zeros((num, num))

    # construction of the grids
    for i in range(num):
        for j in range(num):
            r[j][i], pot_grid[j][i] = PG_maths.GeoPot(coeffs, theta[0][i], phi[j][0], h_order)

    # definiton and construction of the x,y,z grids
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z, pot_grid


def GeoPot_Diff_Grid(coeffs):

    """ parameters :
            - coeffs : harmonics coefficients

        output :
            - x,y,z grids for earth representation
            - geopontential anomalies grid """

    # definition and construction of spherical coords. grids
    theta = np.linspace(0, np.pi, num)
    phi = np.linspace(-np.pi, np.pi, num)
    theta, phi = np.meshgrid(theta, phi)

    # definition of the radius and geopot. grids
    pot_grid = np.zeros((num, num))
    r = np.zeros((num, num))

    # construction of the grids
    for i in range(num):
        for j in range(num):
            r[j][i], pot_grid[j][i] = PG_maths.GeoPotDiff(coeffs, theta[0][i], phi[j][0], h_order)

    # definiton and construction of the x,y,z grids
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z, pot_grid


def Geoid_Grid(coeffs, scaling_factor=10000):

    """ parameters :
            - coeffs : harmonics coefficients

        output :
            - x,y,z grids for earth representation
            - geoid height anomalies grid """

    # definition and construction of spherical coords. grids
    theta = np.linspace(0, np.pi, num)
    phi = np.linspace(-np.pi, np.pi, num)
    theta, phi = np.meshgrid(theta, phi)

    # definition of the radius and geopot. grids
    pot_grid = np.zeros((num, num))
    r = np.zeros((num, num))

    # construction of the grids
    for i in range(num):
        for j in range(num):
            r[j][i], pot_grid[j][i] = PG_maths.GeoPotDiff(coeffs, theta[0][i], phi[j][0], h_order)
            print(r[j][i], pot_grid[j][i])

    # standard gravity for WGS84 model (Helmert's equation)
    g = 9.780327 * (1 + 0.0053024 * np.sin(np.pi/2 - theta)**2 - 0.0000058 * np.sin(np.pi - 2 * theta)**2)

    # definition and construction of the x,y,z grids
    x = (r + pot_grid/g * scaling_factor) * np.sin(theta) * np.cos(phi)
    y = (r + pot_grid/g * scaling_factor) * np.sin(theta) * np.sin(phi)
    z = (r + pot_grid/g * scaling_factor) * np.cos(theta)

    return x, y, z, pot_grid/g

#########################################################################################
#                                           3D Figures                                  #
#########################################################################################


def GeoPot_3D(coeffs, axis=True):

    """ parameters :
            - coeffs : harmonics coefficients
            - axis : to display or not the axis

        output :
            - x,y,z grids for earth representation
            - geoid height anomalies grid """

    x, y, z, pot_grid = GeoPot_Grid(coeffs)

    cmap = plt.get_cmap('jet')
    norm = mcolors.Normalize()
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.plot_surface(x, y, z, facecolors=cmap(norm(pot_grid)), linewidth=0, antialiased=False, rcount=80, ccount=80)
    cbar = plt.colorbar(mappable=m, ax=ax, cmap=cmap, orientation='horizontal', pad=0.10)
    cbar.set_label("Gravitational potential in m^2/s^2")

    if not axis:
        plt.axis("off")

    plt.show()


def GeoPotDiff_3D(coeffs, axis=True):

    x, y, z, pot_grid = GeoPot_Diff_Grid(coeffs)

    cmap = plt.get_cmap('jet')
    norm = mcolors.Normalize()
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.plot_surface(x, y, z, facecolors=cmap(norm(pot_grid)), linewidth=0, antialiased=False, rcount=80, ccount=80)
    cbar = plt.colorbar(mappable=m, ax=ax, cmap=cmap, orientation='horizontal', pad=0.10)
    cbar.set_label("Gravitational potential anomalies in m^2/s^2")

    if not axis:
        plt.axis("off")

    plt.show()


def Geoid_scaled(scaling_factor=10000, axis=True):

    x, y, z, h_diff = Geoid_Grid(scaling_factor)

    cmap = plt.get_cmap('jet')
    norm = mcolors.Normalize()
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.plot_surface(x, y, z, facecolors=cmap(norm(h_diff)), linewidth=0, antialiased=False, rcount=80, ccount=80)
    cbar = plt.colorbar(mappable=m, ax=ax, cmap=cmap, orientation='horizontal', pad=0.10)
    cbar.set_label("geoid height in m (anomalies scaling factor : 10 000", )

    if not axis:
        plt.axis("off")

    plt.show()
