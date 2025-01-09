import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.style
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from statistics import mean

import PG_maths
import PG_globalVars as gv

h_order = 5     # harmonic order
num = 400       # number of points for the grid


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
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)

    return x, y, z, pot_grid


def Geoid_Grid(coeffs, scaling_factor):
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
            # print(r[j][i], pot_grid[j][i])

    # standard gravity for WGS84 model (Helmert's equation)
    g = 9.780327 * (1 + 0.0053024 * np.sin(np.pi / 2 - theta) ** 2 - 0.0000058 * np.sin(np.pi - 2 * theta) ** 2)

    # definition and construction of the x,y,z grids
    x = (r + pot_grid / g * scaling_factor) * np.sin(theta) * np.cos(phi)
    y = (r + pot_grid / g * scaling_factor) * np.sin(theta) * np.sin(phi)
    z = (r + pot_grid / g * scaling_factor) * np.cos(theta)

    return x, y, z, pot_grid / g


def Coastlines_Grid(coeffs, scaling_factor):

    m = Basemap()
    poly_paths = m.drawcoastlines().get_paths()
    plt.close()
    lons, lats = [], []
    for i in range(91):
        poly_path = poly_paths[i]
        coords = np.array([(vertex[0], vertex[1]) for (vertex, code) in poly_path.iter_segments(simplify=False)])
        lon, lat = m(coords[:, 0], coords[:, 1], inverse=True)
        lats.extend(lat.tolist() + [None])
        lons.extend(lon.tolist() + [None])

    lats = np.array(lats, dtype=np.float64) * np.pi / 180
    lons = np.array(lons, dtype=np.float64) * np.pi / 180

    r = np.zeros(len(lats))

    for i in range(len(lats)):
        g = 9.780327 * (1 + 0.0053024 * np.sin(lats[i]) ** 2 - 0.0000058 * np.sin(lats[i]) ** 2)
        r_e, pot = PG_maths.GeoPotDiff(coeffs, np.pi / 2 - lats[i], lons[i], h_order)
        r[i] = 1.005 * r_e + pot / g * scaling_factor

    xs = r * np.cos(lons) * np.cos(lats)
    ys = r * np.sin(lons) * np.cos(lats)
    zs = r * np.sin(lats)

    return xs, ys, zs

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

    ax.plot_surface(x, y, z, facecolors=cmap(norm(pot_grid)), linewidth=0, antialiased=False, rcount=100, ccount=100)
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

    ax.plot_surface(x, y, z, facecolors=cmap(norm(pot_grid)), linewidth=0, antialiased=False, rcount=100, ccount=100)
    cbar = plt.colorbar(mappable=m, ax=ax, cmap=cmap, orientation='horizontal', pad=0.10)
    cbar.set_label("Gravitational potential anomalies in m^2/s^2")

    if not axis:
        plt.axis("off")

    plt.show()


def Geoid_scaled(coeffs, scaling_factor=10000, axis=True):

    x, y, z, h_diff = Geoid_Grid(coeffs, scaling_factor)

    cmap = plt.get_cmap('jet')
    norm = mcolors.Normalize()
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.plot_surface(x, y, z, facecolors=cmap(norm(h_diff)), linewidth=0, antialiased=False, rcount=200, ccount=200)
    cbar = plt.colorbar(mappable=m, ax=ax, cmap=cmap, orientation='horizontal', pad=0.10)
    cbar.set_label("geoid height in m (anomalies scaling factor : 10 000", )

    if not axis:
        plt.axis("off")

    plt.show()


def Geoid_scaled2(coeffs, scaling_factor=10000, coastlines=True):

    x, y, z, h_diff = Geoid_Grid(coeffs, scaling_factor)
    data = [go.Surface(x=x, y=y, z=z, surfacecolor=h_diff, colorscale='jet', opacity=1)]

    if coastlines:

        xs, ys, zs = Coastlines_Grid(coeffs, scaling_factor)

        data.append(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines', line=dict(color='black', width=1)))

    fig = go.Figure(data)
    fig.show()


