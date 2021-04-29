"""
@authors:
# =============================================================================
 Information:
    The purpose of this script is to display various graphs and maps about
    Satellite trajectories
topo: implement earth rotation. precession ?
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

import numpy as np
from numpy import cos, sin, pi

import GH_import       as imp
import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
import GH_earthMap     as emap


# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================
def Plot2D_PosEarthfixed (Pos, Title="No given title"):
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

    FIG, AX = emap.Make_Map()
    AX.plot(Long, Lat, '.', markersize=0.5, transform=ccrs.PlateCarree())
#    AX.scatter(Long, Lat, s=0.1, c='r', latlon=True, alpha=1)
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
        FIG: matplotlib figure object created in this function
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
    plt.clf()

    ax1 = FIG.add_subplot(111, projection = '3d')

    ax1.plot3D( X_1, Y_1, Z_1,
               'ro-', markersize=0.5, alpha=1)

    ax1.plot3D(zero, zero, zero, 'b*')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    plt.suptitle("Position in the terrestrial referential")
    plt.title(Title)
    plt.show(block=False)

    return FIG


def Plot_Acc_Sim_Solv (Time, Acc_sim, Acc_solved, component, title):
    """
    Plot path acceleration of simuated and solved paths
    Input:
        fignum: index of the matplotlib figure to be created
        Time: array of time sampling of the coordinates
        Acc_sim: accelerations from the original harmonics model
        Acc_solved: accelerations from the solved harmonics model
        component = 0, 1, 2 : r, theta, phi
        title: plot title
    """
    FIG = plt.figure()
    plt.clf()

    plt.title(title)
    plt.xlabel("time (s)")
    plt.ylabel("acceleration (?)")

    plt.plot(Time[:Acc_sim.shape[0]], Acc_sim[:,component], "ro-", alpha=0.3, label="simulated")

    plt.plot(Time[:Acc_solved.shape[0]], Acc_solved[:,component], "bo-", alpha=0.3, label="solved")

    plt.legend()
    plt.show(block=False)

    return FIG


def Plot_pos_spe_acc(Pos, Spe, Acc, Time):
    """ Makes a 3x3 set of plots to show all axes of position, speed, and acceleration
    """
    FIG = plt.figure(figsize=(12,5))

    # r
    AX1 = FIG.add_subplot(331)
    plt.plot(Time[:], Pos[:,0], "b")
    plt.title("Altitude ")
    plt.ylabel("r (km)")
    AX1.set_xticklabels('')

    AX2 = FIG.add_subplot(332)
    plt.plot(Time[0:-1], Spe[:,0], "g")
    plt.title("Speed ")
    AX2.set_xticklabels('')

    AX3 = FIG.add_subplot(333)
    plt.plot(Time[0:-2], Acc[:,0], "r")
    plt.title("Acceleration")
    AX3.set_xticklabels('')

    # theta
    AX4 = FIG.add_subplot(334)
    plt.plot(Time[:], Pos[:,1], "b")
    plt.ylabel("theta")
    AX4.set_xticklabels('')

    AX5 = FIG.add_subplot(335)
    plt.plot(Time[0:-1], Spe[:,1], "g")
    AX5.set_xticklabels('')

    AX6 = FIG.add_subplot(336)
    plt.plot(Time[0:-2], Acc[:,1], "r")
    AX6.set_xticklabels('')

    # phi
    AX7 = FIG.add_subplot(337)
    plt.plot(Time[:], Pos[:,2], "b")
    plt.ylabel("phi")

    AX8 = FIG.add_subplot(338)
    plt.plot(Time[0:-1], Spe[:,2], "g")
    plt.xlabel("time (s)")

    AX9 = FIG.add_subplot(339)
    plt.plot(Time[0:-2], Acc[:,2], "r")

    return FIG


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_Plots ():
    file_name = "Polar_400km_EarthFixed_1jour_1sec.e"
    # file_name = "ISS_EarthMJ2000Eq_15jours_60sec.e"
    days = 15
    Pos, Time = imp.Fetch_Pos(file_name, days)
    Title = f"Test_Plots: {file_name} for {days} days"
#    fig1 = Plot3D_Pos(1, Pos, Title)
    fig2 = Plot2D_PosEarthfixed(Pos, Title)


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_Plots()

    print("\nGH_displaySat done")