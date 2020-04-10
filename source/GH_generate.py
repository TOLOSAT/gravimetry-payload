"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to generate data arrays.
    Satellite acceleration, Topology

    These arrays must be in spherical coordinates, in km.

# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from numpy import pi, sin, cos

import GH_import       as imp
import GH_convert      as conv
#import GH_generate     as gen
import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
import GH_terminal     as term
#import GH_basemap      as bmp

# =============================================================================
# FUNCTIONS TO GENERATE ACCELERATION ARRAYS
# =============================================================================
def Gen_Sim_Acc (lmax, HC, HS, Pos):
    """
    Generates simulated acceleration values from known coefficients
    and known positions

    Input:
        lmax: max order desired
        HC: cosine coefficients array
        HS: sine coefficients array
        Pos: line
    Output:
        Acc_sim: simulated acceleration values in spherical coordinates

    """
    print("\nGenerating simulated acclerations, lmax =", lmax, "")

    CS = conv.Make_Line_Coef(lmax, HC, HS)
    print(f"Shape of the Coef array = {CS.shape}")

    M_PotGrad = solv.Get_PotGradMatrix(lmax, Pos) # get M_PotGrad
#    print("shape of M=", M_PotGrad.shape)

    Acc_line = M_PotGrad.dot(CS)
#    print("shape of Acc_line=", Acc_line.shape)

    Acc_sim = conv.Make_Array(Acc_line, 3)
#    print("shape of Acc_sim=", Acc_sim.shape)

    return Acc_sim


# =============================================================================
# FUNCTIONS TO WORK ON TOPOLOGY
# =============================================================================
def Get_Topo_Height (lmax, Lat, Long, HC_topo, HS_topo):
    """
    This function returns the height of Earth's estimated topology at Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    """
    Sum1 = 0
    Pmn, _ = imp.Pol_Legendre(lmax, lmax, cos(Lat)) # I am allowed to write that.

    for l in range(0, lmax+1):
        Sum2 = 0

        for m in range (0,l+1):
            Sum2 += (HC_topo[l,m]*cos(m*Long) + HS_topo[l,m]*sin(m*Long)) * Pmn[m, l] * imp.Normalize(l, m)

        Sum1 += Sum2

    return Sum1

# =============================================================================
# FUNCTIONS TO WORK ON GEOID
# =============================================================================
def Get_Geo_Pot (lmax, R, Lat, Long, HC, HS):
    """
    This function returns the potential at given height/Lat/Long coordinates
    The solution is calculated up to degree lmax in the HC HS model
    """
    GM = 3986004.415E8 # m**3 s**-2 : Earth's standard gravitational parameter
    # wiki says : gm = 6.673*10**-11 * 5.975*10**24 = 398711749999999.94 OR 3.986004418E14
    a = 6378136.3 # m

    Sum1 = 0
    Pmn, _ = imp.Pol_Legendre(lmax, lmax, cos(Lat))

    for l in range (2, lmax+1):
        Sum2 = 0

        for m in range(0, l+1):
            Sum2 += (HC[l,m]*cos(m*Long) + HS[l,m]*sin(m*Long)) * Pmn[m, l] * imp.Normalize(l, m)

        Sum1 += (a/R)**l * Sum2 # OK, the "a" here is weird, idk what it's doing here

    Pot = GM/R*(1 + Sum1)

    return Pot





# =============================================================================
# FUNCTIONS TO GENERATE DATA ARRAYs
# =============================================================================
def Gen_Grid (measure, lmax, HC, HS, tens, lmax_topo=10, HC_topo=[], HS_topo=[]):
    """
    This function generates an array containing Earth's topology at Lat/Long
    coordinates.
    Input:
        measure: "topo", "geopot" or "geoid", the measurement to be mapped on Earth
        lmax: degree to which the topology should be calculated
        HC_topo: Harmonic cosine coefficients to earth's topology
        HS_topo: Harmonic sine coefficients to earth's topology
        tens: how large the array should be
    Output:
        G_Height: array of grid height
        G_long: grid of longitudes
        G_lat: grid of latitudes

    Additions:
        A progress bar for when there are a lot of points
        https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    if (measure == "geopot"):
        HC_topo, HS_topo = imp.Fetch_Topo_Coef()

    size_long = 1 + 36*tens
    size_lat  = 1 + 18*tens
    points = size_long * size_lat

    print(f"Generating {measure} grid for lmax = {lmax}, {points} points")

    Line_long = np.linspace(0, 2*pi, size_long) # 0 to 360 ; must subtract 180
    Line_lat = np.linspace(0, pi, size_lat) # 0 to 180 ; must do 90 - theta
    G_Long, G_Lat = np.meshgrid((Line_long - pi), (pi/2 - Line_lat))

    G_Grid = np.zeros((size_lat, size_long))

    for i in range(0, size_long):
        term.printProgressBar(i+1, size_long)
        Long = Line_long[i]

        for j in range(0, size_lat):
            Lat = Line_lat[j]

            if (measure == "topo"):
                G_Grid[j,i] = Get_Topo_Height(lmax, Lat, Long, HC, HS)
            elif (measure == "geopot"):
                R = imp.Get_Radius(Lat) + Get_Topo_Height(lmax_topo, Lat, Long, HC_topo, HS_topo) # add Earth's radius !!
                G_Grid[j,i] = Get_Geo_Pot(lmax, R, Lat, Long, HC, HS)


    return G_Grid, G_Long*180/pi, G_Lat*180/pi # in degrees now



# =============================================================================
# TEST FUNCTIONS
# =============================================================================

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

    print("\nGH_generate done")

