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
from numpy import pi, sin, cos, tan

import GH_convert     as conv
import GH_import      as imp
import GH_solve       as solv
#import GH_displayCoef as dcoef
#import GH_displaySat  as dsat
#import GH_export      as exp
#import GH_displayTopo as dtopo

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
    This function returns the height of Earth's etimated topology at Lat Long coordinates
    """
    Sum1 = 0
    Pmn, _ = imp.Pol_Legendre(lmax, lmax, cos(Lat)) # I am allowed to write that. 
    
    for l in range(0, lmax):
        Sum2 = 0
        for m in range (0,l):
            Sum2 = Sum2 + imp.Normalize(l, m) * Pmn[m, l] * (HC_topo[l,m]*cos(m*Long) + HS_topo[l,m]*sin(m*Long))
        
        Sum1 = Sum1 + Sum2

    return Sum1


def Gen_Topo (lmax, HC_topo, HS_topo, tens):
    """
    This function generates an array containing Earth's topology at Lat/Long 
    coordinates. 
    Input:
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
    size_long = 1 + 36*tens
    size_lat  = 1 + 18*tens
    points = size_long * size_lat
    
    print(f"Generating topological grid for lmax = {lmax}, {points} points")

    Line_long = np.linspace(0, 2*pi, size_long) # 0 to 360 ; must subtract 180
    Line_lat = np.linspace(0, pi, size_lat) # 0 to 180 ; must do 90 - theta
    G_Long, G_Lat = np.meshgrid((Line_long - pi), (pi/2 - Line_lat))   

    G_Height = np.zeros((size_lat, size_long))
    
    for i in range(0, len(Line_long)) :
        Long = Line_long[i]

        for j in range(0, len(Line_lat)) :
            Lat = Line_lat[j]
            G_Height[j,i] = Get_Topo_Height(lmax, Lat, Long, HC_topo, HS_topo)

    return G_Height, G_Long*180/pi, G_Lat*180/pi


    
# =============================================================================
# TEST FUNCTIONS
# =============================================================================
    
# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    
    print("\nGH_generate done")

