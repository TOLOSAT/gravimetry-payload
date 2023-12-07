"""
@authors:
# =============================================================================
 Information:
    The purpose of this script is to generate data arrays for satellite
    trajectory
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt

import GH_import       as imp
import GH_convert      as conv
import GH_Savitzky_Golay as sg
#import GH_generate     as gen
import GH_solve        as solv
#import GH_displayGeoid as dgeo
import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


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
    print("Generating simulated acclerations, lmax =", lmax, "")

    CS = conv.Make_Line_Coef(lmax, HC, HS)
    print(f"Shape of the Coef array = {CS.shape}")

    M_PotGrad = solv.Get_PotGradMatrix(lmax, Pos) # get M_PotGrad
#    print("shape of M=", M_PotGrad.shape)

    Acc_line = M_PotGrad.dot(CS)
#    print("shape of Acc_line=", Acc_line.shape)

    Acc_sim = conv.Make_Array(Acc_line, 3)
#    print("shape of Acc_sim=", Acc_sim.shape)

    return Acc_sim



def Gen_Acc_2(Pos,Vit,dt):
    x,y,z = Pos.T
    # Filter parameters : use of 20 values and cubic polynomial function
    ax = sg.savitzky_golay(x,20,3,2,dt)
    ay = sg.savitzky_golay(y,20,3,2,dt)
    az = sg.savitzky_golay(z,20,3,2,dt)

    Acc = np.array([ax,ay,az]).T
    Acc = sg.correcRef(Pos,Vit,Acc)

    return Acc

# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def Test_gen_acc ():
    file_name = "Polar_400km_EarthFixed_7jours_5sec.e"
    days = 0.03
    Pos, Time = imp.Fetch_Pos(file_name, days)
    Acc, Spe = Gen_Acc_2(Pos, Time)

    FIG = dsat.Plot_pos_spe_acc(Pos, Spe, Acc, Time)

    return

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

    Test_gen_acc ()

    print("\nGH_generate done")