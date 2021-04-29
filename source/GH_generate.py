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



def Gen_Acc_2(Pos,Vit,t):
    x,y,z = Pos.T
    ax = sg.savitzky_golay(x,20,3,1,t[1]-t[0])
    ay = sg.savitzky_golay(y,20,3,1,t[1]-t[0])
    az = sg.savitzky_golay(z,20,3,1,t[1]-t[0])

    Acc = np.array([ax,ay,az]).T


    Acc = sg.correcRef(Pos,Vit,Acc)



    return Acc












def Gen_Acc (Pos, t):


    """
    calculates the acceleration of each axis using a basic doudle derivation
    Input:
        Pos: trail of postion
        t: time stamp of Pos
    Output:
        Acc: guestimated acceleration values in spherical coordinates
    """
    print("Guestimating acclerations")

    r, th, ph = Pos[:,0], Pos[:,1], Pos[:,2]
    Spe = np.zeros((len(r)-1, 3))
    Acc = np.zeros((len(r)-2, 3))

    for i in range (0, len(t)-1):
        dt = t[i+1] - t[i]
        dr_dt = (r[i+1] - r[i] ) / dt
        dtheta_dt = (th[i+1] - th[i] ) / dt
        dphi_dt = (ph[i+1] - ph[i] ) / dt

        Spe[i,:] = [dr_dt, dtheta_dt, dphi_dt]

    dr, dth, dph = Spe[:,0], Spe[:,1], Spe[:,2]
    for i in range (0, len(t)-2):
        dt = t[i+1] - t[i]
        d2r_dt = (dr[i+1] - dr[i] ) / dt
        d2theta_dt = (dth[i+1] - dth[i] ) / dt
        d2phi_dt = (dph[i+1] - dph[i] ) / dt

        Acc[i,:] = [d2r_dt, d2theta_dt, d2phi_dt]

    return Acc, Spe


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def Test_gen_acc ():
    file_name = "Polar_400km_EarthFixed_7jours_5sec.e"
    days = 0.03
    Pos, Time = imp.Fetch_Pos(file_name, days)
    Acc, Spe = Gen_Acc(Pos, Time)

    FIG = dsat.Plot_pos_spe_acc(Pos, Spe, Acc, Time)

    return

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

    Test_gen_acc ()

    print("\nGH_generate done")