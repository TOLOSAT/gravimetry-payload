"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to serve as the main code for the user.
    Use the user's guide that still has not been written to find out how to use
    this tool.
    Use your python skills to understand how this tool has been coded.

    Change your variables, do your maths, and do your science.
    Let's go to space, y'all

    Process :
    1.  Generate the acceleration values (from orbit data or raw simulation)
            using GH_generate
    2.  Solve for the coefficients using the functions in this script
            using GH_solve
    3.  Generate a Geoid map from solved coefficients
            using GH_displayCoef
    4.  Compare Geoids or coefficients
            using GH_displayCoef
    5.  Plot Satellite positions and accelerations in space
            using GH_displaySat


# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from time import gmtime, strftime

import GH_import       as imp
import GH_convert      as conv
import GH_generate     as gen
import GH_solve        as solv
import GH_displayGeoid as dgeo
import GH_displaySat   as dsat
import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_basemap      as bmp


from GH_import import data_path #= "../data"


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

# =============================================================================
# # ===========================================================================
# =============================================================================
    """ THE THINGS YOU CAN CHANGE AS A USER """

    """ The original satellite path """
    #file_name = "ISS_Earthfixed_1jour_1sec.e"
    #file_name = "ISS_Earthfixed_1jour_60sec.e"
#    file_name = "Polar_400km_EarthFixed_1jour_1sec.e"
    #file_name = "Polar_400km_EarthFixed_15jours_5sec.e"
    file_name = "Polar_400km_EarthFixed_7jours_5sec.e"
    days = 0.5

    """ data solving """
    lmax_gen = 5 # when generating the data

    lmax_solve = 5  # when solving for coefficients

    """ plotting maps of geoids """
    lmax_topo = 6
    tens = 1
    levels = 35

    """ save the plots and coefficients """
    save = False
    save_im_path = "../Rendered/images"
    save_co_path = "../Rendered/coefficients"

# =============================================================================
# # ===========================================================================
# =============================================================================

    time_str = imp.Get_Time()

    HC, HS = imp.Fetch_Coef()
    Pos_sim, Time = imp.Fetch_Pos(file_name, days)

    Acc_sim = gen.Gen_Sim_Acc(lmax_gen, HC, HS, Pos_sim)

    Solved_coef_sim, Acc_solved_sim = solv.Solve_Coef(lmax_solve, Pos_sim, Acc_sim)
    HC_sim, HS_sim = conv.Make_Array_Coef(lmax_solve, Solved_coef_sim)

    Acc_solved_sim = conv.Make_Array(Acc_solved_sim[:-2]) # this "-2" must be replaced with a modulo function to get the highest number thats a multiple of 3, AND smaller than the length of the array

    # plotting path simulation
    title1 = "Simulated and solved acceleration"
    component = 0 #0, 1, 2 : r, theta, phi
    FIG_ACC = dsat.Plot_Acc_Sim_Solv(1, Time, Acc_sim, Acc_solved_sim, component, title1)

    # Mapping the coefficients
    title2 = f"Map of original geoid"
    MAP_GEN = dgeo.Map_Geoid(2, lmax_gen,   HC,     HS,     tens, levels, title2, lmax_topo)

    title3 = f"Map of simulated geoid"
    MAP_SIM = dgeo.Map_Geoid(3, lmax_solve, HC_sim, HS_sim, tens, levels, title3, lmax_topo)


    if save:
        print("Saving plots and coefficients")
#        exp.Store_Array(HC_sim, f"HC_sim {lmax_gen} to {lmax_solve} degrees.txt", save_co_path)
#        exp.Store_Array(HS_sim, f"HS_sim {lmax_gen} to {lmax_solve} degrees.txt", save_co_path)

        exp.Store_Figure(FIG_ACC.number, title1, time_str, save_im_path, 500)
        exp.Store_Figure(MAP_GEN.number, title2, time_str, save_im_path, 500)
        exp.Store_Figure(MAP_SIM.number, title3, time_str, save_im_path, 500)


print("\nUser instructions done")

"""
# =============================================================================
# =============================================================================
# =============================================================================
# # #

Here's whats going wrong:
    The acceleration generation function returns weird values.
    Just look at "Acc_sim" to see that the values are way too large
    What units are we dealing with here ?
    The data diverges as the time samples get longer

Also the solve functions generate completey different coefficients
than the originals, and it seems to diverge when the orders go above 8

Also I have a huge doubt: am I confusing Acceleration with gravity potential?

https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen

# =============================================================================
# =============================================================================
# =============================================================================
"""