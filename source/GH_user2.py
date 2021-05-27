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

todo: replace "acc" by "geopot" lol

future stuff:
    earth rotation
    approx drag and shit (solar radiation, tidal effects)
    
Instruction after executing the code : here are the commands to execute on the terminal : 
    - A = getMat(n), with n a natural integer : we advise n=7.
    - res = np.linalg.lstsq(A, acc)
    - acc_solved = A@res[0] : acc_solved represents the vector of the solved acceleration in spherical coordinates. 
    Then acc_solved = [acc_rad(0),acc_theta(0),acc_phi(0),...]. Then plot the result.
    If you want to plot for example only the radial acceleration, execute this command : 
    - plt.plot([acc_solved[3*i] for i in range(len(acc_solved)//3)])
   
# =============================================================================
"""





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
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


from GH_import import data_path #= "../data"


file_name = "Polar_400km_EarthFixed_15jours_5sec.e"
lmax = 10



Pos,Vit, t = imp.Fetch_Pos_Vit(file_name, 5, spherical=False)

Pos_sphere = conv.cart2sphA(Pos)

acc = gen.Gen_Acc_2(Pos,Vit,t)
acc = conv.cart2sphA(acc)

acc = conv.Make_Line_acc(acc)
getMat = lambda lmax : solv.Get_PotGradMatrix2(lmax, Pos_sphere)

hc, hs = imp.Fetch_Coef()

hc = hc.flatten()
hc = np.sort(hc)


























