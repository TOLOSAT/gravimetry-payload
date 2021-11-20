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
   
Here is the link to the files we use for this code : https://drive.google.com/drive/folders/1XgGn2QoFGJ-u_m4aoL2No-PmxIRraLg6?fbclid=IwAR04xPEqi1h-eGWR_inJfHa4dp8dzG7NPlgHcNKYfz0OT1v0uU7ADew9VR4
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





def GetPartialMatrix(M, comp = 0):
    """ Take only the rows corresponding
    to the Radial part
        - Input : Gravimetry potential matrix
                     Component = 0 , 1 , 2
        
        - Output Gravimetry potential with : - radial part if comp == 0
                                             - theta part if comp == 1
                                             - phi part if comp == 2
        """
    assert type(comp) ==int , "Comp is not an Integer"
    
    return np.array([[M[3*i + comp ,j ] for j in range(len(M[1])) ] for i in range(len(M)//3)])

def main():
    
    file_name = "Polar_400km_EarthFixed_15jours_5sec.e"
    
    """ Define the Legendre polynom degree """
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
    
<<<<<<< Updated upstream
    """Generate the grad Potential matrix"""
    #M = getMat(lmax)
    
    M = np.load(imp.data_path + "\potGradMatrix_Polar_400km_EarthFixed_15jours_5sec.npy")
=======

    
    #M = getMat(lmax)
    
    M = np.load(imp.data_path + "/potGradMatrix_Polar_400km_EarthFixed_15jours_5sec.npy")
    
>>>>>>> Stashed changes
    
    Mradial = GetPartialMatrix(M)
    
    Res = np.linalg.lstsq(M, acc)
    
    accRadial = [acc[3*i] for i in range(len(acc)//3)]
    ResRadial = np.linalg.lstsq(Mradial, accRadial)
    
    acc_solved = M.dot(Res[0])
    
    accRadial_solved = M.dot(ResRadial[0])
    
    acc_solved_R = [acc_solved[3*i] for i in range(len(acc)//3)]
    
    accRadial_solved_R = [accRadial_solved[3*i] for i in range(len(acc)//3)]
    
    plt.figure()
    #plt.plot(acc_solved_R, label="Full Order Matrix")
    plt.plot(accRadial_solved_R[10:], label = "Radial Part Matrix")
    plt.plot(accRadial[10:], label = "theoretical")
    plt.legend()
    plt.show()
    
    return accRadial, acc_solved_R, accRadial_solved_R, M, Mradial



if __name__ == "__main__":
    res = main()




















