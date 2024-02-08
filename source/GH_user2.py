"""

@authors:

# =============================================================================
 Information:

    The purpose of this script is to serve as the main code for the user.
    Use the user's guide that still has not been written to find out how to use
    this tool.
    Use your python skills to understand how this tool has been coded.

    Change your variables, do your maths, and do your science.
    Let's go to space, y'all.

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


from GH_import import data_path #= "..\data"

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

"""
    Là on va utiliser Savitsky-Golay pour déduire des 
    accélérations (accel_sg) à partir de positions (pas sûr du process, on est en sphériques...)
    TRUC A TESTER : géoïde avec coeffs calculés vs avec coeffs donnés...
"""
def main():
    path ="/home/mehdi/Dev/Tolosat/g ravimetry-payload/data"
    file_name = "Polar_400km_EarthFixed_15jours_5sec.e"
    
    lmax = 10 # Legendre polynomial degree

    Pos,Vit,Time,dt = imp.Fetch_Pos_Vit(file_name, 5, path, False) # Read data from the file
    Pos_sph = conv.cart2sphA(Pos)
   
    # Uses Savitzky-Golay filter to get the acceleration
    acc = gen.Gen_Acc_2(Pos,Vit,dt) 
    acc = conv.cart2sphA(acc)
    acc = conv.Make_Line_acc(acc)

    """Generate the grad Potential matrix"""
    getMat = lambda lmax : solv.Get_PotGradMatrix2(lmax, Pos_sph)
    M = getMat(lmax)
   
    # Reference data 
    #M = np.load(imp.data_path + "\potGradMatrix_Polar_400km_EarthFixed_15jours_5sec.npy")



if __name__ == "__main__":
    res = main()





def inv(A,y, epsilon, nIter, x_0, mu, pas = 0.5):
    x = x_0[:]
    
    
    L = []
    
    err = 10.*epsilon + 1
    
    for i in range(nIter):
        if( err < epsilon ):
            return x
        x = x - pas*(2.* A.T@( A@x - y) + mu * ( x - x_0))
        err = np.linalg.norm(A@x - y)
        L.append(err)
    
    plt.plot(L)
        
    return x
        















