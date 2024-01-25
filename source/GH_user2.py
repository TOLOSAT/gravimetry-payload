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
    puis comparer ces mêmes accéls à celles données en résolvant le système linéaire 
    donné par les coeffs du EGM2008 et accel_sg en résultat. Je vois pas trop le but...    
    
    TRUC A TESTER : géoïde avec coeffs calculés vs avec coeffs donnés...
"""
def main():
    path ="../data"
    file_name = "Polar_400km_EarthFixed_15jours_5sec.e"
    
    lmax = 10 # Legendre polynomial degree

    Pos,Vit,Time,dt = imp.Fetch_Pos_Vit(file_name, 5, path, False) # Read data from the file
    Pos_sph = conv.cart2sphA(Pos)
   
    # Uses Savitzky-Golay filter to get the radial acceleration
    acc = gen.Gen_Acc_2(Pos,Vit,dt) 
    acc = conv.cart2sphA(acc)
    acc = conv.Make_Line_acc(acc)    
    accRadial = [acc[3*i] for i in range(len(acc)//3)] # radial acceleration directly give by the SG filter


    # getting the radial acceleration from the full grad potential matrix and the full acceleration
    """Generate the grad Potential matrix"""
    getMat = lambda lmax : solv.Get_PotGradMatrix2(lmax, Pos_sph)
    M = getMat(lmax) # The full grad potential matrix
    Res = np.linalg.lstsq(M, acc) # Solve the linear system with the full matrix      
    acc_solved = M.dot(Res[0]) # acceleration from the solved system
    acc_solved_R = [acc_solved[3*i] for i in range(len(acc)//3)] # radial acceleration from the solved system
   
    # getting the radial acceleration from the partial grad potential matrix and the radial acceleration
    Mradial = GetPartialMatrix(M) # The partial grad potential matrix (only the radial part)  
    ResRadial = np.linalg.lstsq(Mradial, accRadial) # Solve the linear system with the partial matrix  
    accRadial_solved = M.dot(ResRadial[0]) # acceleration from the solved system   
    accRadial_solved_R = [accRadial_solved[3*i] for i in range(len(acc)//3)] # radial component of the acceleration
    
    # Reference data 
    #M = np.load(imp.data_path + "\potGradMatrix_Polar_400km_EarthFixed_15jours_5sec.npy")
    

    plt.figure()
    plt.plot(Time,acc_solved_R, label="Derived from reference data")
    plt.plot(Time[10:],accRadial_solved_R[10:], label = "Derived from partial reference data")
    plt.plot(Time[10:],accRadial[10:], label = "SG calculated data")
    plt.ylabel("Radial acceleration (m/s**2)")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.title("Comparaison des résultats")
    plt.show()
    
    return accRadial, acc_solved_R, accRadial_solved_R, M, Mradial



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
        















