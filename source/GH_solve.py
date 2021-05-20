"""
@authors:
# =============================================================================
 Information:
    The functions in this script are used to solve to the spherical harmonic
    Stokes coefficients.
todo: write acc_matrix generation AND geopot_matrix generation
# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
import numpy.linalg as npl
from numpy import sin, cos

import GH_import       as imp
import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
import GH_terminal     as term
#import GH_harmonics    as harm
import GH_geoMath      as gmath
#import GH_earthMap     as emap


# =============================================================================
# FUNCTIONS FOR Sph Harm SOLVE
# =============================================================================



def Get_PotGradMatrix2 (lmax, Pos): #"R = 6378136.3 m):
    """
    Returns the matrix of the gravitational potential gradient.
    Watch out, it gets big fast.
    Multiplying it with the appropriate column vector of coefficients
    will return the acceleration at the given coordinates.
        *There are no geoid coefficients for l=0, l=1*
        *There are no sine coefficients for m=0*
    Input:
        lmax: max order
        Pos: array of N_points positions in spherical coordinates (r, theta, phi)
        *R: Reference radius in meters*
    Output:
        M_PotGrad: the matrix of the coefficients
    """
    # constants
    R = 6378.1363 # km
    GM = 398600.4418 # km**3 s**-2
    # wiki says : gm = 6.673*10**-11*5.975*10**24 = 398711749999999.94

    N_points = len(Pos) # number of points

    Cos_len = int( (lmax+1)*(lmax+2) /2 ) -3 # c20,c21,c22,c30,c31,c32, ...
    Sin_len = int( (lmax  )*(lmax+1) /2 )  # s21,s22,s31,s32,s33, ...
#    print("cos sin lengths =",Cos_len, ",",Sin_len)
    N_coef = Cos_len + Sin_len

    M_PotGrad = np.ones((N_points * 3, N_coef)) # THE Potential Gradient Matrix
    print(f"Generating BAM of shape = {M_PotGrad.shape}") # BAM =  "Big Ass Matrix"

    for i in range (0, N_points):
        term.printProgressBar(i+1, N_points)

        r, theta, phi = Pos[i] #spherical coordinates at the first point
        Plm_z, Plm_dz = gmath.Pol_Legendre(lmax, lmax, sin(theta))

        j = 0
        k = Cos_len

        for l in range (0, lmax +1):
            for m in range (0, l +1):
                # These equations were found in the GFZ document page 23
                W_r = - GM/r**2 * (R/r)**l * (l+1) * Plm_z[m, l]
                W_phi = -W_r * m * r / (l+1)
                W_theta = GM/r * (R/r)**l * Plm_dz[m, l]

                Sub_mat = np.zeros ((3,1))
                Sub_mat = [ cos(m*phi)*W_r,
                            cos(m*phi)*W_theta,
                            -sin(m*phi)*W_phi] # multiply by: COS_lm_coef

                M_PotGrad [3*i : 3*(i+1), j] = Sub_mat
                j += 1

                # for m of non-null, we get a sine coefficient
                if (m != 0):
                    Sub_mat = np.zeros ((3,1))
                    Sub_mat = [ sin(m*phi)*W_r,
                                sin(m*phi)*W_theta,
                                cos(m*phi)*W_phi] # multiply by: SIN_lm_coef

                    M_PotGrad [3*i : 3*(i+1), k] = Sub_mat
                    k += 1

    return M_PotGrad


def proj_ort(M,y):
    return y - M@np.linalg.lstsq(M,y)[0]




def Solve_Coef (lmax, Pos, Acc):
    """
    Returns the solved for coefficients to the spherical harmonic approximation
    of the Acc accelerations at Pos positions. Uses the least square methods
    Input:
        lmax: maximum degree to be solved for
        Pos: list of N_points positions in spherical coordinates (r, theta, phi)
        Acc: list of N_points accelerations in spherical coordinates (a_r, a_theta, a_phi)
    Output:
        Solved_Coef: line array of solved coefficients
# =============================================================================
# # ISSUES:
        - Diverges beyond lmax = 8
# =============================================================================
    """
    print(f"Solving for coefficients, with lmax = {lmax}")

    Acc_line = conv.Make_Line(Acc)

    M = Get_PotGradMatrix(lmax, Pos) # get M_PotGrad

   # Solved_coef = npl.solve(M.T.dot(M), M.T.dot(Acc_line.T)) #[1:]))
    Solved_coef = npl.lstsq(M, Acc_line)

    Acc_solved = M[:-1, :-1].dot(Solved_coef[:-1]) # change this "-1" to a "-3"?

    return Solved_coef, Acc_solved




# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_Gen_Matrix():
    file_name = "Polar_400km_EarthFixed_1jour_1sec.e"
    days = 0.0001
    Pos, Time = imp.Fetch_Pos(file_name, days)
    Mat = Get_PotGradMatrix(60,Pos)
    print(Mat)
    return Mat


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    np.set_printoptions(linewidth=np.inf, precision=2)

    Mat = TEST_Gen_Matrix()
    #print("\n\nnp.linalg.det(Mat.T @ Mat) =", np.linalg.det(Mat.T @ Mat))


    print("\nGH_solve done\n")