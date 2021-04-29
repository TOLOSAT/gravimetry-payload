# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 20:21:55 2019

@author: Xavier de Labriolle, Antoine Bendimerad

Last edit : 27/02/2020

==============================================================================
 Information : 
     
     This python script uses the spherical coordinates system : 
         r = radius ; [0, np.inf()]
         theta = inclination from the z axis ; [0, pi]
         phi = rotation around the z axis ; [0, 2*pi]
     
         They must be adapted to lat/lon coordinates by adding : 
             -pi/2 to theta
             -pi to phi
     
         Their equivalent in the cartesian coordinates is : sph2cart(r,theta,phi):
             x=r*np.sin(theta)*np.cos(phi)
             y=r*np.sin(theta)*np.sin(phi)
             z=r*np.cos(theta)

    The main source for the mathematics involved in this code is : 
        "Definition of Functionals of the Geopotential 
        and Their Calculation from Spherical Harmonic Models"
        by Franz Barthelmes
        for ICGEM
        Revision: jan. 2013
        Will be refered to as "GFZ" from now on.
    
    The purpose of this script is to solve for spherical harmonic coefficients
        These coefficients must be stored in 2 arrays: 
            HC_coef: solved spherical harmonic cosine coefficients
            HS_coef: solved spherical harmonic sine coefficients
        where: 
            HS_coef(l,m) = SIN_lm_coef
            
==============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
import numpy.linalg as npl
from numpy import pi, sin, cos
from scipy.special import lpmn


# =============================================================================
# VARIABLES
# =============================================================================
GM = 3986004.415*10**8 # m**3 s**-2  
# wiki says : gm = 6.673*10**-11*5.975*10**24 = 398711749999999.94
a_ = 6378136.3 # m
A = 180/pi # degrees / rad   -------------- bad practice
Radius_eq = 6378137 # m
Radius_pol = 6356752.3142 # m




# =============================================================================
# FUNCTIONS FOR ARRAY MANAGEMENT
# =============================================================================
def Make_line (arr):
    """
    Returns the array with all rows appended
    """
    col = len(arr)
    length = col * len(arr[0]) # size(arr)
    line = np.reshape(arr, (1, length))
    return line


def Make_array (line, col = 3):
    """
    Returns the line list in an array of col columns
    """
    length = int(len(line)/col)
    arr = np.reshape(line, (col, length))
    return arr


def Make_array_coef (lmax, CS):
    """
    Returns the arrays of the solved Cosine and Sine coefficients
    
    Input: 
        CS: line array filled in coefficients in such manner :
        CS = [c00,s00,c10,s10,c11,s11,c20,s20,c21,s21,c22,s22,c30,s30, ... ]
                0   1   2   3   4   5   6   7   8   9  10  11  12  13  
    Output:
        HC_coef: solved spherical harmonic cosine coefficients
        HS_coef: solved spherical harmonic sine coefficients
        HS_coef(l,m) = SIN_lm_coef
        
    """
    # Exceptions
    if (lmax == int) or (lmax < 0):
        raise Exception("lmax should be an int and >= 0")
#    if ( != 1): #or (type (CS_solve) != numpy.ndarray):
#        raise Exception("CS is not what you think it is")
    
    HC_coef = np.zeros((lmax+1,lmax+1))
    HS_coef = np.zeros((lmax+1,lmax+1))
    #print("CS size=", len(CS))
    for l in range (0, lmax +1):
        for m in range (0, l +1):
            i = int(((l+1)*l/2 + m) * 2)
            #print("l m i=",l,m,i)
            HC_coef[l, m] = CS[i]
            HS_coef[l, m] = CS[i + 1]
            # end m loop
        # end l loop
    
    return HC_coef, HS_coef
    
def Make_line_coef (lmax, HC, HS):
    """
    Returns the line array filled of Cosine and Sine coefficients
    
    Input: 
        HC: spherical harmonic cosine coefficients
        HS: spherical harmonic sine coefficients
    Output:
        CS: line array filled in coefficients
        
    """   
    # Exceptions
    if (lmax == int) or (lmax < 0):
        raise Exception("lmax should be an int and >= 0")
    if (lmax+1 > len(HC)) or (lmax+1 > len(HS)):
        print("lmax is too large, not enough cofficients")
#    if (size(HC) != size(HS)):
#        print("HC and HS not the same size")

#    HS_sim = HS[:lmax+1, :lmax+1] # resize the coefs we want
#    HC_sim = HC[:lmax+1, :lmax+1] # DON'T NEED TO !!!
    
    N_coef = (lmax+2)*(lmax+1) # /2 *2 for 2 coefficients pet l.m.    
    CS = np.zeros(N_coef)
    
    j=0
    for l in range (0, lmax +1):
        for m in range (0, l +1):
            CS[j] = HC[l, m]
            CS[j + 1] = HS[l, m]
            
            j += 2
            # end m loop
        # end l loop
    
    return CS

"""
def Make_array_coef3 (lmax, CS):    
# =============================================================================
# does not need to be multiplied ! this function should not be used
# =============================================================================    
    # Exceptions
    if (lmax == int) or (lmax < 0):
        raise Exception("lmax should be an int and >= 0")    
    HC_coef = np.zeros((lmax+1,lmax+1))
    HS_coef = np.zeros((lmax+1,lmax+1))
    #print("CS size=", len(CS))
    for l in range (0, lmax +1):
        for m in range (0, l +1):
            i = int(((l+1)*l/2 + m) * 2)
            #print("l m i=",l,m,i)
            HC_coef[l, m] = (CS[3*i] + CS[3*i + 2] + CS[3*i + 4]) /3
            HS_coef[l, m] = (CS[3*i + 1] + CS[3*i + 3] + CS[3*i + 5]) /3
            # end m loop
        # end l loop    
    return HC_coef, HS_coef


def Make_line_coef3 (lmax, HC, HS): 
# =============================================================================
# does not need to be multiplied ! this function should not be used
# =============================================================================  
    # Exceptions
    if (lmax == int) or (lmax < 0):
        raise Exception("lmax should be an int and >= 0")
    if (lmax+1 > len(HC)) or (lmax+1 > len(HS)):
        print("lmax is too large, not enough cofficients")
#    if (size(HC) != size(HS)):
#        print("HC and HS not the same size")
#    HS_sim = HS[:lmax+1, :lmax+1] # resize the coefs we want
#    HC_sim = HC[:lmax+1, :lmax+1] # DON'T NEED TO !!!    
    N_coef = (lmax+2)*(lmax+1) # /2 *2 for 2 coefficients pet l.m.    
    CS = np.zeros(N_coef * 3)    
    j=0
    for l in range (0, lmax +1):
        for m in range (0, l +1):
            for i in range(0, 3): # need to have triple the coefficients
                CS[j] = HC[l, m]
                CS[j + 1] = HS[l, m]                
                j += 2
                # end i loop
            # end m loop
        # end l loop    
    return CS
"""


# =============================================================================
# FUNCTIONS FOR Sph Harm SOLVE
# =============================================================================
def Get_PotGradMatrix(lmax, Pos, R = 6378136.3) : 
    """
    Returns the matrix of the gravitational potential gradient. 
    Watch out, it gets big fast.
    Multiplying it with tha appropriate column vector of coefficients 
    will return the acceleration at the given coordinates.
    
    Input: 
        lmax: max order
        Pos: array of N_points positions in spherical coordinates (r, theta, phi)
        R: Reference radius in meters
    Output: 
        M_PotGrad: the matrix of the coefficients
        
    """
    #exceptions: only if a number is input.
    if (lmax == int) or (lmax < 0):
        raise Exception("lmax should be an int and >= 0")
#    if (len(Pos[0]) != 3) or (type (Pos) != numpy.ndarray):
#        raise Exception("Pos is not a list or has too many columns")

    # starting the function
    print("Get_PotGradMatrix")
    
    # constants
    GM = 3986004.415*10**8 # m**3 s**-2  
    # wiki says : gm = 6.673*10**-11*5.975*10**24 = 398711749999999.94
    
    N_points = len(Pos) # number of points
    N_coef = (lmax+2)*(lmax+1) # /2 *2 for 2 coefficients pet l.m.
       
    M_PotGrad = np.zeros((N_points * 3, N_coef)) #THE Potential Gradient Matrix
    
    for i in range (0, N_points):
        r, theta, phi = Pos[i] #spherical coordinates at the first point
        Plm_z, Plm_dz = lpmn(lmax, lmax, sin(phi))
        
        j = 0
        for l in range (0, lmax +1):
            for m in range (0, l +1):
                """
                These equations were found in the GFZ document page 23
                """
                W_r = - GM/r**2 * (R/r)**l * (l+1) * Plm_z[m, l]
                W_theta = W_r * m * r / (l+1)
                W_phi = GM/r * (R/r)**l * cos(phi)*Plm_dz[m, l]
                Sub_mat = np.zeros ((3,2))
                Sub_mat = [ [cos(m*theta)*W_r,      sin(m*theta)*W_r], 
                            [-sin(m*theta)*W_theta, cos(m*theta)*W_theta],
                            [cos(m*theta)*W_phi,    sin(m*theta)*W_phi] ]
                # multiply by: COS_lm_coef              SIN_lm_coef
                
                M_PotGrad [3*i : 3*(i+1), 2*j : 2*(j+1)] = Sub_mat
                j += 1
                
                # end m loop
            # end l loop
            
            if np.mod(i, 50) == 1 :
                print(i)
        # end i loop
        
    return M_PotGrad


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
# =============================================================================
# # ISSUES:
        - The npl.solve function wants square matrices
        - The sine coefficients for order m=0 are also 0. the matrix is singular 
          and cannot be inverted. 

# # SOLUTIONS: 
        - Put all the cosine before the sine, just like I did last time
        - Find another function that does the job
        
# =============================================================================
# =============================================================================

    """
    # Exceptions
#    if (lmax == int) or (lmax < 0):
#        raise Exception("lmax should be an int and >= 0")
#    if (len(Pos[0]) != 3) or (type(Pos) != np.ndarray):
#        raise Exception("Pos is not an array or has too many columns")
#    if (len(Acc[0]) != 3) or (type(Acc) != np.ndarray):
#        raise Exception("Acc is not an array or has too many columns")
    
    Acc_line = Make_line(Acc)
    
    M = Get_PotGradMatrix(lmax, Pos) # get M_PotGrad

#    Solved_coef = npl.solve(M.T.dot(M), M.T.dot(Acc_line.T)) #[1:]))
#    Solved_coef = npl.solve(M, Acc_line)
    
    
    Acc_solved = M[:-1, :-1].dot(Solved_coef[:-1])

    return Solved_coef, Acc_solved

    
    
    
# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def TEST ():
    """
    Tests the functions in this script.
    """
    A = np.asarray(np.linspace(0, 19,20)) #array [[0,1,2,3,4], [5,6,7,8,9], [10,11,12,13,14],[15,16,17,18,19]]   
    B = Make_array(A, 4)
    print("B = \n", B)    
    C = Make_line(B)
    print("C = \n", C)    
    
    CS = np.array([00,00,
                   10,10,11,11,
                   20,20,21,21,22,22,
                   30,30,31,31,32,32,33,33,
                   40,40,41,41,42,42,43,43,44,44])
    CS3 = np.array([00,00,00,00,00,00,
                   10,10,11,11,10,10,11,11,10,10,11,11,
                   20,20,21,21,22,22,20,20,21,21,22,22,20,20,21,21,22,22,
                   30,30,31,31,32,32,33,33,30,30,31,31,32,32,33,33,30,30,31,31,32,32,33,33])
    print("CS =", len(CS))
    
    HC, HS = Make_array_coef(3, CS)
    print("HS = \n", HS)
    
    CS2 = Make_line_coef(3, HC, HS)
    print("CS = \n",CS2)
    
    # end TEST function

# =============================================================================
# MAIN 
# =============================================================================

#TEST()

"""
Process : 
    1.  Generate the acceleration values (from orbit data or raw simulation)
    2.  Solve for the coefficients using the functions in this script
    3.  Generate a Geoid map from solved coefficients
    4.  Compare Geoids or coefficients
"""


print("\n GH_solve done")
    

