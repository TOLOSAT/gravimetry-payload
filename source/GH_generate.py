"""

@authors:

# =============================================================================
 Information: 
    
    The purpose of this script is to generate acceleration arrays. 
    Random generation, and position derivation
    
    These arrays must be in spherical coordinates, in km.
        
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
#import numpy as np
#from numpy import pi, sin, cos, tan

from GH_convert import *
from GH_import import *
from GH_solve import *

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
    print("\nGenerating simulated acclerations, lmax =", lmax, "\n")
    
    CS = Make_Line_Coef(lmax, HC, HS)
    print("shape of CS =", CS.shape)
    
    M_PotGrad = Get_PotGradMatrix(lmax, Pos) # get M_PotGrad    
#    print("shape of M=", M_PotGrad.shape)
    
    Acc_line = M_PotGrad.dot(CS)
#    print("shape of Acc_line=", Acc_line.shape)
    
    Acc_sim = Make_Array(Acc_line, 3)
#    print("shape of Acc_sim=", Acc_sim.shape)
    
    return Acc_sim


    
# =============================================================================
# TEST FUNCTIONS
# =============================================================================
    
# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    
    print("\nGH_generate done")

