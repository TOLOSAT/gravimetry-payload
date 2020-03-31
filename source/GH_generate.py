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

import GH_convert     as conv
#import GH_import      as imp
#import GH_generate    as gen
import GH_solve       as solv
#import GH_displayCoef as dcoef
#import GH_displaySat  as dsat

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
    print("\nGenerating simulated acclerations, lmax =", lmax, "")
    
    CS = conv.Make_Line_Coef(lmax, HC, HS)
    print(f"Shape of the Coef array = {CS.shape}")
    
    M_PotGrad = solv.Get_PotGradMatrix(lmax, Pos) # get M_PotGrad    
#    print("shape of M=", M_PotGrad.shape)
    
    Acc_line = M_PotGrad.dot(CS)
#    print("shape of Acc_line=", Acc_line.shape)
    
    Acc_sim = conv.Make_Array(Acc_line, 3)
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

