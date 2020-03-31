"""

@authors:

# =============================================================================
 Information: 
    
    The functions in this script are used to convert variables in their nature, 
    or in their shapes.
        
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from numpy import pi, sin, cos

#import GH_convert     as conv
#import GH_import      as imp
#import GH_generate    as gen
#import GH_solve       as solv
#import GH_displayCoef as dcoef
#import GH_displaySat  as dsat

# =============================================================================
# FUNCTIONS FOR COORDINATES MANIPULATION
# =============================================================================
def cart2sph(x,y,z):
    """
    converts carthesian coordinates to spherical
    """
    X2_Y2 = x**2 + y**2
    r = np.sqrt(X2_Y2 + z**2)              # r
    elev = np.arctan2(z, np.sqrt(X2_Y2))   # theta
    az = np.arctan2(y,x)                        # phi
    return r, elev, az

def cart2sphA(pts):
    """
    converts an array of carthesian coordinates to spherical
    """
    Pos = np.array([cart2sph(x,y,z) for x,y,z in pts])
    return Pos

def sph2cart(r,theta,phi):
    """
    converts spherical coordinates to carthesian
    """    
    x=r*cos(theta)*cos(phi)
    y=r*cos(theta)*sin(phi)
    z=r*sin(theta)
    return x, y, z

    
# =============================================================================
# FUNCTIONS FOR ARRAY MANAGEMENT
# =============================================================================
def Make_Line (arr):
    """
    Returns the array with all rows appended
    """
    col = len(arr)
    length = col * len(arr[0]) # size(arr)
    line = np.reshape(arr, (1, length))
    return line


def Make_Array (line, col = 3):
    """
    Returns the line list in an array of col columns
    """
    length = int(len(line)/col)
    arr = np.reshape(line, (length, col))
    return arr


def Make_Array_Coef (lmax, CS):
    """
    Returns the arrays of the solved Cosine and Sine coefficients
    
    Input: 
        CS: line array filled in coefficients in such manner :
        CS = [c00,c10,c11,c20,c21,c22, ... s11,s21,s22,s31,s32,s33 ... ]
            There is no sine coef for order l degree m=0
    Output:
        HC_coef: solved spherical harmonic cosine coefficients
        HS_coef: solved spherical harmonic sine coefficients
        HS_coef(l,m) = SIN_lm_coef
  
    """
#    Cos_len = int( (lmax+1)*(lmax+2) /2 ) # c00,c10,c11,c20,c21,c22, ... 
#    Sin_len = int( (lmax  )*(lmax+1) /2 ) # s11,s21,s22,s31,s32,s33, ...
    
    HC_coef = np.zeros( (lmax+1,lmax+1) ) 
    HS_coef = np.zeros( (lmax+1,lmax+1) ) 
    
    j = 0
    for l in range (0, lmax +1):
        for m in range (0, l +1):
            HC_coef[l, m] = CS[j] # Get the Cosine coefs out first
            j += 1
    # end loops
    # Normally, at this point, j == Cos_len
    
    for l in range (1, lmax +1):
        for m in range (1, l +1):  
            HS_coef[l, m] = CS[j] # Get the Sine coefs out next
            j += 1
    # end loops    
    
    return HC_coef, HS_coef


def Make_Line_Coef (lmax, HC, HS):
    """
    Returns the line array filled of Cosine and Sine coefficients
    
    Input: 
        HC: spherical harmonic cosine coefficients
        HS: spherical harmonic sine coefficients
    Output:
        CS: line array filled in coefficients
        
    """   
    Cos_len = int( (lmax+1)*(lmax+2) /2 ) # c00,c10,c11,c20,c21,c22, ... 
    Sin_len = int( (lmax  )*(lmax+1) /2 ) # s11,s21,s22,s31,s32,s33, ...
    
    N_coef = Cos_len + Sin_len
    CS = np.zeros(N_coef)
    
    j=0
    for l in range (0, lmax +1):
        for m in range (0, l +1):
            CS[j] = HC[l, m] # Write in the Cosine coefs 
            j += 1
    # end  loops
    # Normally, at this point, j == Cos_len 

    for l in range (1, lmax +1):
        for m in range (1, l +1):
            CS[j] = HS[l, m] # Write in the Sine coefs            
            j += 1
    # end loops    
    
    return CS
    
    
    
# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def TEST_Line_Array ():
    """
    Tests the functions in this script.
    """
    A = np.asarray(np.linspace(0, 19,20)) 
    B = Make_Array(A, 4)
    print("B = \n", B, "\n")    
    C = Make_Line(B)
    print("C = \n", C, "\n")    
    
    CS = np.array([5,
                   10,11,
                   20,21,22,
                   30,31,32,33,
                   40,41,42,43,44, # Cos coeffs
                   110,
                   210,220,
                   310,320,330,
                   410,420,430,440]) # Sin coeffs
                    
    print("CS shape =", CS.shape)
    
    HC, HS = Make_Array_Coef(4, CS)
    print("HC = \n", HC, "\n")
    print("HS = \n", HS, "\n")
    
    lmax = 4
    Cos_len = int( (lmax+1)*(lmax+2) /2 ) # c00,c10,c11,c20,c21,c22, ... 
    Sin_len = int( (lmax  )*(lmax+1) /2 ) # s11,s21,s22,s31,s32,s33, ...
    print("cos sin lengths =",Cos_len, ",",Sin_len)
    
    CS2 = Make_Line_Coef(lmax, HC, HS)
    print("CS =",CS2.shape, "\n", CS2)
    return CS2


# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    CS = TEST_Line_Array()
        
    print("\nGH_convert done")

