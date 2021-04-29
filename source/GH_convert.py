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

#import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


# =============================================================================
# FUNCTIONS FOR COORDINATES MANIPULATION
# =============================================================================
def thph2lola (theta, phi):
    """ converts radians to geographical long/lat (geodetic)
    """
    Lat = (pi/2 - theta) * 180/pi
    Long = (phi - pi) * 180/pi
    return Long, Lat
def lola2thph (Lat, Long):
    """ converts geographical long/lat to radians
    """
    theta = pi/2 - Lat*pi/180
    phi = Long*pi/180 + pi
    return theta, phi


def cart2sph (x,y,z):
    """ converts carthesian coordinates to spherical
    """
    radius    = np.sqrt(x**2 + y**2 + z**2)         # r
    elevation = np.arctan2(z, np.sqrt(x**2 + y**2)) # theta
    azimuth   = np.arctan2(y,x)                     # phi
    return radius, elevation, azimuth
def sph2cart222 (r,theta,phi): # cannot find where and when this function was used
    """ converts spherical coordinates to carthesian
    """
    x=r*cos(theta)*cos(phi)
    y=r*cos(theta)*sin(phi)
    z=r*sin(theta)
    return x, y, z
def sph2cart (r,theta,phi):
    """ converts spherical coordinates to carthesian, ISO convention
    """
    x=r*sin(theta)*sin(phi)
    y=r*sin(theta)*cos(phi)
    z=r*cos(theta)
    return x, y, z


def cart2sphA (pts):
    """ converts an array of carthesian coordinates to spherical
    """
    Pos = np.array([cart2sph(x,y,z) for x,y,z in pts])
    return Pos


def sph2cart_Grid(G_Grid, G_Long, G_Lat):
    """ returns the grid in a 3D x,y,z plottable format
    """
    X, Y, Z = sph2cart(G_Grid, pi/2-G_Lat*pi/180, G_Long*pi/180+pi)
    return X, Y, Z


def geodes2geocen (Lat_gd):
    """converts geodetic (or geographic) latitude into geocentric latitude
    """
    a = 6378137 # m
    f = 1/298.257223563  # Some constants
    b = a * (1-f) # m
#    Lat_gc = np.arctan( (b/a)**2 * np.tan(Lat_gd))
    Lat_gc = np.arctan(np.tan(Lat_gd) * (1-f)**2)
    return Lat_gc
def geocen2geodes (Lat_gd):
    """converts geocentric latitude into geodetic (or geographic) latitude
    """
    f = 1/298.257223563  # ellipsoid flattening
    Lat_gc = np.arctan(np.tan(Lat_gd) / (1-f)**2)
    return Lat_gc



# =============================================================================
# FUNCTIONS FOR ARRAY MANAGEMENT
# =============================================================================
def Make_Line (arr):
    """ Returns the array with all rows appended
    """
    return np.reshape(arr, (1, arr.size))

def Make_Line_acc(acc):
    L = len(acc[:,0])
    line = np.zeros(3*L)
    j = 0
    for i in acc:
        line[3*j: 3*(j+1)] = i
        j +=1
    return line


def Make_Array (line, col = 3):
    """ Returns the line list in an array of col columns
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
            There are no sine coeffs for degree m=0
            There are no coeffs for order l=0, l=1
    Output:
        HC_coef: solved spherical harmonic cosine coefficients
        HS_coef: solved spherical harmonic sine coefficients
            To fetch use: HS_coef(l,m) = {SIN_lm_coef}
    """
#    Cos_len = int( (lmax+1)*(lmax+2) /2 ) # c00,c10,c11,c20,c21,c22, ...
#    Sin_len = int( (lmax  )*(lmax+1) /2 ) # s11,s21,s22,s31,s32,s33, ...

    HC_coef = np.zeros( (lmax+1,lmax+1) )
    HS_coef = np.zeros( (lmax+1,lmax+1) )

    j = 0
    for l in range (2, lmax+1):
        for m in range (0, l+1):
            HC_coef[l, m] = CS[j] # Get the Cosine coefs out first
            j += 1
    # Normally, at this point, j == Cos_len
    for l in range (2, lmax+1):
        for m in range (1, l+1):
            HS_coef[l, m] = CS[j] # Get the Sine coefs out next
            j += 1

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
    Cos_len = int( (lmax+1)*(lmax+2) /2 ) -3 # c20,c21,c22, ...
    Sin_len = int( (lmax  )*(lmax+1) /2 ) -1 # s21,s22,s31,s32,s33, ...

    N_coef = Cos_len + Sin_len
    CS = np.zeros(N_coef)

    j=0
    for l in range (2, lmax +1):
        for m in range (0, l +1):
            CS[j] = HC[l, m] # Write in the Cosine coefs
            j += 1
    # Normally, at this point, j == Cos_len
    for l in range (2, lmax +1):
        for m in range (1, l +1):
            CS[j] = HS[l, m] # Write in the Sine coefs
            j += 1

    return CS



# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def TEST_Line_Array ():
    """ Tests the functions in this script.
    """
    A = np.asarray(np.linspace(0, 19,20))
    B = Make_Array(A, 4)
    print("B = \n", B, "\n")
    C = Make_Line(B)
    print("C = \n", C, "\n")

    CS = np.array([20,21,22,
                   30,31,32,33,
                   40,41,42,43,44,
                   50,51,52,53,54,55,# Cos coeffs
                   210,220,
                   310,320,330,
                   410,420,430,440,
                   510,520,530,540,550]) # Sin coeffs

    print("CS shape =", CS.shape)

    HC, HS = Make_Array_Coef(5, CS)
    print("HC = \n", HC, "\n")
    print("HS = \n", HS, "\n")

    lmax = 5
    Cos_len = int( (lmax+1)*(lmax+2) /2 ) -3 # c00,c10,c11,c20,c21,c22, ...
    Sin_len = int( (lmax  )*(lmax+1) /2 ) -1 # s11,s21,s22,s31,s32,s33, ...
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