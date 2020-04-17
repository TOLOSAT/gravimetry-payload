"""

@authors:

# =============================================================================
 Information:

    The functions in this script are used to import and fetch values and arrays
    from text files and from the system

# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from time import gmtime, strftime

#import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_basemap      as bmp
#import GH_harmonics    as harm
#import GH_geoMath      as gmath


from GH_convert import cart2sphA

# =============================================================================
# GLOBAL VARIABLES
# =============================================================================
data_path = "../data"



# =============================================================================
# FUNCTIONS - SYSTEM
# =============================================================================
def Get_Time (format_="%Y%m%d_%H%M%S"):
    """Returns the time in the given string format"""
    time = strftime(format_, gmtime())
    return time



# =============================================================================
# FUNCTIONS TO FETCH FILES
# =============================================================================
def Fetch_Pos (file_name, days=0.7, data_path="../data"):
    """
    Imports coordinates from file_name text file (generated from GMAT)

    Input:
        file_name: well, the file's name! remove all header text
        days: what time duration the outplut file should correspond to
              regardless of the sampling rate
        data_path: path to go and fetch the file
    Output:
        Pos: The position of the satellite in spherical coordinates
        Time: Associated time sampling of each position

    """
    Eph = np.loadtxt(f"{data_path}/{file_name}")
    t = np.array(Eph[:,0]) #time in seconds
    x = np.array(Eph[:,1]) #  \
    y = np.array(Eph[:,2]) #  | cordinates, in km
    z = np.array(Eph[:,3]) # /
    dt = np.int(t[1]*100)/100
    L = np.int(days*(86400/dt))
    # convert coord system and shorten array if needed
    pts = np.transpose(np.array([x,y,z]))
    if L >= len(pts):
        L = len(pts) # this is not necessary in python
    Pos = cart2sphA(pts[:L])
    Time = t[:L]
    print(f"Importing Pos file with {L} coordinates.")
    return Pos, Time


def Fetch_Coef ():
    """
    Returns the spherical harmonic coefficients for Earth's Geopotential
    Data originally extracted from : EGM2008_to2190_ZeroTide.txt
    These coef are already normalized
    These files exist with a degree up to lmax = 2190
    """
    data_path = "../data"
    HC = np.loadtxt(f"{data_path}/GeoPot_Coef_cos_deg30.txt")
    HS = np.loadtxt(f"{data_path}/GeoPot_Coef_sin_deg30.txt")
    return HC, HS


def Fetch_Topo_Coef ():
    """
    Returns the spherical harmonic coefficients for Earth's Topography
    Data originally extracted from : Coeff_Height_and_Depth_to2190_DTM2006.txt
    These coef are already normalized
    These files exist with a degree up to lmax = 2190
    """
    data_path = "../data"
    HC_topo = np.loadtxt(f"{data_path}/Height_Coef_cos_deg49.txt")
    HS_topo = np.loadtxt(f"{data_path}/Height_Coef_sin_deg49.txt")
    return HC_topo, HS_topo


# =============================================================================
# TEST FUNCTIONS
# =============================================================================



# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    HC, HS = Fetch_Coef()
    HC_topo, HS_topo = Fetch_Topo_Coef ()

    
    print("\nGH_import done")

