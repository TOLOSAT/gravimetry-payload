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
import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


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
def Fetch_Pos (file_name, days=0.7, data_path="../data", spherical = True ):
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
    if spherical :
        Pos = cart2sphA(pts[:L])
    else:
        Pos = pts[:L]
    Time = t[:L]
    return Pos, Time




def Fetch_Pos_Vit (file_name, days=0.7, data_path="../data", spherical = True ):
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
    vx = np.array(Eph[:,4])
    vy = np.array(Eph[:,5])
    vz = np.array(Eph[:,6])

    dt = np.int(t[1]*100)/100
    L = np.int(days*(86400/dt))
    # convert coord system and shorten array if needed
    pts = np.transpose(np.array([x,y,z]))

    ptsVit = np.transpose(np.array([vx,vy,vz]))

    if L >= len(pts):
        L = len(pts) # this is not necessary in python
    if spherical :
        Pos = cart2sphA(pts[:L])
        Vit = cart2sphA(ptsVit[:L])
    else:
        Pos = pts[:L]
        Vit = ptsVit[:L]
    Time = t[:L]
    return Pos,Vit, Time


def Fetch_Coef (data="subset"):
    """
    Returns the spherical harmonic coefficients for Earth's Geopotential
    Data originally extracted from : EGM2008_to2190_ZeroTide.txt
    These coef are already normalized
    A subset is returned unless full coefficient matrix is specified
    """
    data_path = "../data"
    if (data == "full"):
        HC = np.loadtxt(f"{data_path}/GeoPot_Coef_cos_deg2190.txt")
        HS = np.loadtxt(f"{data_path}/GeoPot_Coef_sin_deg2190.txt")
    else:
        HC = np.loadtxt(f"{data_path}/GeoPot_Coef_cos_deg30.txt")
        HS = np.loadtxt(f"{data_path}/GeoPot_Coef_sin_deg30.txt")
    return HC, HS


def Fetch_Topo_Coef (data="subset"):
    """
    Returns the spherical harmonic coefficients for Earth's Topography
    Data originally extracted from : Coeff_Height_and_Depth_to2190_DTM2006.txt
    These coef are already normalized
    A subset is returned unless full coefficient matrix is specified
    """
    data_path = "../data"
    if (data == "full"):
        HC_topo = np.loadtxt(f"{data_path}/Height_Coef_cos_deg2190.txt")
        HS_topo = np.loadtxt(f"{data_path}/Height_Coef_sin_deg2190.txt")
    else:
        HC_topo = np.loadtxt(f"{data_path}/GeoPot_Coef_cos_deg30.txt")
        HS_topo = np.loadtxt(f"{data_path}/GeoPot_Coef_sin_deg30.txt")
    return HC_topo, HS_topo


def Load_GLl (detail="zeros"):
    """
    Should be used with exp.Store_temp_GLl()
    You MUST move the files you're interested in, they are in Rendered
    Load a G_Grid, G_Long, G_Lat
    """
    GLl_path = "../Rendered/grid"
    G_Grid = np.loadtxt(f"{GLl_path}/{detail} G_Grid")
    G_Long = np.loadtxt(f"{GLl_path}/{detail} G_Long")
    G_Lat  = np.loadtxt(f"{GLl_path}/{detail} G_Lat")
    return G_Grid, G_Long, G_Lat


def Load_gridget_xmin(shape = (61,61), name="pyOUTPUT.txt"):
    """
    This function is to be used along with py_gridget_xmin.py
    It renders a temporary 3-column file with the grid data
    This function imports that file, and returns the grid and meshgrid
    """
#    F77_path = "../../Fortran NGA"
    F77_path = "../Rendered/temp"
    GLl_path = "../Rendered/grid"

    RAW = np.loadtxt(f"{F77_path}/{name}")
    print(f"Importing {name}: RAW length={len(RAW)}")
#    shape = (181,361) # if length  = 65160
#    shape = (181,360) # if length  = 65160
#    shape = (19,19)
#    if (len(RAW)==703): shape = (19,37)
#    if (len(RAW)==65341): shape = (181,361)
    G_Lat  = np.flip(np.reshape(RAW[:,0], shape), 0)
    G_Long = np.reshape(RAW[:,1], shape)
    G_Grid = np.flip(np.reshape(RAW[:,2], shape), 0)

    # resize them ?
    return G_Grid, G_Long, G_Lat



# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_load_temp():
#    exp.TEST_store_temp()
    A, B, C = Load_GLl()
    print(A); print(B); print(C)



# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
#    HC, HS = Fetch_Coef()

#    HC_topo, HS_topo = Fetch_Topo_Coef ("full")

#    TEST_load_temp()

    G_Grid, G_Long, G_Lat = Load_gridget_xmin()

#    exp.Store_temp_GLl(G_Grid, G_Long, G_Lat, "EGM2008 1h")

    print("\nGH_import done")