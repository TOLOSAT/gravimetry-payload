"""
@authors:
# =============================================================================
 Information:
    The purpose of this script is to calculate the sums from spherical
    harmonic coefficients, and generate grids of data mapped out over the
    surface of the Earth
    Generally used variables:
        R, Lat, Long  = coordinates in the geographic (geodesic) CRS
        R, theta, phi = coordinates in the geocentric CRS - ISO convention
        	R     <=> R    = radius in km ; [0, inf()]
        	theta <=> LAT  = inclination around the z axis in radians; [0, pi]
        	phi   <=> LONG = azimuth around the z axis in radians; [0, 2*pi]
        a_e, R_e, GM_e = relative to the reference ellipsoid
        a_g, R_g, GM_g = relative to the geopotential model
        lmax   = maximum degree to calculate to for the geopotantial
        HC, HS = Geopotential stokes coefficients
        lmax_topo, HC_topo, HS_topo = same, but for topography.
        limits = [Western_long, Eastern_long, Southern_lat, Northern_lat]
        mins = The grid resolution in arc minutes
debug:
    the Sph Harm canot be computed for degrees above 154.
    for l = 155 and beyond, it does not work
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin, cos
from time import sleep

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
import GH_earthMap     as emap



# =============================================================================
# FUNCTIONS TO GENERATE DATA ARRAYs
# =============================================================================
def init_grid (mins=0, limits=np.array([-180, 180, -90, 90])):
    """
    Initiates the grid variables based on the number of points wanted
    within the given limits
    gris step is in minutes, 60 min = 1 degree
    """
    dim = limits * pi/180

    if ( (mins < 1) or (mins > 600) ):
        size_long = 5
        size_lat  = 5
        Line_theta = np.linspace(dim[0], dim[1], size_long)
        Line_phi   = np.linspace(dim[2], dim[3], size_lat)
    else:
        step = mins/60 * pi/180
        Line_theta = np.arange(dim[0], dim[1]+step/2, step)
        Line_phi   = np.arange(dim[2], dim[3]+step/2, step)
        size_long = len(Line_theta)
        size_lat  = len(Line_phi)

    G_theta, G_phi = np.meshgrid(Line_theta, Line_phi)
    G_Grid = np.zeros((size_lat, size_long))

    return G_Grid, G_theta, G_phi


def Gen_Grid (mins, Get_FUNCTION, in_args, limits=np.array([-180, 180, -90, 90])):
    """
    This function generates a grid of the desired spherical harmonic model
    at Lat/Long coordinates
    Input:
        mins: hthe grid resolution in arc minutes
        Get_FUNCTION: the callable function that must be used
        *in_args: the arguments to the callable function besides R, phi, theta
        limits: the geographical limits to the Long/lat map
    Output:
        G_Grid: grid of Get_FUNCTION(R,phi,theta,*in_args)
        G_Long: grid of longitudes, [mins] step, within bounraries [limits]
        G_Lat:  same for latitudes
    """
    G_Grid, G_theta, G_phi = init_grid(mins, limits)
    print(f"Making a grid with \"{Get_FUNCTION.__name__}()\", with {G_Grid.size} points\n",end="\r")

    it=0
    for j in range(0, G_phi.shape[0]):
        phi = pi/2 - G_phi[j][0]
        R_e = gmath.Get_Ellipsoid_Radius(phi)

        for i in range(0, G_theta.shape[1]):
            term.printProgressBar(it+1, G_phi.size); it+=1 # print(f"\rLong =  {theta*pi/180-180} ;Lat {90-phi*pi/180}",end="\r")
            theta = G_theta[0][i]+ pi
            G_Grid[j,i] = Get_FUNCTION(R_e, phi, theta, *in_args)

    return G_Grid, G_theta*180/pi, G_phi*180/pi # in degrees, L





# =============================================================================
# FUNCTIONS TO CALCULATE SPHERICAL HARMONIC SUMS
# =============================================================================
def Get_Topo_Height (R_e, phi, theta,   lmax_topo, HC_topo, HS_topo):
    """
    This function returns the height of Earth's estimated topography at phi/theta
    coordinates
    The solution is calculated up to degree lmax in the [HC_topo,HS_topo] model
    """
    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax_topo, lmax_topo, cos(phi)) # I am allowed to write that.
#    P_lm = gmath.ALF_norm_gcb(lmax_topo, lmax_topo, phi).T
    for l in range (0, lmax_topo+1):
        Sum2 = 0
        for m in range (0, l+1):
#            print(f"lm= {l} {m}")
#            print(f"\rl={l} ; m={m}",end="\r")
            Sum2 += (HC_topo[l,m]*cos(m*theta) + HS_topo[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
#            Sum2 += (HC_topo[l,m]*cos(m*theta) + HS_topo[l,m]*sin(m*theta)) * P_lm[m, l]
        Sum1 += Sum2

    return Sum1


def Get_Geo_Pot (R_e, phi, theta,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the potential at given height/phi/theta coordinates
    The solution is calculated up to degree lmax in the HC HS model
    """
    cst = gmath.Constants()

    R_t = R_e #+ Get_Topo_Height (R_e, phi, theta,    lmax_topo, HC_topo, HS_topo)
    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(phi))
#    LPNM = gmath.ALF_norm_gcb(lmax, lmax, phi)
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, 1): # l+1
            Sum2 += (HC[l,m]*cos(m*theta) + HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize1(l, m)
#            Sum2 += (HC[l,m]*cos(m*theta) + HS[l,m]*sin(m*theta)) * LPNM[l, m]
        Sum1 += (cst.a_g/R_t)**l * Sum2

    geopot = cst.GM_g/R_t*(1 + Sum1)

    return geopot


def Get_Geoid_Height (R_e, phi, theta,    lmax, HC, HS):
    """
    This function returns the potential at given height/phi/theta coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The cosine coefficients for even l and m=0 are corrected to remove the
    reference ellipsoid from the results
    Equations come from the geoid cook book
    """
    cst = gmath.Constants()
    g_0 = gmath.Get_Normal_Gravity(phi)
    phi_gc = conv.geodes2geocen(phi)

#    R = cst.a_g
#    g_0 = cst.g
    phi_gc = phi

    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(phi_gc) )
#    LPNM = gmath.ALF_norm_gcb(lmax, lmax, phi_gc)
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            HC_lm = CorrCos_lm(l, m, HC[l,m])
            Sum2 += (HC_lm*cos(m*theta) + HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
#            Sum2 += (HC[l,m]*cos(m*theta) + HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
#            Sum2 += (HC_lm*cos(m*theta) + HS[l,m]*sin(m*theta)) * LPNM[l, m]

        Sum1 +=  (cst.a_g/R_e)**l * Sum2

    Geo_H = cst.GM_g * Sum1 / (R_e*g_0)

    return Geo_H


'''
def Get_Geoid_Height2 (R_e, phi, theta,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the potential at given height/phi/theta coordinates
    The solution is calculated up to degree lmax in the HC HS model
    Equations from the GFZ handbook, eq.116
    """
    c = gmath.Constants()
    a_g=c.a_g; GM_g=c.GM_g; # g_e=c.g
    G=c.G; ro=c.ro;
    g_e = gmath.Get_Normal_Gravity(phi)
#    phi_gc = conv.geodes2geocen(phi)
    phi_gc = phi
    Sum_geo = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, sin(phi_gc) )
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*theta)+HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum_geo +=  (a_g/R_e)**l * Sum2
    Sum_topo = 0
    P_lm, _ = gmath.Pol_Legendre(lmax_topo, lmax_topo, cos(phi))
    for l in range (0, lmax_topo+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC_topo[l,m]*cos(m*theta) + HS_topo[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum_topo += Sum2
    Geo_H = GM_g/(R_e*g_e) * Sum_geo  -  2*pi*G*ro/g_e * (R_e*Sum_topo)**2
    return Geo_H
'''


def Get_acceleration (R_e, phi, theta,    lmax, HC, HS):
    """
    This function returns the acceleration at given height/phi/theta coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The acceleration is the gratient of the geopotential, and is calculated
    over a distance "d"
    """
    c = gmath.Constants()
    a_g=c.a_g; GM_g=c.GM_g;

    d = 1 # m
    R_1 = R_e - d/2
    R_2 = R_e + d/2

    Sum1 = 0
    Sum2 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(phi))
    for l in range (2, lmax+1):
        Sum3 = 0
        for m in range (0, l+1):
            Sum3 += (HC[l,m]*cos(m*theta) + HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum1 += (a_g/R_1)**l * Sum3
        Sum2 += (a_g/R_2)**l * Sum3

    GP_1 = GM_g/R_1*(1 + Sum1)
    GP_2 = GM_g/R_2*(1 + Sum2)
    local_acc = (GP_1 - GP_2) / d
    return local_acc


def Get_acceleration2 (R_e, phi, theta,    lmax, HC, HS):
    """
    This function returns the acceleration at given height/phi/theta coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The equations come from my own calculation of the derivative of the
    potential along the r axis, with respect to r.
    """
    c = gmath.Constants()
    a_g=c.a_g; GM_g=c.GM_g;

    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(phi))
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*theta) + HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum1 += (a_g/R_e)**l * -l/R_e * Sum2

    W_ar = -GM_g/R_e**2 * Sum1
    return W_ar


def Get_acceleration3 (R_e, phi, theta,    lmax, HC, HS):
    """
    This function returns the acceleration at given height/phi/theta coordinates
    The solution is calculated up to degree lmax in the HC HS model
    The equations come from the GFZ handbook :
        the derivative of the potential along the r axis
    """
    c = gmath.Constants()
    a_g=c.a_g; GM_g=c.GM_g;

    Sum1 = 0
    P_lm, _ = gmath.Pol_Legendre(lmax, lmax, cos(phi))
    for l in range (2, lmax+1):
        Sum2 = 0
        for m in range (0, l+1):
            Sum2 += (HC[l,m]*cos(m*theta) + HS[l,m]*sin(m*theta)) * P_lm[m, l] * gmath.Normalize(l, m)
        Sum1 += (a_g/R_e)**l * (l+1) * Sum2

    W_ar = GM_g/R_e**2 *Sum1
    return W_ar


def Get_isopot (R_e, phi, theta,    W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    This function returns the Height at given phi/theta coordinates at which the
    geopotential is equal to the given W_0
    The solution is calculated up to degree lmax in the HC HS model
    The approach is a dichotomy within a width of "de" m error
    """
    de = 10
    grad = -10 # 9.81

    arg_after = [phi, theta, lmax, HC, HS, lmax_topo, HC_topo, HS_topo]
    R_iso = gmath.dichotomy_grad (Get_Geo_Pot, [], R_e, arg_after, W_0, de, grad)
    Height = R_iso - R_e
    return Height # , R_iso



# =============================================================================
# SUB FUNCTIONS BUT STILL HARMONICS
# =============================================================================
def CorrCos_lm (l, m, HC_lm):
    """ corrects the HC_lm coef if the conditions are correct - GCB """
    lmax_corr=20
    if ( (m==0) and (l<=lmax_corr) and (l%2 == 0) ): # print(f"corr lm {l} {m}")
        HC_lm -= Cosine_Correction2(l)
    return HC_lm
def Cosine_Correction2 (N):
    """ raw data from the hsynth fortran code """
    C = np.zeros(21)
    C[ 2] =  -0.484166774985E-03
    C[ 4] =   0.790303733524E-06
    C[ 6] =  -0.168724961164E-08
    C[ 8] =   0.346052468485E-11
    C[10] =  -0.265002226377E-14
    C[12] =  -0.410790140988E-16
    C[14] =   0.447177356743E-18
    C[16] =  -0.346362564558E-20
    C[18] =   0.241145603096E-22
    C[20] =  -0.160243292770E-24
    return C[N]


def Cosine_Correction (N):
    """
    Returns the corrected Cosine coefficient
    The maths in this book comes from the Geoid cook book and the fortran code
    I dissasembled
    """
    n = N/2
    c=gmath.Constants()
    GM_e = c.GM_e;
    GM_g = c.GM_e
    a_e = c.a_g
    a_g = c.a_g
    E = c.E
    m = c.m
    e_2 = c.e_2

    q0 = 1/2*((1 + 3/e_2**2) * np.arctan(e_2 - 3/e_2))
    CA_ME2 = 1/3*(1 - 2/15*m*e_2/q0)

    J_2n = (-1)**(n+1) * 3*(E/a_e)**(2*n) / ((2*n+1)*(2*n+3))
    J_2n = J_2n * (1 - n + 5*n*CA_ME2)

    Corr_2n = J_2n / np.sqrt(4*n + 1) * GM_e/GM_g * (a_e/a_g)**n

    return Corr_2n



# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_gen_grid():
    mins = 600
    limits = np.array([-180, 180, -90, 90])

#    G_Grid, G_theta, G_phi = init_grid(mins, limits)
    G_Grid, G_theta, G_phi =  Gen_Grid(mins, void, [], limits)

    FIG = plt.figure()
    AX = FIG.add_subplot("111")
    data = AX.contourf(G_theta, G_phi, G_Grid)
    _ = plt.colorbar(mappable=data, ax=AX)
    return G_theta, G_phi, G_Grid
def void(r, phi, theta):
    return 0


def Math_calc_geopot_basic(z):
    """ some function needed in TEST_plot_radius """
    G = 6.673E-11
    M = 5.975E24
    a = 6.378E6
    P = G*M*(1/a + 1/(a+float(z)))
    return P

def TEST_plot_radius(fignum, lmax, Lat, Long, HC, HS, lmax_topo, HC_topo, HS_topo):
    """
    Plots the geopotential at given coordinates from Earth's center to the
    surface, and beyond into space.
    """
    R_earth = gmath.Get_Ellipsoid_Radius(Lat*pi/180)
    Rs = np.linspace(95, 105, 100)
    G_Pot = np.zeros(Rs.shape)
    G_Pot_Basic = np.zeros(Rs.shape)

    for i in range (len(Rs)):
        G_Pot[i]       = Get_Geo_Pot (Rs[i]*R_earth, Lat*pi/180, Long*pi/180, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)
        G_Pot_Basic[i] = Get_Geo_Pot (Rs[i]*R_earth, Lat*pi/180, Long*pi/180, 2, HC, HS, lmax_topo, HC_topo, HS_topo)
#        G_Pot_Basic[i] = Math_calc_geopot_basic(Rs[i]*R_earth)

    plt.figure(fignum)
    plt.plot(Rs*R_earth, G_Pot, label=f"{Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot - G_Pot_Basic, label=f"{Lat}-{Long}; {lmax}")
#    plt.plot(Rs, G_Pot_Basic, label=f"basic {Lat}-{Long}; {lmax}")
    plt.title("geopotential against the radius of the Earth, loop lmax")
    plt.xlabel("Distance from the center to the surface of the Earth (in %)")
    plt.ylabel("local value of the geopotential (m^2/s^2)")
    plt.legend(fontsize = 8)
    plt.show(block=False)


def TEST_plotGeoPot_radius():
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    plt.figure()
    plt.clf()
    for i in range (2, 15):
        TEST_plot_radius(1, i*2, 50, 50, HC, HS, 10, HC_topo, HS_topo)


def TEST_lmax_loop_long_line():
    """ plots the geoid height at the equator, around the world """
    plt.figure()
    plt.clf()
    plt.grid(True)

    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
#    lmax_topo = 10

    lmaxs = np.arange(3, 25, 2)
    for lmax in lmaxs:
        Lat = pi/180 * 40
        R = gmath.Get_Ellipsoid_Radius(Lat)
        Longs = np.linspace(0, 2*pi, 91)

        Geo_H = np.zeros(len(Longs))

        for i in range(len(Longs)):
            Long = Longs[i]
#            Geo_H[i] = Get_acceleration   (R, Lat, Long,    lmax, HC, HS); title_spec="Acceleration"
#            Geo_H[i] = Get_Topo_Height   (R, Lat, Long,    lmax_topo, HC_topo, HS_topo); title_spec="Topography height"
#            Geo_H[i] = Get_Geo_Pot       (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="GeoPot"
            Geo_H[i] = Get_Geoid_Height  (R, Lat, Long,    lmax, HC, HS); title_spec="Geoid height"
#            Geo_H[i] = Get_Geoid_Height2 (R, Lat, Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="Geoid height"

        Longs = (Longs-pi) * 180/pi
        plt.plot(Longs, Geo_H, label=f"lx={lmax}")

    plt.suptitle(f"{title_spec} at equator (m) vs Longitude; loop lmax")
    plt.legend()


def TEST_lmax_loop_lat_line():
    """ plots the geoid height at the equator, around the world """
    plt.figure()
    plt.clf()
    plt.grid(True)

    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
#    lmax_topo = 10

#    lmaxs = [155]
#    lmaxs = np.arange(5, 180, 10)
#    lmaxs = np.array([5, 15, 35, 60, 150, 600])
    lmaxs = np.arange(1, 25, 3)
    for lmax in lmaxs:
        Long =  80
        Lats = np.linspace(0, pi, 91)
        print(f"making lmax = {lmax}")
        Geo_H = np.zeros(len(Lats))

        for i in range(len(Lats)):
            Lat = Lats[i]
#            print(f"\tmaking lat = {(90-Lat*180/pi):0.2f}")
            R = gmath.Get_Ellipsoid_Radius(Lat)

            Geo_H[i] = Get_acceleration2 (R, Lat, pi/180 *Long,    lmax, HC, HS); title_spec="Acceleration"
#            Geo_H[i] = Get_Topo_Height   (R, Lat, pi/180 *Long,    lmax, HC_topo, HS_topo); title_spec="Topography height"
#            Geo_H[i] = Get_Geo_Pot       (R, Lat, pi/180 *Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="GeoPot"
#            Geo_H[i] = Get_Geoid_Height  (R, Lat, pi/180 *Long,    lmax, HC, HS); title_spec="Geoid height"
#            Geo_H[i] = Get_Geoid_Height2 (R, Lat, pi/180 *Long,    lmax, HC, HS, lmax_topo, HC_topo, HS_topo); title_spec="Geoid height"

        Lats = (pi/2-Lats) * 180/pi
        plt.plot(Lats, Geo_H, label=f"lx={lmax}")

    plt.suptitle(f"{title_spec} at long{Long} vs Latitude; loop lmax")
    plt.legend()
    return Geo_H


def Get_isopot_average ():
    """ Returns the geopotential average at the surface of the ellipsoid """
#    HC, HS = imp.Fetch_Coef()
#    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
#
#    lmax_av = 29; lmax_topo_av = 48; mins_av = 10
#    Grid, _, _ = Gen_Grid (mins_av, Get_Geo_Pot, [lmax_av, HC, HS, lmax_topo_av, HC_topo, HS_topo])
#    mm = np.amin(Grid)
#    MM = np.amax(Grid)
#    W_0 = np.average(Grid)
#    print(f"lmax_av = {lmax_av}; lmax_topo_av={lmax_topo_av}; mins={mins_av}; \n\tmm={mm}; \n\tMM={MM}; \n\tW_0={W_0}")
#
#    mm=62499132.77072437
#    MM=62680791.364166744
    W_0=62601662.83663869
    return W_0


def TEST_Get_isopot ():
    HC, HS = imp.Fetch_Coef()
    HC_topo, HS_topo = imp.Fetch_Topo_Coef()
    Lat =  11
    Long =  10
    lmax = 5
    lmax_topo = 10

    W_0 = Get_isopot_average()
    R_e = gmath.Get_Ellipsoid_Radius(Lat)
    R_iso, height = Get_isopot(R_e, pi/180*Lat, pi/180*Long, W_0, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)

    print(f"Average potential at {Lat} {Long} is at R={R_iso}, H={height}")


def TEST_Cosine_corr():
    HC, HS = imp.Fetch_Coef()
    lmax = 10
    CCos_2n0 = np.zeros(lmax)

    print("Ref. Gravity Potential Even Zonal Terms (C-form)")
    print("Subtracted From the Input Coefficients:")

    for l in range(1,lmax+1):
        CCos_2n0[l-1] = Cosine_Correction(2*l, HC[2*l, 0])
        print(f"C({2*l},0) =\t{CCos_2n0[l-1]}")


def TEST_high_lmax():
    HC, HS = imp.Fetch_Coef("full")
    HC_topo, HS_topo = imp.Fetch_Topo_Coef("full")
    lmax =154
    Lat = 60
    Long = 60
#    lmax_topo = 10

    theta, phi = conv.lola2thph(Lat,Long)
    R = gmath.Get_Ellipsoid_Radius(Lat)

    val = Get_Topo_Height   (R, phi, theta, lmax, HC_topo, HS_topo)
#    val = Get_Geo_Pot       (R, phi, theta, lmax, HC, HS, lmax_topo, HC_topo, HS_topo)
#    val = Get_Geoid_Height  (R, phi, theta, lmax, HC, HS)

    return val


def TEST_ellipsoid_corr():
    """
    replace range(0, l+1) by (0, 1) in m loop of Get_Geoid_Height
    replace if(mins<=0): statement in of init_grid
    """
    HC, HS = imp.Fetch_Coef()
    lmax = 29
    c_cos = np.zeros((lmax+1, lmax+1))
    c_sin = np.zeros((lmax+1, lmax+1))
    for i in range(2, 11):
#        c_cos[i, 0] = Cosine_Correction2(i)
#        if (i%2==1): c_cos[i, 0] = HC[i, 0]
        if (i%2==0): c_cos[i, 0] = Cosine_Correction(i)


    mins = 0
    limits = np.array([-180, 180, -90, 90])

    G_Grid, G_Long, G_Lat =  Gen_Grid (mins, Get_Geoid_Height,
#                                        [lmax, c_cos, c_sin],
                                        [lmax, HC, HS],
                                         limits)
#    FIG, AX = emap.Make_Map(limits=limits)#proj = ccrs.Mollweide)
#    CBAR = emap.Plot_contourf(G_Grid, G_Long, G_Lat, AX)
    FIG, AX = emap.Make_Map_3D()
    _ = emap.Plot_surface(G_Grid, G_Long, G_Lat, AX)
    AX.set_zlabel("Geoid Height (m)",rotation=90)
    plt.title("Ellipsoid removal residual")
    return c_cos, G_Grid




# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

#    TEST_plotGeoPot_radius()
    g = TEST_lmax_loop_lat_line()
#    g = TEST_lmax_loop_long_line()
#    TEST_Get_isopot ()
#    TEST_Cosine_corr()

#    val = TEST_high_lmax()

#    th, ph, gr = TEST_gen_grid()
#    cc, gg = TEST_ellipsoid_corr()


    print("\nGH_harmonics done")