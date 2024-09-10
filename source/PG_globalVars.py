import numpy as np

g = 9.80665                         # m/s^2 reference acceleration

# WGS84 reference ellipsoid model
ref_e = "WGS 84"                    # https://en.wikipedia.org/wiki/World_Geodetic_System
GM_e = 3986004.418E8                # m^3/s^2 : standard gravitational parameter
wo = 7292115E-11                    # rad/s : angular velocity of Earth
a_e = 6378137.00 # m : equatorial radius or ellipsoid model
f = 1/298.257223563 # flat parameter

b_e = a_e * (1-f)                                   # m : polar radius
E = np.sqrt(a_e**2 - b_e**2)                        # linear eccentricity
e_1 = E/a_e
e_2 = E/b_e
m = wo**2 * a_e*2 * b_e / GM_e                      # just to simplify the code
g_a = GM_e/(a_e*b_e) * (1 - 3/2*m - 3/14*e_2*m)     # m/s^2 : gravity acc. at equator
g_b = GM_e/(a_e**2) * (1 - m - 3/7*e_2*m)           # m/s^2 : gravity acc. at poles

# EGS2008 potential model
ref_g = "EGS2008"
a_g = 6378136.3                     # m : Reference radius for the potential model
GM_g = 3986004.415E8                # m^3/s^2 : standard gravitational parameter in the potential model


coeffs = [1, 0,
          0, 0, 0, 0,
          -0.10826360229840E-02, 0, -0.24140000522221E-09, 0.15430999737844E-08,  0.15745360427672E-05, -0.90386807301869E-06,
          0.25324353457544E-05, 0, 0.21927988018965E-05, 0.26801189379726E-06, 0.30901604455583E-06, -0.21140239785975E-06, 0.10055885741455E-06, 0.19720132389889E-06,
          0.16193312050719E-05, 0, -0.50872530365024E-06, -0.44945993508117E-06, 0.78412230752366E-07, 0.14815545694714E-06, 0.59215743214072E-07, -0.12011291831397E-07, -0.39823957404129E-08, 0.65256058113396E-08,
          0.22771610163688E-06, 0, -0.53716510187662E-07, -0.80663463828530E-07, 0.10559053538674E-06, -0.52326723987632E-07, -0.14926153867389E-07, -0.71008771406986E-08, -0.22979123502681E-08, 0.38730050770804E-09, 0.43047675045029E-09, -0.16482039468636E-08]

