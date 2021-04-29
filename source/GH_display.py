# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 01:27:06 2020

@author: Xavier de Labriolle, Antoine Bendimerad

Last edit : 29/02/2020

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
    
    The purpose of this script is to Display various graphs, plots to show:
        "True" and "solved" Geoids
        Topology ? - maybe make a separate file
        Differences in coefficients
        
        
==============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import matplotlib.pyplot as plt

import numpy as np
from numpy import pi, sin, cos

from GH_solve2 import *
from GH_generate import *


# =============================================================================
# DISPLAY FUNCTIONS  
# =============================================================================
def Plot_Array_Diff(HS_nm_slv, HC_nm_slv, fig_num = 6):
    print("plotting coeff difference")
    
    
    #resize the official coef
    HC, HS = Fetch_Coef()
    HS_nm_sz = HS[:len(HS_nm_slv), :len(HS_nm_slv)]
    HC_nm_sz = HC[:len(HC_nm_slv), :len(HC_nm_slv)]

    #subtract calculated coeffs
    HS_nm_sz -= HS_nm_slv
    HC_nm_sz -= HC_nm_slv
    
    fig_HC = plt.figure(fig_num)
    plt.clf()
    plt.suptitle("Harmonic coeff diference between official and solved; degree: "+str(len(HS_nm_sz)-1))
    
    for n in range (0, len(HC_nm_sz)):
        Ms_n = np.arange(0, n+1)
        
        HC_ni = HC_nm_sz[n, :n+1]
        HS_ni = HS_nm_sz[n, :n+1]
        
        plt.subplot(211)
        plt.plot(Ms_n, HC_ni,'-*', label='n='+str(n))
        
        plt.subplot(212)
        plt.plot(Ms_n, HS_ni,'-*', label='n='+str(n))
    
    plt.subplot(211)
    plt.ylabel("COSINE coeff diff")
    plt.grid(True)
#    plt.xlabel("order m of derivation (log)")
#    plt.ylabel("value of HC_nm")
    plt.legend(loc = 'upper right', title = 'Degree n', fontsize = 5)
    
    plt.subplot(212)
    plt.ylabel("SINE coeff diff")
    plt.grid(True)
    plt.xlabel("order m of derivation (log)")
#    plt.ylabel("value of HS_nm")
#    plt.legend(loc = 'lower right', title = 'Degree n', fontsize = 5)
    
    plt.show()
    

# =============================================================================
# MAIN 
# =============================================================================

"""
Process : 
    1.  Generate the acceleration values (from orbit data or raw simulation)
    2.  Solve for the coefficients using the functions in this script
    3.  Generate a Geoid map from solved coefficients
    4.  Compare Geoids or coefficients
"""

print("\nGH_display done")
    

