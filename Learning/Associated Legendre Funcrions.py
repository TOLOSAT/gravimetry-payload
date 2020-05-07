# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:07:08 2020

@author: Xavier
"""
import numpy as np
from numpy import cos, sin
from scipy.special import lpmn
import math



def LPMNS(M, N, X):
    """ 
    trying to understand this Fortran script:
    https://github.com/scipy/scipy/blob/master/scipy/special/specfun/specfun.f
    line 1491
C       ========================================================
C       Purpose: Compute associated Legendre functions Pmn(x)
C                and Pmn'(x) for a given order
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(n) --- Pmn(x)
C                PD(n) --- Pmn'(x)
C       ========================================================
    """
#    lines 1503 - 1507
    PM = np.zeros(N+1)
    PD = np.zeros(N+1)
#    lines 1508 - 1525
    if (abs(X) == 1.0):
        for k in range (0, N+1):
            if (M == 0):
                PM[k] = 1
                PD[k] = 0.5 * k * (k+1)
                if (X < 0.0):
                    PM[k] = (-1)**k * PM[k]
                    PD[k] = (-1)**(k+1) * PD[k]
            elif (M == 1.0): # line 1517
                PD[k] = 1E300
            elif (M == 2):
                PD[k] = 0.25 * (k+2) * (k+1) * k * (k-1)
                if (X < 0):
                    PD[k] = (-1)**(k+1) * PD[k]
#    lines 1526 - 1534
    X0 = abs(1 - X*X)
    PM0 = 1
    PMK = PM0
    for k in range (1, M+1):
        PMK = (2*k+1) * np.sqrt(X0) * PM0
        PM0 = PMK
    PM1 = (2*M + 1) * X * PM0
    PM[M] = PMK
    PM[M+1] = PM1
#    lines 1535 - 1540
    for k in range (M+2, N+1):
        PM2 = ( (2*k-1)*X*PM1 - (k+M-1)*PMK ) / (k - M)
        PM[k]= PM2
        PMK = PM1
        PM1 = PM2
    PD[0]=( (1-M)*PM[1] - X*PM[0] ) / (X*X - 1)    
#    lines 1541 - 1547
    for k in range (1, N+1):
        PD[k] = ( k*X*PM[k] - (k+M)*PM[k-1] ) / (X*X - 1)
    for k in range (1, N+1):
        PM[k] = (-1)**M * PM[k]
        PD[k] = (-1)**M * PD[k]
    return PM, PD

def LPMNA(M, N, X):
    """ just making useful arrays as I like them """
    PM = np.zeros((N+1, M+1))
    PD = np.zeros((N+1, M+1))
    for m in range (0, M+1):
        PM[:,m], PD[:,m] = LPMNS(m, N, X)
    return PM, PD



def Pol_Legendre (N, M, X):
    """
    function from scipy directly
    """
    Plm_z, Plm_dz = lpmn(M, N, X)
    
    for n in range(0, Plm_z.shape[0]):
        for m in range (0, n+1):
            Plm_z[m, n] = Plm_z[m, n] * Normalize(n, m)
    
    return Plm_z.T, Plm_dz.T # use Plm_z[m, l]
 

def Normalize (l, m):
    """
    Returns the normalization coefficient of degree l and order m
    Equation obtained from
    """
    k = 2
    if (m == 0) : k = 1
    P1 = math.factorial(l - m)
    P2 = k*(2*l + 1)
    P3 = math.factorial(l + m)
    N = np.sqrt(P1*P2/P3)
    return N
def ALF_norm_gcb (N, M, phi):
    """
    returns an array[m+1,n+1] of the values of the Associated Legendre Function
    of all integer degrees l and order m, at point x
    Array is normalized, equations from the geoid cook book
    This method is called the "standard forward colums method" explained in: 
    https://link.springer.com/article/10.1007/s00190-002-0216-2
    """
#    phi = conv.geodes2geocen(phi_gd)
    t = sin(phi)
    u = cos(phi)
    
    POL = np.zeros((N+1, M+1))
    POL[0,0] = 1
    POL[1,0] = t 
    POL[1,1] = u * np.sqrt(3)
#    POL[2,0] = 3/2*t**2-1/2 * Normalize(2,0)
#    POL[2,1] = 3*u*t * Normalize(2,1)
#    POL[2,2] = 3*u**2 * Normalize(2,2)
    
    for n in range(2, M+1):
        POL[n,n] = u*np.sqrt((2*n+1)/(2*n))*POL[n-1,n-1]
        
    a_nm = lambda n, m : np.sqrt( (2*n+1)*(2*n-1) / ((n-m)*(n+m)) )
    b_nm = lambda n, m : np.sqrt( (2*n+1)*(n+m-1)*(n-m-1) / ((n-m)*(n+m)*(2*n-3)) )
    
    for n in range (1, N+1): # m
        POL[n, 0] = a_nm(n,0) * t * POL[n-1,0]

    for m in range (1, M+1): # m
        for n in range (m+1, N+1): # n
            POL[n, m] = a_nm(n,m)*t*POL[n-1,m] - b_nm(n,m)*POL[n-2,m]
    
    return POL



    
# =============================================================================

N = 70
M = N - 1
X = 1

#PM, PD         = LPMNS(M, N, cos(X))
PM_for, PD_for = LPMNA(M, N, cos(X))
PM_spy, PD_spy = Pol_Legendre (N, M, cos(X))
#Err = PM_spy - PM_for

PM_gcb = ALF_norm_gcb(N, M, X)

t = sin(X)
u = cos(X)

POL = np.zeros((N+1, M+1))
POL[0,0] = 1
POL[1,0] = t 
POL[1,1] = u * np.sqrt(3)
POL[2,0] = (3/2*t**2 - 1/2) * Normalize(2,0)
POL[2,1] = 3*u*t * Normalize(2,1)
POL[2,2] = 3*u**2 * Normalize(2,2)

aa = PM_gcb - POL




