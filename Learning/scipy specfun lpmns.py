# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:07:08 2020

@author: Xavier
"""
import numpy as np
from scipy.special import lpmn

def LPMNS(M, N, X):
    """ 
    trying to understand this fortran script:
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
    return PM.T, PD.T


def Pol_Legendre (l, m, x):
    """
    function from scipy directly
    """
    Plm_z, Plm_dz = lpmn(m, l, x)
    return Plm_z, Plm_dz # use Plm_z[m, l]
 
    
# =============================================================================

m = 5
l = 7
x = 1.1

PM, PD     = LPMNS(m, l, x)
PM_for, PD_for = LPMNA(m, l, x)
PM_sp, PD_sp = lpmn (m, l, x)
Err = PM_sp - PM_for









