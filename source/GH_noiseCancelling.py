"""
@authors:
# =============================================================================
 Information:
    This program filters the signal to erase the noise
todo:
    Code the noise cancelling
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
import scipy as sc
import scipy.signal as sig
import math as m
import cmath as cm
from time import gmtime, strftime
from pykalman import KalmanFilter

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
import GH_Savitzky_Golay as sg

# =============================================================================
# GLOBAL VARIABLES
# =============================================================================

# =============================================================================
# FUNCTIONS - NOISE FILTERING
# =============================================================================

def polynomialFilter(Nmeasures, t):
    """Computes the square matrix from the coefficients of the polynomial filter F.
    It comes from the Savitzky-Golay filter.
    Nmeasures : number of measures
    t : time data"""

    F = sg.savitzky_golay_mat(Nmeasures, 8, 9, 1, t[1] - t[0])
    return F




def cholesky(F):
    """Performs a Cholesky decomposition of the covariance matrix FFt, which must
    be a symmetric and positive definite matrix. The function
    returns the lower variant triangular matrix inverse, W."""

    FFt = np.dot(F, np.transpose(F))

    n = len(FFt)

    # Create zero matrix for L
    T = [[0.0] * n for i in xrange(n)]

    # Perform the Cholesky decomposition
    for i in xrange(n):
        for k in xrange(i+1):
            tmp_sum = sum(T[i][j] * T[k][j] for j in xrange(k))

            if (i == k): # Diagonal elements
                T[i][k] = np.sqrt(FFt[i][i] - tmp_sum)
            else:
                T[i][k] = (1.0 / T[k][k] * (FFt[i][k] - tmp_sum))
    W = sc.linalg.inv(T)

    return W


def karmanFilter(acc, M):
    """Generates the Kalman Filter linear transformation matrix
    acc : acceleration data
    M : gradient matrix"""



def antiNoise(W, acc, M):
    """W : the linear transformation matrix
    acc : acceleration data
    M : gradient matrix
    We deriive the filter and apply the least squares with the equation acc = M.x + epsilon
    with epsilon the error and x the Snm/Cnm vector
    Applying the linear transformation, we got acc* = M*.x + epsilon*.
    Applying the filter, we got acc** = M**.x + epsilon**"""








"""Next : implement an alogrithm to delete the noise : derive the noise from the noise variance matrix with the help of W"""









