import numpy as np
import scipy.signal as sp
from scipy.special import lpmn
import PG_globalVars as gv
import statsmodels.api as sm

G = 1
M = 1
R = 1

SG_window_size = 5  #window_size : the length of the filter window (i.e. the number of coefficients) for SG filter
SG_order = 2        #order : the order of the polynomial used to fit the samples for SG filter
ar_order = 4


# COMMENT : If we wanted to do this properly, we would do an analysis of the error introduced by the Savitzky-Golay filter,
#           and we would compare it to the error introduced by the finite difference method
# https://personal.math.ubc.ca/~cbm/aands/abramowitz_and_stegun.pdf
# https://www.mat.univie.ac.at/~westra/associatedlegendrefunctions.pdf


def SG_filt_cartesian(data, d_order):
    """ parameters :
            - data: positions over time in cartesians coordinates

        output :
            - the approximated accelerations in cartesians coordinates // Bezdek 2014"""

    [x, y, z, time] = data

    acc_x = sp.savgol_filter(x, SG_window_size, SG_order, d_order)
    acc_y = sp.savgol_filter(y, SG_window_size, SG_order, d_order)
    acc_z = sp.savgol_filter(z, SG_window_size, SG_order, d_order)

    return acc_x, acc_y, acc_z


def SG_d2_matrix_cartesian(data):
    [x, y, z, time] = data
    time_step = time[1] - time[0]

    SGM = np.zeros((3*len(x), 3*len(x)))
    coeffs = sp.savgol_coeffs(SG_window_size, SG_order, deriv=2, delta=time_step)

    for i in range(len(x)):
        for k in range(SG_window_size):

            if len(x) > i + k - (SG_window_size - 1) // 2 >= 0:
                SGM[i][(i + k - (SG_window_size - 1) // 2)] = coeffs[k]                         # construction of the X lines
                SGM[i + len(x)][(i + k - (SG_window_size - 1) // 2)] = coeffs[k]        # construction of the Y lines
                SGM[i + 2 * len(x)][(i + k - (SG_window_size - 1) // 2)] = coeffs[k]    # construction of the Z lines

    return SGM


def linear_transformation_matrix(data):

    """ parameters :
            - data : GPS data

        output :
            - linear transformation matrix"""

    SGM = SG_d2_matrix_cartesian(data)
    cov_matrix = np.matmul(SGM, np.transpose(SGM))
    T = np.linalg.cholesky(cov_matrix)

    return np.linalg.inv(T)


def vector_from_spherical_to_cartesian(vector, theta, phi):

    """ parameters :
            - vector: R3 vector in spherical frame
            - theta, phi : sphericals coordinates

        output :
            - cartesian version of the R3 vector"""

    assert (len(vector) == 3)

    cartesian_vector = [np.sin(theta) * np.cos(phi) * vector[0]
                        + np.cos(theta) * np.cos(phi) * vector[1]
                        - np.sin(phi) * vector[2],
                        np.sin(theta) * np.sin(phi) * vector[0]
                        + np.cos(theta) * np.sin(phi) * vector[1]
                        + np.cos(phi) * vector[2],
                        np.cos(theta) * vector[0]
                        - np.sin(theta) * vector[1]
                        ]

    return cartesian_vector


def gradient_Vnm_c_spherical(n, m, r, theta, phi):

    """ parameters :
            - n, m : Legendre Polynome parameters
            - r, theta, phi : spherical coordinates

        output :
            - sphrerical gradient of potential Vnm(c) // Bezdek 2014
    """

    #full and semi normalisation
    normalization_factor = np.sqrt((2 * n + 1) * np.math.factorial(n - m) / (2 * np.math.factorial(n + m)))
    normalization_factor_schmidt = np.sqrt(2 * np.math.factorial(n - m) / np.math.factorial(n + m))

    #calculation of Legendre associated functions with parameters in [0,n] X [0,m], dpmn is the array of the derivatives values
    pnm, dpnm = lpmn(m, n, np.cos(theta))

    #gradient computation
    Gradient = [
        - pnm[m, n] * np.cos(m * phi) * G * M * (n + 1) * R ** n / (r ** (n + 2)) * normalization_factor * (-1)**m,
        - 1 / r * np.cos(m * phi) * np.sin(theta) * dpnm[m, n] * G * M * R ** n / (r ** (n + 1)) * normalization_factor * (-1)**m,
        - 1 / (r * np.sin(theta)) * np.sin(m * phi) * pnm[m, n] * np.cos(m * phi) * G * M * R ** n / (r ** (n + 1)) * normalization_factor * (-1)**m
    ]

    return Gradient


def gradient_Vnm_s_spherical(n, m, r, theta, phi):
    """ parameters :
            - n, m : Legendre Polynome parameters
            - r, theta, phi : spherical coordinates

        output :
            - sphrerical gradient of potential Vnm(s) // Bezdek 2014"""

    # full and semi normalisation
    normalization_factor = np.sqrt((2 * n + 1) * np.math.factorial(n - m) / (2 * np.math.factorial(n + m)))
    normalization_factor_schmidt = np.sqrt(2 * np.math.factorial(n - m) / np.math.factorial(n + m))

    # calculation of Legendre associated functions with parameters in [0,n] X [0,m], dpmn is the array of the derivatives values
    pnm, dpnm = lpmn(m, n, np.cos(theta))

    #Gradient computation
    Gradient = [
        - pnm[m, n] * np.sin(m * phi) * G * M * (n + 1) * R ** n / (r ** (n + 2)) * normalization_factor * (-1)**m,
        - 1 / r * np.sin(m * phi) * np.sin(theta) * dpnm[m, n] * G * M * R ** n / (r ** (n + 1)) * normalization_factor * (-1)**m,
        + 1 / (r * np.sin(theta)) * np.cos(m * phi) * pnm[m, n] * np.cos(m * phi) * G * M * R ** n / (r ** (n + 1)) * normalization_factor * (-1)**m,
    ]

    return Gradient


def linear_regression_solution(accelerations, data, order):

    """ parameters :
            -  accelerations : the approximated accelerations in cartesians coordinates [[acc_x],[acc_y], [acc_z]
            -  data : GPS data in spherical coordinates
            _  order of the spherical approximation

        output :
            - array containing the coefficient of the spherical harmonics representation // Bezdek 2014 + wikipedia Geopot. model
            ex: [C_00, S_00, C_10, S_10, C_11, S_11]"""

    [r, theta, phi, time] = data

    # test to assure we have coherents number of R3 vectors
    assert (len(accelerations[0]) == len(r))

    n_coeff = (order + 1) * (order + 2)  # number of harmonics coefficients
    n_mesured_values = len(accelerations) // 3  # number of acc mesured

    # matrix of the linear problem :  Acc = design_matrix * Coeff
    design_matrix = np.zeros((len(accelerations), n_coeff))

    c = 0
    for n in range(order+1):
        for m in range(n+1):
            vnm_c = [[], [], []]
            vnm_s = [[], [], []]

            #SG filtering of gradients data
            for i in range(n_mesured_values):
                v_c = vector_from_spherical_to_cartesian(gradient_Vnm_c_spherical(n, m, r[i], theta[i], phi[i]), theta[i], phi[i])
                v_s = vector_from_spherical_to_cartesian(gradient_Vnm_s_spherical(n, m, r[i], theta[i], phi[i]), theta[i], phi[i])
                for k in range(3):
                    vnm_c[k].append(v_c[k])
                    vnm_s[k].append(v_s[k])

            vnm_c_f = SG_filt_cartesian(vnm_c, 0)
            vnm_s_f = SG_filt_cartesian(vnm_s, 0)

            # construction of the matrix
            for j in range(n_mesured_values):
                for l in range(3):
                    design_matrix[j + l * n_mesured_values / 3][c] = vnm_c_f[l][j]
                    design_matrix[j + l * n_mesured_values / 3][c+1] = vnm_s_f[l][j]

    # construction of the accelerations vector and linear transformation matrix
    acc_vector = accelerations[0] + accelerations[1] + accelerations[2]
    lt_matrix = linear_transformation_matrix(data)

    #linear transformation of the problem
    transformed_design_matrix = np.matmul(lt_matrix, design_matrix)
    transformed_acc_vector = np.matmul(lt_matrix, acc_vector)

    return np.linalg.lstsq(transformed_design_matrix, transformed_acc_vector)


def GeoPot(coefficients, theta, phi, order):

    """ parameters :
            - coefficients : harmonics coefficients
            - theta, phi : spherical coordinates
            - order : order of the harmonic approx

        output :
            - geopotential """

    latitude = np.pi/2 - theta
    r_e = np.sqrt(((gv.a_e**2 * np.cos(latitude))**2 + (gv.b_e**2 * np.sin(latitude))**2)/((gv.a_e * np.cos(latitude))**2 + (gv.b_e * np.sin(latitude))**2))
    gpot = 0
    i = 0

    for n in range(order+1):
        for m in range(n+1):

            normalization_factor = np.sqrt((2 * n + 1) * np.math.factorial(n - m) / (2 * np.math.factorial(n + m)))

            pnm, dpnm = lpmn(m, n, np.cos(theta))
            gpot += gv.GM_e * gv.a_e**n / r_e**(n+1) * pnm[m, n] * (coefficients[i] * np.cos(m * phi) + coefficients[i+1] * np.sin(m * phi)) * (-1)**m

            i += 2

    return r_e, gpot


def GeoPotDiff(coefficients, theta, phi, order):
    """ parameters :
             - coefficients : harmonics coefficients
             - theta, phi : spherical coordinates
             - order : order of the harmonic approx

         output :
             - geopotential anomalies """

    latitude = np.pi/2 - theta
    r_e = np.sqrt(((gv.a_e**2 * np.cos(latitude))**2 + (gv.b_e**2 * np.sin(latitude))**2)/((gv.a_e * np.cos(latitude))**2 + (gv.b_e * np.sin(latitude))**2))
    gpot = 0
    i = 8

    for n in range(2, order+1):
        for m in range(n+1):

            if n != 2 or m != 0:

                normalization_factor = np.sqrt((2 * n + 1) * np.math.factorial(n - m) / (2 * np.math.factorial(n + m)))

                pnm, dpnm = lpmn(m, n, np.cos(theta))
                gpot += gv.GM_e * gv.a_e**n / r_e**(n+1) * pnm[m, n] * (coefficients[i] * np.cos(m * phi) + coefficients[i+1] * np.sin(m * phi)) * (-1)**m

                i += 2

    return r_e, gpot


def linear_transformation_matrix_2(accelerations):
    x_cov = sm.tsa.acovf(accelerations[0],nlag=len(accelerations[0]))
    y_cov = sm.tsa.acovf(accelerations[1], nlag=len(accelerations[1]))
    z_cov = sm.tsa.acovf(accelerations[2], nlag=len(accelerations[2]))

    M = np.zeros((3*len(accelerations[0]),3*len(accelerations[0])))
    pass



#########################################################################################
#                                  NOT USED FUNCTIONS                                   #
#########################################################################################


def sav_filt_spherical(data, window_size, order):

    [r, theta, phi, time] = data
    data_len = len(r)

    r_dot = sp.savgol_filter(r, window_size, order, 1)
    theta_dot = sp.savgol_filter(theta, window_size, order, 1)
    phi_dot = sp.savgol_filter(phi, window_size, order, 1)

    r_dot_dot = sp.savgol_filter(r, window_size, order, 2)
    theta_dot_dot = sp.savgol_filter(theta, window_size, order, 2)
    phi_dot_dot = sp.savgol_filter(phi, window_size, order, 2)

    r_sq_times_theta_dot = [r[i] ** 2 * theta_dot[i] ** 2 for i in range(data_len)]
    time_deriv_r_sq_times_theta_dot = sp.savgol_filter(r_sq_times_theta_dot, window_size, order, 1)
    r_sq_times_sin_sq_theta_times_phi_dot = [r[i] ** 2 * (np.sin(theta[i]) ** 2) * phi_dot[i] for i in range(data_len)]
    time_deriv_r_sq_times_sin_sq_theta_times_phi_dot = sp.savgol_filter(r_sq_times_sin_sq_theta_times_phi_dot, window_size, order, 1)

    accel_r = [r_dot_dot[i] - r[i] * (theta_dot[i] ** 2 + (phi_dot[i] ** 2) * (np.sin(theta[i]) ** 2)) for i in range(data_len)]
    accel_theta = [time_deriv_r_sq_times_theta_dot[i] / r[i] - r[i] * np.sin(theta[i]) * np.cos(theta[i]) * (phi_dot[i] ** 2) for i in range(data_len)]
    accel_phi = [time_deriv_r_sq_times_sin_sq_theta_times_phi_dot[i] / (r[i] * np.sin(theta[i])) - 2 * r[i] * np.cos(theta[i]) * phi_dot[i] * theta_dot[i] for i
                 in range(data_len)]

    return accel_r, accel_theta, accel_phi
