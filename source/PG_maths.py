import numpy as np
import scipy.signal as sp

# Approx accel : usage of a Savitzky-Golay filter to smooth the data
# data : Positions over time, r, theta, phi according to WGS84
# window_size : the length of the filter window (i.e. the number of coefficients)
# order : the order of the polynomial used to fit the samples
# returns : the approximated accelerations in the r, theta, phi coordinates
# COMMENT : If we wanted to do this properly, we would do an analysis of the error introduced by the Savitzky-Golay filter
#           and we would compare it to the error introduced by the finite difference method
#           
def sav_filt(data, window_size, order):
    [r, theta, phi] = data
    data_len = len(r)

    r_dot = sp.savgol_filter(r, window_size, order, 1)
    theta_dot = sp.savgol_filter(theta, window_size, order, 1)
    phi_dot = sp.savgol_filter(phi, window_size, order, 1)

    r_dot_dot = sp.savgol_filter(r, window_size, order, 2)
    theta_dot_dot = sp.savgol_filter(theta, window_size, order, 2)
    phi_dot_dot = sp.savgol_filter(phi, window_size, order, 2)

    r_sq_times_theta_dot = [r[i]**2 * theta_dot[i]**2 for i in range(data_len)]
    time_deriv_r_sq_times_theta_dot = sp.savgol_filter(r_sq_times_theta_dot, window_size, order, 1)
    r_sq_times_sin_sq_theta_times_phi_dot = [r[i]**2 * (np.sin(theta[i])**2) * phi_dot[i] for i in range(data_len)]
    time_deriv_r_sq_times_sin_sq_theta_times_phi_dot = sp.savgol_filter(r_sq_times_sin_sq_theta_times_phi_dot, window_size, order, 1)

    
    accel_r = [r_dot_dot[i] - r[i]*(theta_dot[i]**2 + (phi_dot[i]**2)*(np.sin(theta[i])**2)) for i in range(data_len)]
    accel_theta = [time_deriv_r_sq_times_theta_dot[i]/r[i] - r[i]*np.sin(theta[i])*np.cos(theta[i])*(phi_dot[i]**2) for i in range(data_len)]
    accel_phi = [time_deriv_r_sq_times_sin_sq_theta_times_phi_dot[i]/(r[i]*np.sin(theta[i])) - 2*r[i]*np.cos(theta[i])*phi_dot[i]*theta_dot[i] for i in range(data_len)]

    return accel_r, accel_theta, accel_phi


