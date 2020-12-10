#!/usr/bin/env python
# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from math import factorial


def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] / (rate**deriv) * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


# In[43]:


t = np.linspace(-4, 4, 500)
y = np.sin(t) + np.random.normal(0, 0.025, t.shape) 
d2ysg = savitzky_golay(y, window_size=50, order=4, deriv = 2, rate = t[0] - t[1])

import matplotlib.pyplot as plt
plt.plot(t, y, label='Noisy signal')
#plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
plt.plot(t, d2ysg, 'r', label='Filtered signal')
plt.legend()
plt.show()

