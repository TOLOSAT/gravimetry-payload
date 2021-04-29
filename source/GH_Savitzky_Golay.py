#!/usr/bin/env python
# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from math import factorial
from scipy.linalg import toeplitz


import GH_import       as imp
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


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    ###-----------------###
    #Use the Savitzky-Golay method to determine the accelerations from a cloud of positions (y)
    ###-----------------###
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


def savitzky_golay_mat(Nmeasures, order, window_size, deriv = 0, rate = 1):

    order_range = range(order+1)
    half_window = (window_size - 1)//2

    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] / (rate**deriv) * factorial(deriv)
    print(m)

    index = (half_window  + min(Nmeasures - 1,  half_window))

    print(m[0:2])

    M = toeplitz(np.concatenate((m[half_window : index + 1 ],np.zeros(max((0,Nmeasures - half_window - 1 ) )) )))

    return M







def correcRef(Pos,Vit, Acc, omegaTerre = 7292115E-11):
    '''Corrects acceleration for coriolis and centrifugal forces in terrestrial
    referential'''
    AccCorrec = np.copy(Acc)
    AccCorrec[:,0] += - 2*omegaTerre*Vit[:,1] - omegaTerre**2*Pos[:,0]
    AccCorrec[:,1] +=   2*omegaTerre*Vit[:,0] - omegaTerre**2*Pos[:,1]

    return AccCorrec


def testCorrecRef(file_name, days=0.7, data_path="../data"):
    '''temporary code that tests correction of coriolis etc and  savitsky Golay'''
    Eph = np.loadtxt(f"{data_path}/{file_name}")
    t = np.array(Eph[::20,0]) #time in seconds
    x = np.array(Eph[::20,1]) #  \
    y = np.array(Eph[::20,2]) #  | cordinates, in km
    z = np.array(Eph[::20,3]) # /
    dt = np.int(t[1]*100)/100
    L = np.int(days*(86400/dt))
    # convert coord system and shorten array if needed
    pts = np.transpose(np.array([x,y,z]))
    if L >= len(pts):
        L = len(pts) # this is not necessary in python
    Time = t[:L]


    Vx = savitzky_golay(x, 50, 2 , 1, dt)
    Vy = savitzky_golay(y, 50, 2 , 1, dt)
    Vz = savitzky_golay(z, 50, 2 , 1, dt)
    Vit = np.transpose(np.array([Vx,Vy,Vz]))

    Ax = savitzky_golay(x, 50, 2 , 2, dt)
    Ay = savitzky_golay(y, 50, 2 , 2, dt)
    Az = savitzky_golay(z, 50, 2 , 2, dt)


    Acc = np.transpose(np.array([Ax,Ay,Az]))

    AccCor = correcRef(pts,Vit,Acc)



    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x = x[::10]
    y = y[::10]
    z = z[::10]

    Ax = Ax[::10]
    Ay = Ay[::10]
    Az = Az[::10]

    AccCor = AccCor[::10,:]

    ax.quiver(x,y,z,Ax,Ay,Az,   length=2500, normalize = True)
    ax.quiver(x,y,z,AccCor[:,0],AccCor[:,1],AccCor[:,2],colors = 'purple',   length=2500, normalize = True)

    plt.show()




# In[43]:


# t = np.linspace(-4, 4, 500)
# y = np.sin(t) + np.random.normal(0, 0.025, t.shape)
# d2ysg = savitzky_golay(y, window_size=250, order=5, deriv = 2, rate = t[0] - t[1])
#
# import matplotlib.pyplot as plt
# plt.plot(t, y, label='Noisy signal')
# #plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
# plt.plot(t, d2ysg, 'r', label='Filtered signal')
# plt.legend()
# plt.show()