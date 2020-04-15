# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 01:02:16 2020

@author: Xavier

This script is an attempt to make a dichotomy approach to find the "z_e" input 
to a "f" function xith " *[x,y], z" inputs
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin, cos
from time import sleep

# =============================================================================
def f (x, y):
    Z = x*cos(pi*x) + y*sin(pi*y)
    return Z

def G (x, y, z):
    w = 0.7*x + 1.2*y + 5.6*z
    return w

def dichotomy_grad (f, in_first, z_e, in_after, w_0, de, g):
    """
    f :    function for f(*in_first,z,*in_after) = w
    in_first, in_after : f function variables that come before and after z_e
    z_e :     value to be tested
    w_0:    target value for w
    de :    delta error  
    g :     gradient
    """
    w_i = f(*in_first, z_e, *in_after)
    di = w_0 - w_i
    z_i = z_e
    
    c = 0
    while (abs(di) >= de):
        c+=1
        
        z_i += di/g
        w_i = f(*in_first, z_i, *in_after)
        di = w_0 - w_i
        
#        sleep(1); print(f"w_0={w_0}; z_i={z_i}; w_i={w_i}; di={di}; add={di/g}; "); 
#    print(f"dichotomy_grad: {c} steps")
    return z_i

    
def dichotomy_step (f, in_first, z_e, in_after, w_0, de, st):
    """ see dichotomy_grad """
    w_i = f(*in_first, z_e, *in_after)
    di_1 = w_0 - w_i
    di_2 = w_0 - w_i
    z_i = z_e
    
    c=0
    while (abs(di_1) >= de):  
        c+=1
        
        di_1 = di_2
        w_i = f(*in_first, z_i, *in_after)
        di_2 = w_0 - w_i    
        
        if (di_1*di_2 < 0): 
            st = st/2
        
        z_i += st * di_1/abs(di_1)   
        
#        sleep(1); print(f"w_0={w_0}; z_i={z_i}; w_i={w_i}; di_1={di_1}; di_2={di_2}; st={st}"); 
#    print(f"dichotomy_step: {c} steps")
    return z_i

   
def TEST_dichotomy():
    x=1; y=2; z=2; w=45; de=0.001; gr=10; st=3    
    z_i = dichotomy_grad(G, [x,y], z, [], w, de, gr)
    print (f"G(x,y,z_i)={w}: z_i = {z_i} ; G(x,y,z_i)={G(x,y,z_i)}")    
    z_i = dichotomy_step(G, [x,y], z, [], w, de, st)
    print (f"G(x,y,z_i)={w}: z_i = {z_i} ; G(x,y,z_i)={G(x,y,z_i)}")




# =============================================================================
# 
# =============================================================================
TEST_dichotomy()



'''
w=0; de=0.001; gr=6 

sz = 20
Xs = np.linspace(-5, 5, sz)
Ys = np.linspace(-5, 5, sz)
GXs, GYs = np.meshgrid(Xs, Ys)

GZs = np.zeros((sz, sz))

for i in range(0, len(Xs)):
    x = Xs[i]

    for j in range(0, len(Ys)):
        y = Ys[j]
        
#        GZs[j,i] = f(x,y)
#        GZs[j,i] = g(x,y,1+x)
        GZs[j,i] = dichotomy_grad(G, x, y, 1, w, de, gr)
    
    print(f"x={x} done")

plt.figure()
plt.clf()
plt.contourf(GXs, GYs, GZs, 50, cmap=plt.get_cmap("jet"))
plt.colorbar()
'''







