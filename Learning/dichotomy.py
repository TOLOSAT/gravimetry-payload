# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 01:02:16 2020

@author: Xavier
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin, cos
from time import sleep

# =============================================================================
def f (x, y):
    Z = x*cos(pi*x) + y*sin(pi*y)
    return Z

def g (x, y, z):
    w = 0.7*x + 1.2*y + 0.5*z*x
    return w

def dichotomy_grad (f, x, y, z, w_0, de, g):
    """
    f :    function for f(x,y,z) = w
    x, y :  position
    z :     value to be tested
    w_0:    target value for w
    de :    delta error  
    g :     gradient
    """
    w_e = f(x, y, z)
    di = w_e - w_0
    z_i = z
    
    while ( abs(di) >= de):
        print(f"di={di}; ")
        z_i += di/g
        w_i = f(x, y, z_i)
        di = w_i - w_0

    return z_i

    
def dichotomy_step (f, x, y, z, w_0, de, st):
    """
    f :    function for f(x,y,z) = w
    x, y :  position
    z :     initial value to be tested
    w_0:    target value for w
    de :    delta error  
    st :    step size
    """
    w_e = f(x, y, z)
    g = 10
    di_1 = w_e - w_0
    di_2 = w_e - w_0
    z_i = z
    
    while ( abs(di_1) >= de):        
        di_1 = di_2
        w_e = f(x, y, z_i)
        di_2 = w_e - w_0    
        
        if (di_1*di_2 < 0): 
            st = st/2
        
        z_i -= st * di_1/abs(di_1)    
#        print(f"z_i={z_i}; w_0={w_0}; w_e={w_e}; di_1={di_1}; di_2={di_2}; st={st}"); 

    return z_i


   
def TEST_dichotomy():
    x=1; y=2; z=2; w=0; de=0.01; gr=10; st=1
#    z_i = dichotomy_grad(g, 1, 2, 3, -4, 0.1, 10)
    z_i = dichotomy_step(g, x, y, z, w, de, st)
    print (f"g(x,y,z_i)={w}, z_i = {z_i} ; g(x,y,{z_i})={g(1,2,z_i)}")



def MAIN():
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
            GZs[j,i] = dichotomy_grad(g, x, y, 1, 0, 0.01, 10)
            print(f"{x} {y} done")
    
    plt.figure()
    plt.clf()
    plt.contourf(GXs, GYs, GZs, 50, cmap=plt.get_cmap("jet"))
    plt.colorbar()




# =============================================================================
# 
# =============================================================================
TEST_dichotomy()

#MAIN()








