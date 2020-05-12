# -*- coding: utf-8 -*-
"""
Created on Tue May 12 22:36:50 2020

@author: Xavier
"""
import numpy as np


f = np.array(
    [["#","#","#"],
     ["#"," "," "],
     ["#","#","#"],
     ["#"," "," "],
     ["#"," "," "]])

for ii in range (0,2):
    print(np.flip(f, ii ))


