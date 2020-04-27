# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 01:27:00 2020

@author: Xavier
"""
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

FIG = plt.figure(figsize=(2,2))
AX = FIG.add_subplot("111")
for i in range (1,11):
    TEXT_BOX = AnchoredText(f"{i}", loc=i, prop={'size': 8}, frameon=True)
    AX.add_artist(TEXT_BOX)
plt.title("AnchoredText locations")
plt.tight_layout()