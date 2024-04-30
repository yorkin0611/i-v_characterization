#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 16:29:43 2021

@author: yorkin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

xdata = np.linspace(0, 4, 50)
y = func(xdata, 2.5, 1.3, 0.5)
ydata = y + 0.2 * np.random.normal(size=len(xdata))

popt, pcov = curve_fit(func, xdata, ydata,p0=[2,1,1])

plt.ion()
plt.plot(xdata,ydata,'o')
xplot = np.linspace(0,4,100)
plt.plot(xplot,func(xplot,*popt))