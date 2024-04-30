#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 23 01:37:14 2021

@author: yorkin
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import derivative

# create 1000 equally spaced points between -10 and 10
x = np.linspace(-10, 10)

# calculate the y value for each element of the x vector
f = lambda x: 3*x**2 + 2*x - 1
y = f(x)
y1 = 6*x + 2

y_x = derivative(f, x)

plt.figure()
plt.plot(x, y, label="Parabola 3x$^2$+2x-1")
plt.plot(x, y_x, label="1st derivative of Parabola - diff(x)/diff(y)")
#plt.plot(x, y1, label="Linear function 6x+2")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Parabola and it's derivatives")
plt.legend()