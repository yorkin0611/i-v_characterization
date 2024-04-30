#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 23 02:08:33 2021

@author: yorkin
"""

# importing the library
import matplotlib.pyplot as plt
from scipy.misc import derivative
import numpy as np
  
# defining the function
def function(x):
    return 2*x*x*x+x+3
  
# calculating its derivative
def deriv(x):
    return derivative(function, x)
  
# defininf x-axis intervals
y = np.linspace(-6, 6)
  
# plotting the function
plt.plot(y, function(y), color='purple', label='Function')
  
# plotting its derivative
plt.plot(y, deriv(y), color='green', label='Derivative')
  
# formatting
plt.legend(loc='upper left')
plt.grid(True)