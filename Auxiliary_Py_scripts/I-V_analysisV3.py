# I-V_analysis Version 3 is intended to make stuff as automatized as possible.
# For alternative input/output check Version 2.

# Some commonly used libraries.
import csv, math, os 
# Library required to read/load data from various file types.
import numpy as np
#Library for plotting data.
from matplotlib import pyplot as plt

# Several libraries and functions for numerical analysis.
from scipy.interpolate import CubicSpline
from scipy.optimize import newton as nwt
from scipy.optimize import fminbound

# Libraries which provide file selection interface.
from tkinter import Tk
from tkinter.filedialog import askopenfile as getfile

# This function is required to extract just the input filename, not the whole directory path.
from os.path import basename as no_path

Tk().withdraw()

darkIV = getfile(title="Select dark IV measurement file:")  

lightIV = getfile(title="Select light IV measurement file:")

device_area = 0.04 # Define the device area in cm^2.

# Defines type of delimiter (e.g. "\t" for tab, "," for csv), while 'skiprows' indicates 
# how many initial rows will not be read/loaded (e.g. name of individual columns).
fid1 = np.loadtxt(darkIV, delimiter="\t", skiprows=1) 
V1 = fid1[0:, 0]   # Starting with 0th row, it reads everything in the 0th column.
I1 = fid1[0:, 1]   # Starting with 0th row, it reads everything in the 1st column.
J1 = I1 * 1000. / device_area # Converts the measured current I into current density J, based on defined device area. 

# In principle, this pairs voltage and current density J back into (x,y) coordinate pair 
# to define a graph line from which characteristic values (Voc, Jsc, etc.) can be extracted.
Vcspline1 = CubicSpline(V1, J1)
#print(cspline1)

fid2 = np.loadtxt(lightIV, delimiter="\t", skiprows=1)
V2 = fid2[0:, 0]
I2 = fid2[0:, 1]
J2 = I2 * 1000. / device_area
Vcspline2 = CubicSpline(V2, J2)
#print(cspline2)

# Estimating open-circuit voltage Voc. Package "scipy.optimize.newton" is used. 
Voc = nwt(Vcspline2, 0.5, fprime=None, args=(), tol=1e-06, maxiter=1000, fprime2=None)
print("Voc =", Voc, "V")

#Estimating short-circuit current Jsc. Here cubic spline is defined with (y,x) pairs instead.
Jcspline2 = CubicSpline(J2, V2)
Jsc = -nwt(Jcspline2, -5, fprime=None, args=(), tol=1e-06, maxiter=1000, fprime2=None)
print("Jsc =", Jsc, "mA/cm2")

#Estimating Vmax, Pmax & FF
Pwr2 = V2*J2
#print(Pwr2)

Pcspline = CubicSpline(V2, Pwr2)
Vmax = fminbound(Pcspline, 0, Voc)
#print("Vmax =", Vmax, "V")

Pmax = abs(Vmax*Vcspline2(Vmax))
#print("Pmax =", Pmax, "mW/cm2")

FF = Pmax/(Voc*Jsc)
print("FF =", FF)

PCE = Voc*Jsc*FF
print("PCE =", PCE, "%")

plt.plot(V1, J1, label="Dark I-V")
plt.plot(V2, J2, label="Light I-V")
plt.xlabel("Voltage V (V)")
plt.ylabel("Current density J (mA/cm$^2$)")
plt.title("I-V characterization")
plt.legend()
plt.savefig(no_path(lightIV.name)+".png")

JVdark = "JV_"+no_path(darkIV.name)+".txt"
A1 = np.vstack((V1,J1)).T # or use A1 = np.column_stack((V1,J1))
np.savetxt(JVdark, A1, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy

JVlight = "JV_"+no_path(lightIV.name)+".txt"
A2 = np.vstack((V2,J2)).T # or use A2 = np.column_stack((V2,J2))
np.savetxt(JVlight, A2, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy


# Writes a separate text file with Voc, Jsc, FF and PCE.
fp3 = open("Parm_"+no_path(lightIV.name)+".txt", "w")
#fp3.write("Voc = " + str(Voc) + " V \n")
fp3.write(f"Voc = {Voc:.6f} V \n") # Voc rounded up to six decimals.
fp3.write(f"Jsc = {Jsc:.6f} mA/cm2 \n")
fp3.write(f"FF = {FF:.6f} \n")
fp3.write(f"PCE = {PCE:.6f} %")
#fp3.write(f"\n\nVmax = {Vmax:.6f} V")
#fp3.write(f"\nPmax = {Pmax:.6f} mW/cm2")
fp3.close()