# Dark I-V Version 1 should work on explicit version of diode equation,
# in which voltage V is function of current density J, i.e. V(J).

# Some commonly used libraries.
#import csv, math, os, matplotlib, scipy
import numpy as np # Library required to read/load data from various file types.
from matplotlib import pyplot as plt #Library for plotting data.

# Several libraries and functions for numerical analysis.
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Libraries which provide file selection interface.
from tkinter import Tk
from tkinter.filedialog import askopenfile as getfile

# This function is required to extract just the input filename, not the whole directory path.
#from os.path import basename as no_path

def inv_shockley_equation(Jred, n, J0, Rs):
    return (n*kB*T)/e * (np.log(Jred/J0 + 1))  + Jred*Rs

Tk().withdraw()
darkIV = getfile(title="Select dark IV measurement file:")  
#darkIV = "Dev1P1_dark"  

V_limit = 0.8
device_area = 0.04 # Define the device area in cm^2.

e = 1.602176e-19 # Electric charge
kB = 1.380649e-23 # Boltzmann constant
T = 300.0 # Temperature in Kelvin

fid1 = np.loadtxt(darkIV, delimiter="\t", skiprows=1)
# Define a delimiter (e.g. "\t" for tab, "," for csv), while skiprows indicates 
# how many initial rows will not be loaded (e.g. name of individual columns).

V1 = fid1[0:, 0]
I1 = fid1[0:, 1]
J1 = I1 * 1000. / device_area
#Jabs = abs(J1)

darkspline = interp1d(V1, J1, kind="cubic", assume_sorted=False) 
# In principle, this pairs voltage and current density J back into (x,y) coordinate pair 
# to define a graph line from which characteristic values (Voc, Jsc, etc.) can be extracted.

Vred = V1[V1 >= V_limit]
Jred = (darkspline(Vred))

init_vals = [1.2, 1e-7, 1e-3]

fits, covar = curve_fit(inv_shockley_equation, Jred, Vred, p0 = init_vals)

n, J0, Rs = fits[0], fits[1], fits[2]

print("Dark J-V \n", np.vstack((Vred,Jred)).T, "\n")
print("Ideality factor n =", n)
print("Saturation current J0 =", J0, "mA/cm2")
print("Series resistance Rs =", Rs, "m-Ohms/cm2")

plt.figure(1, dpi=200)
plt.plot(V1, abs(J1), "k", marker="o", ls="", label="Dark I-V")
plt.plot(inv_shockley_equation(Jred, n, J0, Rs), Jred, "r", label="Dark I-V Fit")
plt.xlabel("Voltage V (V)")
plt.ylabel("Current density J (mA/cm$^2$)")
plt.yscale("log")
plt.title("I-V characterization")
plt.legend()

#plt.show()
#plt.savefig(no_path(lightIV.name)+".png")

#JVdark = "JV_"+no_path(darkIV.name)+".txt"
#A1 = np.vstack((V1,J1)).T # or use A1 = np.column_stack((V1,J1))
#np.savetxt(JVdark, A1, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy instruction

# ====== Currently unused stuff!!! ======
#fp1 = open(JVdark, "w") #saving via loop
#fp1.write("Voltage (V) \t Current density (mA/cm2) \n")
#for i in range (0, len(V1)):
#    fp1.write(str(V1[i]) + "\t"+ str(J1[i]) + "\n")
#fp1.close()

#fp2 = open(JVlight, "w") #saving via loop
#fp2.write("Voltage (V) \t Current density (mA/cm2) \n")
#for i in range (0, len(V2)):
#    fp2.write(str(V2[i]) + "\t"+ str(J2[i]) + "\n")
#fp2.close()

# plt.figure(2, dpi=300)
# plt.plot(V1, J1, label="Dark I-V")
# plt.plot(V2, J2, label="Light I-V")
# plt.axis([-0.1, 1, -10, 1])
# plt.grid(True)
# plt.xlabel("Voltage V (V)")
# plt.ylabel("Current density J (mA/cm$^2$)")
# plt.title("I-V characterization")
# plt.legend()
# plt.savefig(no_path(lightIV.name)+"_zoomed.png")