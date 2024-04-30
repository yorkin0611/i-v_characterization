'''
I-V_analysis Version 3 is intended to make stuff as automatized as possible.
For alternative input/output check Version 2.
'''

# Some commonly used libraries.
#import csv, math, os, matplotlib, scipy
import numpy as np # Library required to read/load data from various file types.
from matplotlib import pyplot as plt #Library for plotting data.

# Libraries and functions for numerical analysis.
from scipy.interpolate import interp1d
from scipy.optimize import newton as nwt
from scipy.optimize import fminbound

# Libraries which provide file selection interface.
from tkinter import Tk
from tkinter.filedialog import askopenfile as getfile

# This function is required to extract just the input filename,
# not the whole directory path.
from os.path import basename as no_path

Tk().withdraw()
darkIV = getfile(title="Select dark IV measurement file:")  
lightIV = getfile(title="Select light IV measurement file:")

device_area = 0.04 # Define the device area in cm^2.

'''Define a delimiter (e.g. "\t" for tab, "," for csv), while skiprows indicates 
how many initial rows will not be loaded (e.g. name of individual columns).'''
fid1 = np.loadtxt(darkIV, delimiter="\t", skiprows=1)

V1 = fid1[0:, 0]
I1 = fid1[0:, 1]
J1 = I1 * 1000. / device_area 
# Converts the measured current I into current density J, based on defined device area.

# Prints out the loaded two-column array V-J.
#print("Dark J-V \n", np.vstack((V1,J1)).T, "\n")

fid2 = np.loadtxt(lightIV, delimiter="\t", skiprows=1)
V2 = fid2[0:, 0]
I2 = fid2[0:, 1]
J2 = I2 * 1000. / device_area

#print("Light J-V \n", np.vstack((V2,J2)).T, "\n")

 
# In principle, this pairs voltage V and current density J back into (x,y)
# pair to define a graph line from which Voc, Jsc, etc. can be found.
# darkspline = interp1d(V1, J1, kind="cubic", assume_sorted=False) 

Current_spline = interp1d(V2, J2, kind="cubic", assume_sorted=False) 
#Jsc_init = Current_spline(0)
Jsc = -Current_spline(0)

# Same thing, but coupling as (J,V) -> (x,y) pairs instead.
Voltage_spline = interp1d(J2, V2, kind="cubic", assume_sorted=False) 
Voc_init = Voltage_spline(0)
# It seems that code line 'Voc = Voltage_spline(0)' would also be enough
# to estimate Voc, but this Python code mimicks the original MATLAB script on purpose.


# === Estimating open-circuit voltage Voc and Jsc. ===
Voc = nwt(Current_spline, Voc_init, fprime=None, args=(), tol=1e-06, maxiter=2000, fprime2=None)
#Jsc = -nwt(Voltage_spline, Jsc_init, fprime=None, args=(), tol=1e-06, maxiter=2000, fprime2=None)


# === Estimating Vmax, Pmax & FF ===
Pwr2 = V2*J2

# (V, Pwr2) -> (x,y). Voltage-power curve to determine Pmax.
Pcspline = interp1d(V2, Pwr2, kind="cubic", assume_sorted=False)

# Because the curve is "upside down", we look for the "minimal" power between 0 and Voc.
Vmax = fminbound(Pcspline, 0, Voc)
Pmax = Vmax*Current_spline(Vmax)

FF = abs(Pmax/(Voc*Jsc))
PCE = Voc*Jsc*FF

print("Voc =", Voc, "V")
print("Jsc =", Jsc, "mA/cm2")
print("FF =", FF)
print("PCE =", PCE, "%")

plt.figure(1, dpi=200)
plt.plot(V1, abs(J1), marker="o", ls="", label="Dark I-V")
plt.plot(V2, abs(J2), marker="o", ls="", label="Light I-V")
plt.xlabel("Voltage V (V)")
plt.ylabel("Current density J (mA/cm$^2$)")
plt.yscale("log")
plt.title("I-V characterization")
plt.legend()
#plt.savefig(no_path(lightIV.name)+".png")

# plt.figure(1, dpi=200)
# plt.plot(V1, abs(J1), "k", marker="o", ls="", label="Dark I-V")
# plt.plot(inv_shockley_equation(Jred, n, J0, Rs), Jred, "r", label="Dark I-V Fit")
# plt.xlabel("Voltage V (V)")
# plt.ylabel("Current density J (mA/cm$^2$)")
# plt.yscale("log")
# plt.title("I-V characterization")
# plt.legend()

JVdark = "JV_"+no_path(darkIV.name)+".txt"
A1 = np.vstack((V1,J1)).T # or use A1 = np.column_stack((V1,J1))
#saving via numpy instruction for arrays
np.savetxt(JVdark, A1, delimiter="\t", header="V (V) \t J (mA/cm2)") 

JVlight = "JV_"+no_path(lightIV.name)+".txt"
A2 = np.vstack((V2,J2)).T # or use A2 = np.column_stack((V2,J2))
np.savetxt(JVlight, A2, delimiter="\t", header="V (V) \t J (mA/cm2)") 

# This block of code creates a text file in which parameters are written down.
fp3 = open("ParmPy_"+no_path(lightIV.name)+".txt", "w")
fp3.write(f"Voc = {Voc:.6f} V \n")
fp3.write(f"Jsc = {Jsc:.6f} mA/cm2 \n")
fp3.write(f"FF = {FF:.6f} \n")
fp3.write(f"PCE = {PCE:.6f} % \n")
fp3.close()

# ====== Currently unused stuff!!! ======
'''
fp1 = open(JVdark, "w") #saving via loop
fp1.write("Voltage (V) \t Current density (mA/cm2) \n")
for i in range (0, len(V1)):
    fp1.write(str(V1[i]) + "\t"+ str(J1[i]) + "\n")
fp1.close()

fp2 = open(JVlight, "w") #saving via loop
fp2.write("Voltage (V) \t Current density (mA/cm2) \n")
for i in range (0, len(V2)):
    fp2.write(str(V2[i]) + "\t"+ str(J2[i]) + "\n")
fp2.close()

#more graphs in one image
plt.figure(1, dpi=200)
plt.subplot(2,2,1)
#plt.title("I-V char")
plt.plot(V1, J1, 'k', label="Dark I-V")
plt.plot(V2, J2, 'r', label="Light I-V")
plt.legend()
#plt.xlabel("V (V)")
plt.ylabel("J (mA/cm$^2$)")
plt.axhline(y=0, color="black", linewidth=0.5)
plt.axvline(x=0, color="black", linewidth=0.5)
plt.minorticks_on()

plt.subplot(2,2,2)
#plt.title("I-V char (log-lin)")
plt.plot(V1, abs(J1), 'k', label="Dark I-V")
plt.plot(V2, abs(J2), 'r', label="Light I-V")
plt.yscale("log")
#plt.xlabel("V (V)")
#plt.ylabel("J (mA/cm$^2$)")
plt.minorticks_on()

plt.subplot(2,2,3)
#plt.plot(V1, abs(J1), label="Dark I-V")
#plt.title("P-V char")
plt.plot(V2, Pwr2, 'r', label="Light P-V")
plt.legend()
plt.xlabel("V (V)")
plt.ylabel("P (mW/cm$^2$)")
plt.minorticks_on()

plt.tight_layout()
plt.savefig(lightIV+"_many.png")

#zooming in Jsc-Voc region
plt.figure(2, dpi=300)
plt.plot(V1, J1, label="Dark I-V")
plt.plot(V2, J2, label="Light I-V")
plt.axis([-0.1, 1, -10, 1])
plt.grid(True)
plt.xlabel("Voltage V (V)")
plt.ylabel("Current density J (mA/cm$^2$)")
plt.title("I-V characterization")
plt.legend()
plt.savefig(no_path(lightIV.name)+"_zoomed.png")
'''