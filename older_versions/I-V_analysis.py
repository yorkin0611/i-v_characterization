import csv, math, os
import numpy as np

from os.path import basename as filename

from matplotlib import pyplot as plt

from scipy.interpolate import CubicSpline
from scipy.optimize import newton as nwt
from scipy.optimize import fminbound

from tkinter import Tk
from tkinter.filedialog import askopenfile as getfile

Tk().withdraw()
#darkIV = getfile(title="Select dark IV file")  # "DarkIV-1" Insert the name of dark I-V curve file, including the extenstion, if there is one.
darkIV = "DarkIV-1"
lightIV = "LightIV-1"  # Insert the name of light I-V curve file, including the extenstion, if there is one.

device_area = 0.04     # Define the device area in cm^2.

fid1 = np.loadtxt(darkIV, delimiter="\t", skiprows=1)   # Defines type of delimiter (e.g. "\t" for tab, "," for csv), while skiprows indicates how many initial rows will not be loaded (e.g. name of individual columns).
V1 = fid1[0:, 0]   # Starting with 0th row, it reads everything in the 0th column.
I1 = fid1[0:, 1]   # Starting with 0th row, it reads everything in the 1st column.
J1 = I1 * 1000. / device_area # Converts the measured current I into current density J, based on defined device area. 
#Vcspline1 = CubicSpline(V1, J1)
#print(cspline1)

fid2 = np.loadtxt(lightIV, delimiter="\t", skiprows=1)
V2 = fid2[0:, 0]
I2 = fid2[0:, 1]
J2 = I2 * 1000. / device_area
Vcspline2 = CubicSpline(V2, J2) # In principle, this pairs voltage and current density J back into (x,y) coordinate pair to define a graph line from which characteristic values (Voc, Jsc, etc.) can be extracted.
#print(cspline2)

# Estimating open-circuit voltage Voc.
# Package "scipy.optimize.newton" is used  
Voc = nwt(Vcspline2, 0.5, fprime=None, args=(), tol=1e-06, maxiter=1000, fprime2=None)
print("Voc =", Voc, "V")

#Estimating short-circuit current Jsc
Jcspline2 = CubicSpline(J2, V2)
Jsc = -nwt(Jcspline2, -5, fprime=None, args=(), tol=1e-06, maxiter=1000, fprime2=None)
print("Jsc =", Jsc, "mA/cm2")

#Estimating Vmax, Pmax & FF
Pwr2 = V2*J2
#print(Pwr2)

Pcspline = CubicSpline(V2, Pwr2)
Vmax = fminbound(Pcspline, 0, Voc)
#print("Vmax =", Vmax, "V")
Pmax = Vmax*Vcspline2(Vmax)

FF = abs(Pmax/(Voc*Jsc))
print("FF =", FF)

PCE = Voc*Jsc*FF
print("PCE =", PCE, "%")

plt.plot(V1, J1, label="Dark I-V")
plt.plot(V2, J2, label="Light I-V")
plt.xlabel("Voltage V (V)")
plt.ylabel("Current density J (mA/cm$^2$)")
plt.title("I-V characterization")
plt.legend()
plt.savefig(lightIV+".png")

JVlight = "JV_"+lightIV+".txt"
A2 = np.vstack((V2,J2)).T # or use A2 = np.column_stack((V2,J2))
np.savetxt(JVlight, A2, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy instruction

JVdark = "JV_"+filename(darkIV.name)+".txt"
A1 = np.vstack((V1,J1)).T # or use A1 = np.column_stack((V1,J1))
np.savetxt(JVdark, A1, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy instruction
print(JVdark)

fp3 = open("Parm_"+lightIV+".txt", "w") #saving via loop
#fp3.write("Voc = " + str(Voc) + " V \n")
fp3.write(f"Voc = {Voc:.6f} V \n") # Voc is first rounded up to six decimals, then converted to a string/letter.
fp3.write(f"Jsc = {Jsc:.6f} mA/cm2 \n")
fp3.write(f"FF = {FF:.6f} \n")
fp3.write(f"PCE = {PCE:.6f} %")
#fp3.write(f"\n\nVmax = {Vmax:.6f} V")
#fp3.write(f"\nPmax = {Pmax:.6f} mW/cm2")
fp3.close()

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