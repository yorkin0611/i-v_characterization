import csv, math, os
import numpy as np

from os.path import basename as filename

from matplotlib import pyplot as plt

from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from scipy.optimize import newton as nwt
from scipy.optimize import fminbound

from tkinter import Tk
from tkinter.filedialog import askopenfile as getfile

darkIV = "Dev1P1_dark"
lightIV = "Dev1P1_light"  # Insert the name of light I-V curve file, including the extenstion, if there is one.

device_area = 0.04     # Define the device area in cm^2.

fid1 = np.loadtxt(darkIV, delimiter="\t", skiprows=1)   # Defines type of delimiter (e.g. "\t" for tab, "," for csv), while skiprows indicates how many initial rows will not be loaded (e.g. name of individual columns).
V1 = fid1[0:, 0]
I1 = fid1[0:, 1]
#V1 = np.sort(fid1[0:, 0])   # Starting with 0th row, it reads everything in the 0th column.
#I1 = np.sort(fid1[0:, 1])   # Starting with 0th row, it reads everything in the 1st column.
J1 = I1 * 1000. / device_area # Converts the measured current I into current density J, based on defined device area.
#print("Dark J-V \n", np.vstack((V1,J1)).T, "\n")

fid2 = np.loadtxt(lightIV, delimiter="\t", skiprows=1)
V2 = fid2[0:, 0]
I2 = fid2[0:, 1]
#V2 = np.sort(fid2[0:, 0])
#I2 = np.sort(fid2[0:, 1])
J2 = I2 * 1000. / device_area
#print("Light J-V \n", np.vstack((V2,J2)).T, "\n")

#darkspline = interp1d(V1, J1, kind="cubic", assume_sorted=False) # In principle, this pairs voltage and current density J back into (x,y) coordinate pair to define a graph line from which characteristic values (Voc, Jsc, etc.) can be extracted.
Current_spline = interp1d(V2, J2, kind="cubic", assume_sorted=False) 
#Jsc_init = Current_spline(0)
Jsc = -Current_spline(0)

Voltage_spline = interp1d(J2, V2, kind="cubic", assume_sorted=False) # Same thing, but coupling as (y,x) pairs instead.
Voc_init = Voltage_spline(0)
#Voc = Voltage_spline(0)

#print("Voc initial guess",Voc_init, "V")
#print("Jsc initial guess",Jsc_init, "mA/cm2")
#print()

# Estimating open-circuit voltage Voc and Jsc.
Voc = nwt(Current_spline, Voc_init, fprime=None, args=(), tol=1e-06, maxiter=2000, fprime2=None)
#Jsc = -nwt(Voltage_spline, Jsc_init, fprime=None, args=(), tol=1e-06, maxiter=2000, fprime2=None)
#print("Fractional difference between Voc-init and Voc:", ((Voc_init-Voc)/Voc)*100,"%")

#Estimating Vmax, Pmax & FF
Pwr2 = V2*J2
Pcspline = interp1d(V2, Pwr2, kind="cubic", assume_sorted=False)
Vmax = fminbound(Pcspline, 0, Voc)
Pmax = Vmax*Current_spline(Vmax)

FF = abs(Pmax/(Voc*Jsc))
PCE = Voc*Jsc*FF

print("Voc =", Voc, "V")
print("Jsc =", Jsc, "mA/cm2")
print("FF =", FF)
print("PCE =", PCE, "%")
#print("Vmax =", Vmax, "V")
#print("Pmax =", Pmax, "V")

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

JVdark = "JV_"+darkIV+".txt"
A1 = np.vstack((V1,J1)).T # or use A1 = np.column_stack((V1,J1))
np.savetxt(JVdark, A1, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy instruction

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