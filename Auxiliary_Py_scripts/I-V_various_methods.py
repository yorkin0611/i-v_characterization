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
#V1 = fid1[0:, 0]
#I1 = fid1[0:, 1]
V1 = np.sort(fid1[0:, 0])   # Starting with 0th row, it reads everything in the 0th column.
I1 = np.sort(fid1[0:, 1])   # Starting with 0th row, it reads everything in the 1st column.
J1 = I1 * 1000. / device_area # Converts the measured current I into current density J, based on defined device area.
DarkInputArray = np.vstack((V1,J1))
#print("Dark J-V \n", DarkInputArray, "\n")

fid2 = np.loadtxt(lightIV, delimiter="\t", skiprows=1)
#V2 = fid2[0:, 0]
#I2 = fid2[0:, 1]
V2 = np.sort(fid2[0:, 0])
I2 = np.sort(fid2[0:, 1])
J2 = I2 * 1000. / device_area
LightInputArray = np.vstack((V2,J2))
#print("Light J-V \n", LightInputArray, "\n")

# Group 1, 1d interpolation
V1_1dspline = interp1d(V1, J1, kind="cubic", assume_sorted=False) #1d interpolated dark
V2_1dspline = interp1d(V2, J2, kind="cubic", assume_sorted=False) #1d interpolated light (V,J)
J2_1dspline = interp1d(J2, V2, kind="cubic", assume_sorted=False) #1d interpolated light (J,V)

Voc_init1 = J2_1dspline(0)
Jsc_init1 = V2_1dspline(0)
print("Initial guesses by scipy-interp1d:")
print("Voc initial guess",Voc_init1, "V")
print("Jsc initial guess",Jsc_init1, "mA/cm2")
print()

# Group 2, cubic spline interpolation
Vcspline1 = CubicSpline(V1, J1) # In principle, this pairs voltage and current density J back into (x,y) coordinate pair to define a graph line from which characteristic values (Voc, Jsc, etc.) can be extracted.
Vcspline2 = CubicSpline(V2, J2) 
Jcspline2 = CubicSpline(J2, V2) # Same thing, but coupling as (y,x) pairs instead.

Voc_init2 = Jcspline2(0)
Jsc_init2 = Vcspline2(0)

print("Initial guesses by scipy-CubicSpline:")
print("Voc initial guess",Voc_init2, "V")
print("Jsc initial guess",Jsc_init2, "mA/cm2")
print()

# Difference in interpolations between Groups 1 and 2
#print("Difference between initial Voc-1 and Voc-2 is", ((Voc_init1-Voc_init2)/Voc_init2)*100,"%")
#print("Difference between initial Jsc-1 and Jsc-2 is", Jsc_init1-Jsc_init2,"mA/cm2")

# Estimating open-circuit voltage Voc and Jsc.
Voc1 = nwt(V2_1dspline, Voc_init1, fprime=None, args=(), tol=1e-08, maxiter=2000, fprime2=None)
Jsc1 = -nwt(J2_1dspline, Jsc_init1, fprime=None, args=(), tol=1e-08, maxiter=2000, fprime2=None)

Voc2 = nwt(Vcspline2, Voc_init2, fprime=None, args=(), tol=1e-08, maxiter=2000, fprime2=None)
Jsc2 = -nwt(Jcspline2, Jsc_init2, fprime=None, args=(), tol=1e-08, maxiter=2000, fprime2=None)

# print("Voc1 =", Voc1, "V")
# print("Voc2 =", Voc2, "V")

# print("Difference between estimated Voc-1 and Voc-2 is", ((Voc1-Voc2)/Voc2)*100,"%")
# print("Jsc1 =", Jsc1, "mA/cm2")

# Estimating Vmax, Pmax & FF
Pwr2 = V2*J2

P_1dspline = interp1d(V2, Pwr2, kind="cubic", assume_sorted=False)
Vmax1 = fminbound(P_1dspline, 0, Voc1)
Pmax1 = Vmax1*V2_1dspline(Vmax1)

FF1 = abs(Pmax1/(Voc1*Jsc1))
PCE1 = Voc1*Jsc1*FF1

print("Voc1 =", Voc1, "V")
print("Jsc1 =", Jsc1, "mA/cm2")
print("FF1 =", FF1)
print("PCE1 =", PCE1, "%")
print()

Pcspline = CubicSpline(V2, Pwr2)
Vmax2 = fminbound(Pcspline, 0, Voc2)
Pmax2 = Vmax2*Vcspline2(Vmax2)

FF2 = abs(Pmax2/(Voc2*Jsc2))
PCE2 = Voc2*Jsc2*FF2

print("Voc2 =", Voc2, "V")
print("Jsc2 =", Jsc2, "mA/cm2")
print("FF2 =", FF2)
print("PCE2 =", PCE2, "%")
print()

print("Difference between estimated Voc-1 and Voc-2 is", ((Voc1-Voc2)/Voc2)*100,"%")
print("Difference between estimated FF1 and FF2 is", ((FF1-FF2)/FF2)*100,"%")
print("Difference between estimated PCE1 and PCE2 is", ((PCE1-PCE2)/PCE2)*100,"%")

#-----------------------------------------------------------------#
# plt.plot(V2, J2, label="Light I-V scipy-interp1d")
# plt.plot(V2, J2, label="Light I-V scipy-CubicSpline")
# plt.xlabel("Voltage V (V)")
# plt.ylabel("Current density J (mA/cm$^2$)")
# plt.title("Comparison of various interpolation methods")
# plt.legend()
# plt.savefig("Interp_Comps.png")

# JVlight = "JV_"+lightIV+".txt"
# A2 = np.vstack((V2,J2)).T # or use A2 = np.column_stack((V2,J2))
# np.savetxt(JVlight, A2, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy instruction

# JVdark = "JV_"+darkIV+".txt"
# A1 = np.vstack((V1,J1)).T # or use A1 = np.column_stack((V1,J1))
# np.savetxt(JVdark, A1, delimiter="\t", header="V (V) \t J (mA/cm2)") #saving via numpy instruction

# fp3 = open("Parm_"+lightIV+".txt", "w") #saving via loop
# #fp3.write("Voc = " + str(Voc) + " V \n")
# fp3.write(f"Voc = {Voc:.6f} V \n") # Voc is first rounded up to six decimals, then converted to a string/letter.
# fp3.write(f"Jsc = {Jsc:.6f} mA/cm2 \n")
# fp3.write(f"FF = {FF:.6f} \n")
# fp3.write(f"PCE = {PCE:.6f} %")
# #fp3.write(f"\n\nVmax = {Vmax:.6f} V")
# #fp3.write(f"\nPmax = {Pmax:.6f} mW/cm2")
# fp3.close()

# #fp1 = open(JVdark, "w") #saving via loop
# #fp1.write("Voltage (V) \t Current density (mA/cm2) \n")
# #for i in range (0, len(V1)):
# #    fp1.write(str(V1[i]) + "\t"+ str(J1[i]) + "\n")
# #fp1.close()

# #fp2 = open(JVlight, "w") #saving via loop
# #fp2.write("Voltage (V) \t Current density (mA/cm2) \n")
# #for i in range (0, len(V2)):
# #    fp2.write(str(V2[i]) + "\t"+ str(J2[i]) + "\n")
# #fp2.close()