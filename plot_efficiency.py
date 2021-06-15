import os
import re
import numpy as np
import matplotlib.pyplot as plt

thickness = 1e-4

files = os.listdir("result_arrays/")
h=10 #cm

eff_arr_1 = []
eff_arr_2 = []
eff_arr_5 = []
d_arr_1 = []
d_arr_2 = []
d_arr_5 = []

for file in files:
    arr = np.load("result_arrays/" + file)
    res = re.findall(r"[-+]?\d*\.\d+|\d+", file)
    thick = float(res[0])
    d = float(res[1])
    E = float(res[2])
    if thick==thickness:
        if E==1:
            print("E: ", E, "d: ", d, "t: ", thick)
            eff_arr_1.append(len(arr)/(d*np.pi**2))
            d_arr_1.append(h/d)
        elif E==2:
            print("E: ", E, "d: ", d, "t: ", thick)
            eff_arr_2.append(len(arr)/(d*np.pi**2))
            d_arr_2.append(h/d)

        else:
            print("E: ", E, "d: ", d, "t: ", thick)
            eff_arr_5.append(len(arr)/(d*np.pi**2))
            d_arr_5.append(h/d)


d_arr_1 = np.array(d_arr_1)
d_arr_2 = np.array(d_arr_2)
d_arr_5 = np.array(d_arr_5)
eff_arr_1 = np.array(eff_arr_1)
eff_arr_2 = np.array(eff_arr_2)
eff_arr_5 = np.array(eff_arr_5)

index_1 = np.argsort(d_arr_1)
index_2 = np.argsort(d_arr_2)
index_5 = np.argsort(d_arr_5)

d_arr_1 = d_arr_1[index_1]
d_arr_2 = d_arr_2[index_2]
d_arr_5 = d_arr_5[index_5]

eff_arr_1 = eff_arr_1[index_1]
eff_arr_2 = eff_arr_2[index_2]
eff_arr_5 = eff_arr_5[index_5]




plt.plot(d_arr_1, eff_arr_1, ".-", label = "Neutron Energy = 1 MeV")
plt.plot(d_arr_2, eff_arr_2,  ".-", label = "Neutron Energy = 2 MeV")
plt.plot(d_arr_5, eff_arr_5,  ".-", label = "Neutron Energy = 5 MeV")
title_string = f"Thickness: {thickness} cm"
plt.legend(title=r'$\bf{{{}}}$'.format(title_string.replace(' ', r'\;')), fancybox=True)
plt.xlabel("h/d (cm)", fontsize=13)
plt.ylabel(r"# Protons / Collimator Area (cm$^{-2}$)", fontsize=13)

plt.savefig(f"result_plots/efficiency_t_{thickness}.png", bbox_inches="tight")