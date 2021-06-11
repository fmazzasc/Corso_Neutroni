import os
import re
import numpy as np
import matplotlib.pyplot as plt

files = os.listdir("result_arrays/")
h=10 #cm
for file in files:
    arr = np.load("result_arrays/" + file)
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.hist(arr, bins=100)
    res = re.findall(r"[-+]?\d*\.\d+|\d+", file)
    ax.text(0.1, 0.86, f" Thickness: {float(res[0])} cm \n Neutron Energy: {float(res[2])} MeV \n h/d: {h/float(res[1])}",  transform=ax.transAxes, fontsize=10)
    plt.xlabel("Proton Energy (MeV)", fontsize=12)
    plt.ylabel("Counts", fontsize=12)
    plt.savefig(f"result_plots/t_{float(res[0])}_E_{float(res[2])}_d_{float(res[1])}.png", bbox_inches="tight")
    plt.close()



