import numpy as np
import matplotlib.pyplot as plt
from helpers import *
import time

# a = np.array([3.75,1.7])
# print(a -3)

def find_closer_stopping_power(Ep_array, stopping_power):
    stopping_power_array = np.ones(len(Ep_array))

    energy_array = stopping_power[:,0]
    stopping_power = stopping_power[:,1]
    for index, Ep in enumerate(Ep_array): 
        diff_array = np.abs(energy_array  - Ep)
        min_index = np.argmin(diff_array)
        stopping_power_array[index] = stopping_power[min_index]
    return stopping_power_array

data = np.loadtxt("Stopping_power.txt")
closer = find_closer_stopping_power([3.75, 1.2], data)
print(closer)


# thick = 1e-2
# rho_pol = 0.93

# data = np.loadtxt("Stopping_power.txt")

# energy_in = data[:,0]
# energy = energy_in[energy_in>0]
# stopping_power = data[:,1][energy_in>0]
# energy_loss_exp = stopping_power*rho_pol*thick

# energy_loss_theory = compute_energy_loss(energy, thick, 0, 0)[0]
# fwhm_theory = compute_energy_loss(energy, thick, 0, 0)[1]

# for index in range(len(energy_loss_theory)):
#     energy_loss_theory[index] = landau_smearing(energy_loss_theory[index], fwhm_theory[index])

# plt.plot(energy, energy_loss_exp, "r-", label="exp")
# # plt.plot(energy, energy_loss_theory, "b-", label="theory")
# plt.legend()
# plt.savefig("en_loss.png")

