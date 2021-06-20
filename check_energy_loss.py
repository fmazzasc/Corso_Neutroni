import numpy as np
import matplotlib.pyplot as plt
from helpers import *
import time

# a = np.array([3.75,1.7])
# print(a -3)




thick = 1e-2
rho_pol = 0.93

data = np.loadtxt("Stopping_power.txt")

energy_in = data[:,0]
energy = energy_in[energy_in>0]
stopping_power = data[:,1][energy_in>0]
energy_loss_exp = stopping_power*rho_pol*thick

energy_loss_theory = compute_energy_loss(energy, thick, 0, 0)[0]
energy_loss_theory_landau = energy_loss_theory.copy()

fwhm_theory = compute_energy_loss(energy, thick, 0, 0)[1]

for index in range(len(energy_loss_theory)):
    energy_loss_theory_landau[index] = landau_smearing(energy_loss_theory_landau[index], fwhm_theory[index])

plt.plot(energy, energy_loss_exp, "r-", label="Experimental")
plt.plot(energy, energy_loss_theory_landau, "b--", label="Landau-Vavilov theory + smearing")
plt.plot(energy, energy_loss_theory, "g-", label="Landau-Vavilov theory")
plt.xlabel("Proton energy (MeV)", fontsize=12)
plt.ylabel(f"Energy loss in {thick} cm (Polyethilene)", fontsize=12)
plt.yscale("log")
plt.legend()
plt.savefig("result_plots/en_loss.png", bbox_inches="tight")

