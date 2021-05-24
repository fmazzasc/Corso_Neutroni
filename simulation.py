import numpy as np
import matplotlib.pyplot as plt
from helpers import *



n = int(1e8)
t = 1e-2### cm
h = 100 ### cm
d = 10  ### cm
En = 1 #MeV
collimator_x_width = 20 
target_density = 0.93 * (6.022e23)/(28)    ###g/cm**3
do_landau_smearing = True
detector_smearing = 0.02

print("d/h: ", d/h)
print("Thickness: ", t)
print("N neutrons generated: ", n)
print("Neutron energy: ", En)
print("Detector energy resolution: ", detector_smearing)

print("--------------------------")


cross_section_H_array = np.loadtxt("cross_sections/H.txt")
cross_section_C_array = np.loadtxt("cross_sections/C.txt")
cross_section_H = find_closer_cross_section(En, cross_section_H_array)
cross_section_C = find_closer_cross_section(En, cross_section_C_array)
mu = compute_mu(cross_section_H, cross_section_C, target_density)
free_paths = -np.log(np.random.rand(n))/mu
print("Mean free path: ", np.mean(free_paths))


x_coord = collimator_x_width*np.random.rand(n) - collimator_x_width/2
is_scattered = free_paths<=t
scattered_proton_x = x_coord[is_scattered]
scattered_proton_y = free_paths[is_scattered]
n_protons = len(scattered_proton_x)
print("N protons converted: ", n_protons)
print("--------------------------")

prob_h_scat = 4*cross_section_H/(4*cross_section_H + 2*cross_section_C)
prob_c_scat = 1-prob_h_scat
is_h_scattering = np.random.rand(n_protons) <= prob_h_scat
scattered_proton_theta = compute_theta(n_protons, is_h_scattering)
scattered_proton_energy = En*np.cos(scattered_proton_theta)**2

is_proton_detected = np.zeros(n_protons, dtype=bool)
for index in range(n_protons):
    scat_x = scattered_proton_x[index]
    scat_y = scattered_proton_y[index]
    scat_theta = scattered_proton_theta[index]
    if scat_x < -d/2 -h/d*t:
        continue
    if scat_x > d/2 + h/d*t:
        continue
    if -d/2 -h/d*t<=scat_x<-d/2:
        if (-scat_x-d/2)/(t-scat_y) <= np.tan(scat_theta) <= (d/2 - scat_x)/(h+t-scat_y):
            is_proton_detected[index] = 1
        continue
    if -d/2<=scat_x<d/2:
        if -((d/2 + scat_x)/(h+t-scat_y))<np.tan(scat_theta)<(d/2 - scat_x)/(h + t - scat_y):
            is_proton_detected[index] = 1
        continue
    if d/2 <= scat_x <= d/2 + h/d*t:
        if -(d/2 + scat_x)/(h+t-scat_y)<np.tan(scat_theta)<-(scat_x - d/2)/(t - scat_y):
            is_proton_detected[index] = 1
        continue

print("Max tan: ", np.max(np.tan(scattered_proton_theta[is_proton_detected])))
print("N protons detected (before energy loss correction): ", len(scattered_proton_theta[is_proton_detected]))
print("--------------------------")

mpv_energy_loss = compute_energy_loss(scattered_proton_energy[is_proton_detected], t, scattered_proton_y[is_proton_detected], scattered_proton_theta[is_proton_detected])[0]
fwhm_energy_loss = compute_energy_loss(scattered_proton_energy[is_proton_detected], t, scattered_proton_y[is_proton_detected], scattered_proton_theta[is_proton_detected])[0]

energy_loss = mpv_energy_loss
if do_landau_smearing:
    for index in range(len(mpv_energy_loss)):
        energy_loss[index] = landau_smearing(mpv_energy_loss[index], fwhm_energy_loss[index])
print("Mean energy loss: ", np.mean(energy_loss))

proton_final_energy = scattered_proton_energy[is_proton_detected] - energy_loss
proton_final_energy = proton_final_energy[proton_final_energy>0]

if detector_smearing > 0:
    proton_final_energy = np.random.normal(loc=proton_final_energy, scale=detector_smearing)

print("N protons detected (after energy loss correction): ", len(proton_final_energy))
print("--------------------------")

plt.hist(proton_final_energy, bins=30)
plt.xlabel("Proton energy (MeV)")
plt.ylabel("Counts")
plt.show()
