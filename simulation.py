import numpy as np
import matplotlib.pyplot as plt
from helpers import *
import time


start_time = time.time()

test_mode = False

n = int(1e9) if not test_mode else int(1e9)
t_list = [1e-4, 1e-3, 1e-2, 1e-1] if not test_mode else [1e-2] #cm
h = 10 ### cm  
d_list = [1, 0.1, 0.01] if not test_mode else [1] ### cm 
En_list = [1, 2, 5] if not test_mode else [2] #MeV

target_density = 0.93 * (6.022e23)/(28)    ###g/cm**3

do_energy_loss = True
do_exp_energy_loss = True
do_landau_smearing = False if do_exp_energy_loss else True
detector_smearing = 0.01

counter = 1
for t in t_list:
    for d in d_list:
        for En in En_list:
            print(f"Combination: {counter}/{36}")
            counter += 1
            print("--------------------------")
            print("h/d: ", h/d)
            print("Thickness(cm): ", t)
            print("N neutrons generated: ", n)
            print("Neutron energy(MeV): ", En)
            print("Detector energy resolution(MeV): ", detector_smearing)
            print("--------------------------")
            collimator_x_width = 2*(d/2 + d/h*t)
            cross_section_H_array = np.loadtxt("cross_sections/H.txt")
            cross_section_C_array = np.loadtxt("cross_sections/C.txt")
            cross_section_H = find_closer_cross_section(En, cross_section_H_array)
            cross_section_C = find_closer_cross_section(En, cross_section_C_array)
            cross_section_C = 0
            mu = compute_mu(cross_section_H, cross_section_C, target_density)
            free_paths = -np.log(np.random.rand(n))/mu
            print("Mean free path (cm): ", np.mean(free_paths))

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
            scattered_proton_theta, scattered_proton_energy = compute_theta_and_energy(n_protons, En, is_h_scattering)
            # plt.figure()
            # plt.hist(scattered_proton_theta, bins=100)
            # plt.savefig("theta.png")
            # plt.close()
            # plt.figure()
            # plt.hist(scattered_proton_energy, bins=100)
            # plt.savefig("en.png")
            # plt.close()


            is_proton_detected = np.zeros(n_protons, dtype=bool)
            for index in range(n_protons):
                scat_x = scattered_proton_x[index]
                scat_y = scattered_proton_y[index]
                scat_theta = scattered_proton_theta[index]
                if scat_x < -d/2 -d/h*t:
                    continue
                if scat_x > d/2 + d/h*t:
                    continue
                if -d/2 -d/h*t<=scat_x<-d/2:
                    if (-scat_x-d/2)/(t-scat_y) <= np.tan(scat_theta) <= (d/2 - scat_x)/(h+t-scat_y):
                        is_proton_detected[index] = 1
                    continue
                if -d/2<=scat_x<d/2:
                    if -((d/2 + scat_x)/(h+t-scat_y))<np.tan(scat_theta)<(d/2 - scat_x)/(h + t - scat_y):
                        is_proton_detected[index] = 1
                    continue
                if d/2 <= scat_x <= d/2 + d/h*t:
                    if -(d/2 + scat_x)/(h+t-scat_y)<np.tan(scat_theta)<-(scat_x - d/2)/(t - scat_y):
                        is_proton_detected[index] = 1
                    continue
            if np.sum(is_proton_detected)==0:
                print("No protons detected, moving to next run")
                continue
            print("Max tan: ", np.max(np.tan(scattered_proton_theta[is_proton_detected])))
            print("N protons detected (before energy loss correction): ", len(scattered_proton_theta[is_proton_detected]))
            print("--------------------------")

            if do_energy_loss==True:
                if do_exp_energy_loss==False:
                    mpv_energy_loss = compute_energy_loss(scattered_proton_energy[is_proton_detected], t, scattered_proton_y[is_proton_detected], scattered_proton_theta[is_proton_detected])[0]
                    fwhm_energy_loss = compute_energy_loss(scattered_proton_energy[is_proton_detected], t, scattered_proton_y[is_proton_detected], scattered_proton_theta[is_proton_detected])[1]
                else:
                    proton_path = (t-scattered_proton_y[is_proton_detected])/np.abs(np.cos(scattered_proton_theta[is_proton_detected]))
                    data = np.loadtxt("Stopping_power.txt")
                    stopping_power = find_closer_stopping_power(scattered_proton_energy[is_proton_detected], data)
                    mpv_energy_loss = stopping_power*0.93*proton_path

            else:
                mpv_energy_loss = np.zeros(len(scattered_proton_energy[is_proton_detected]))
                fwhm_energy_loss = np.zeros(len(scattered_proton_energy[is_proton_detected]))
                do_landau_smearing = False

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
            arr_name = f"thick_{t}_d_{d}_En_{En}"
            np.save("result_arrays/" + arr_name, proton_final_energy)


stop_time = time.time()
print(f"Final time: {stop_time - start_time}")