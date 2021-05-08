import numpy as np
import matplotlib.pyplot as plt



n = 100000

t = 1e-2 ### cm
h = 100 ### cm
d = 10  ### cm

En = 1 #MeV


collimator_x_width = 20 


target_density = 0.93 * (6.022e23)/(28)    ###g/cm**3


cross_section_H_array = np.loadtxt("cross_sections/H.txt")
cross_section_C_array = np.loadtxt("cross_sections/C.txt")


def find_closer_cross_section(En, en_cross_section_array):
    energy_array = en_cross_section_array[:,0]*10**-6
    cross_section_array = en_cross_section_array[:,1]
    diff_array = np.abs(energy_array - En)
    min_index = np.argmin(diff_array)
    return cross_section_array[min_index]

def compute_mu(cross_section_H, cross_section_C, target_density):
    cross_section_H = cross_section_H *10**(-24)
    cross_section_C = cross_section_C *10**(-24)
    mu = target_density*4*cross_section_H + target_density*2*cross_section_C ## cm-1
    return mu


    


cross_section_H = find_closer_cross_section(En, cross_section_H_array)
cross_section_C = find_closer_cross_section(En, cross_section_C_array)
mu = compute_mu(cross_section_H, cross_section_C, target_density)

free_paths = -np.log(np.random.rand(n))/mu
x_coord = collimator_x_width*np.random.rand(n) - collimator_x_width/2
is_scattered = free_paths<=t


scattered_proton_y = free_paths[is_scattered]
scattered_proton_x = x_coord[is_scattered]


