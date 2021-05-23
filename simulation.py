import numpy as np
import matplotlib.pyplot as plt



n = int(1e8)

t = 1e-2### cm
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


def compute_theta(n_protons, is_h_scattered):
    phi = 2*np.pi*np.random.rand(n_protons)
    is_positive = np.random.rand(n_protons)
    theta = np.ones(n_protons)
    for index in range(len(theta)):
        A=1
        if is_h_scattered[index]==False:
            A=12
        
        theta[index] = np.arccos((A*np.cos(phi[index]) + 1)/np.sqrt(A**2 + 2*A*np.cos(phi[index]) + 1))
        theta[index] = theta[index] if is_positive[index]<0.5 else -theta[index]
    return theta

def compute_energy_loss(scattered_proton_energy, thickness, scat_y, theta):
    sqr = lambda x: x*x
    rho = 0.93 #g/cm**3
    Z = (4*1 + 6*2)/3
    A = (4*1 + 2*12)/3
    proton_path = (thickness-scat_y)/np.abs(np.cos(theta)) #cm

    M = 939 # Mass of heavy particle in MeV
    m_e = 0.511 # Mass of electron in MeV
    K = 0.307075 # constant K in MeV cm mol^-1
    p_proton = np.sqrt(sqr(scattered_proton_energy) + scattered_proton_energy*M*2)
    epsilon = (K * rho * Z * proton_path) / (2 * A)
    gamma = lambda p: np.sqrt(1 + sqr(p / M))
    beta = lambda p: np.sqrt(1 - 1 / sqr(gamma(p)))
    beta_gamma = lambda p: p / M
    I = 0.000057 # Polyethilene Ionisation energy in MeV

    E_mpv = lambda p:  epsilon * (np.log( ( 2 * m_e * sqr( beta_gamma(p) ) ) / I ) + np.log(epsilon / (I * sqr(beta(p)))) +  0.2 - sqr(beta(p)))/sqr(beta(p))
    fwhm = 4*epsilon/sqr(beta(p_proton))
    mpv = E_mpv(p_proton)

    return mpv, fwhm


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

print("N protons detected: ", len(scattered_proton_theta[is_proton_detected]))
print("d/h: ", d/h)
# print("theta: ", scattered_proton_theta[is_proton_detected])
print("Max tan: ", np.max(np.tan(scattered_proton_theta[is_proton_detected])))

energy_loss = compute_energy_loss(scattered_proton_energy[is_proton_detected], t, scattered_proton_y[is_proton_detected], scattered_proton_theta[is_proton_detected])[0]
print("Mean energy loss: ", np.mean(energy_loss))

proton_final_energy = scattered_proton_energy[is_proton_detected] - energy_loss
proton_final_energy = proton_final_energy[proton_final_energy>0]
plt.hist(proton_final_energy)
plt.xlabel("Proton energy (MeV)")
plt.ylabel("Counts")
plt.show()



