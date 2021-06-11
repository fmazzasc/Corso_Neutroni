import numpy as np
import matplotlib.pyplot as plt
import pylandau



def find_closer_cross_section(En, en_cross_section_array):
    energy_array = en_cross_section_array[:,0]*10**-6
    cross_section_array = en_cross_section_array[:,1]
    diff_array = np.abs(energy_array - En)
    min_index = np.argmin(diff_array)
    return cross_section_array[min_index]

def find_closer_stopping_power(Ep_array, stopping_power):
    stopping_power_array = np.ones(len(Ep_array))
    energy_array = stopping_power[:,0]
    stopping_power = stopping_power[:,1]
    for index, Ep in enumerate(Ep_array): 
        diff_array = np.abs(energy_array  - Ep)
        min_index = np.argmin(diff_array)
        stopping_power_array[index] = stopping_power[min_index]
    return stopping_power_array

def compute_mu(cross_section_H, cross_section_C, target_density):
    cross_section_H = cross_section_H *10**(-24)
    cross_section_C = cross_section_C *10**(-24)
    mu = target_density*4*cross_section_H + target_density*2*cross_section_C ## cm-1
    return mu


def compute_theta_and_energy(n_protons, E_neutr, is_h_scattered):
    cos_phi = 2*np.random.rand(n_protons) -1
    is_positive = np.random.rand(n_protons)
    theta = np.ones(n_protons)
    cos_th = np.ones(n_protons)
    energy = np.ones(n_protons)
    for index in range(len(theta)):
        A = 1
        if is_h_scattered[index]==False:
            A=12
        cos_theta = (A*cos_phi[index] + 1)/np.sqrt(A**2 + 2*A*cos_phi[index] + 1)
        cos_th[index] = cos_theta
        theta[index] = np.arccos(cos_theta)
        theta[index] = theta[index] if is_positive[index]<0.5 else -theta[index]
        energy[index] = 4*A/((1+A)**2)*cos_theta**2*E_neutr
    return theta, energy


def compute_energy_and_theta(n_protons, E_neutr, is_h_scattered):
    is_positive = np.random.rand(n_protons)
    theta = np.ones(n_protons)
    energy = np.ones(n_protons)
    for index in range(len(theta)):
        # is_h_scattered[index] = False
        A=1
        if is_h_scattered[index]==False:
            A=12
        energy[index] = np.random.rand()*E_neutr*(4*A/(1+A)**2)
        theta[index] = np.arccos(np.sqrt(energy[index]/(4*A/(1+A)**2)/E_neutr))
        # theta[index] = theta[index] if is_positive[index]<0.5 else -theta[index]
    return theta, energy



def compute_energy_loss(scattered_proton_energy, thickness, scat_y, theta):
    sqr = lambda x: x*x
    rho = 0.93 #g/cm**3
    Z = (4*1 + 6*2)/10
    A = (4*1 + 2*12)/10
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

def landau_smearing(mpv, fwhm):
    eta_in = fwhm/4
    if eta_in < 1:
        expo = np.abs(np.floor(np.log10(eta_in)))
        mpv = mpv*10**expo
        fwhm = fwhm*10**expo
        eta = eta_in*10**expo
    else:
        eta=eta_in
    lower_x = np.max(mpv - fwhm, 0)
    upper_x = mpv + fwhm*2
    is_smeared = 0
    count = 0
    # print(mpv, eta)
    while is_smeared==0 and count<100:
        x = np.random.rand(1)*(upper_x - lower_x) + lower_x
        y_land = pylandau.landau(x, mpv=mpv, eta=eta, A=1)
        y = np.random.rand(1)
        if y<y_land:
            is_smeared=1
        count+=1
    if eta_in<1:
        x /=10**expo
    return float(x)
