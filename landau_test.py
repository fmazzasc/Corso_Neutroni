
import numpy as np
import matplotlib.pyplot as plt
import ROOT



En = 10
scat_y = 0.12
theta = 0


def compute_energy_loss(scattered_proton_energy, scat_y, theta):
    sqr = lambda x: x*x
    rho = 0.93 #g/cm**3
    Z = (4*1 + 6*2)/3
    A = (4*1 + 2*12)/3
    proton_path = (1-scat_y)/np.cos(theta) #cm
    print("path:", proton_path)
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

mpv,fwhm = compute_energy_loss(En, scat_y, theta)
print(mpv, fwhm)
land_distr = ROOT.TF1("fa3","TMath::Landau(x,[0],[1])", 0, 200)
land_distr.SetParameters(mpv, fwhm)

en_list = []
i=0
while i<1e5:
    en_list.append(land_distr.GetRandom())
    i+=1

# print(en_list)
plt.hist(en_list, bins=500)
plt.show()


# E_MPV_LVB = lambda p:  epsilon * ( log( ( 2 * m_e * sqr( beta_gamma(p) ) ) / I ) + log( epsilon / ( I * sqr(beta(p)) ) ) +  0.2 - sqr( beta(p) ) - delta_full( X(p) ) )  / sqr(beta(p)) 