# Plot Landau distributions with FWHM and MPV
import numpy as np
import matplotlib.pyplot as plt
import pylandau


def landau_smearing(mpv, fwhm):  # http://stackoverflow.com/questions/10582795/finding-the-full-width-half-maximum-of-a-peak
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

mpv = 0.1
fwhm = 0.08
values = []
for i in range(int(1e4)):
    values.append(landau_smearing(mpv, fwhm))

# print(values)
plt.hist(values, bins = 50)
plt.show()
