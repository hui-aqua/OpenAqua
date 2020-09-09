"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
"""

import numpy as np

import matplotlib.pyplot as plt


def jonswap(omega, hs, tp, gamma=3.3, gamma_auto=False):
    """
    Generate JONSWAP spectrum
    The Jonswap wave spectrum is expected to be a reasonable model for:
    3.6 < Tp/sqrt(hs) < 5
    # ref: DNVGL-RP-C205 ver.2018,pp22
    :param omega: numpy.ndarray | Array of frequencies
    :param hs: float  |  significant wave height [m]
    :param tp: float  |  peak wave period [s]
    :param gamma: float | peak shape parameter (default: 3.3)
    :param gamma_auto: Boolean (True or False)  |  The value of gamma will be calculated automatically
    :return: numpy.ndarray  |     Array of shape omega with wave energy densities
    """
    # default values
    sigma_low = 0.07
    sigma_high = 0.09
    if gamma_auto:
        if tp / np.sqrt(hs) <= 3.6:
            gamma = 5.0
        elif tp / np.sqrt(hs) >= 5:
            gamma = 1.0
        else:
            gamma = np.exp(5.75 - 1.15 * tp / np.sqrt(hs))

    # Pierson-Moskowitz
    omega_p = float(2 * np.pi / tp)
    pm_spectra = 5.0 / 16.0 * pow(hs, 2) * pow(omega_p, 4) * pow(omega, -5) * np.exp(-1.25 * pow(omega / omega_p, -4))

    # JONSWAP
    a_gamma = 1 - 0.287 * np.log(gamma)
    sigma = np.ones(omega.shape) * sigma_low
    sigma[omega > omega_p] = sigma_high
    jonswap_spectra = a_gamma * pm_spectra * pow(gamma, np.exp(-0.5 * pow((omega - omega_p) / (sigma * omega_p), 2)))
    return jonswap_spectra


if __name__ == "__main__":
    time_max = 600
    fre_range = np.linspace(2 * np.pi / time_max, 3, 1000)

    plt.figure()
    plt.plot(fre_range, jonswap(fre_range, 4, 8, gamma=5), label="jonswap_gama=5")
    plt.plot(fre_range, jonswap(fre_range, 4, 8, gamma=2), label="jonswap_gama=2")
    plt.plot(fre_range, jonswap(fre_range, 4, 8, gamma=1), label="jonswap_gama=1")
    plt.plot(fre_range, jonswap(fre_range, 4, 8, gamma=3.3), label="jonswap_gama=3.3")
    plt.xlabel("omega (rad)")
    plt.ylabel("S(omega)")
    plt.xlim(0, 3)
    plt.ylim(0, 6)
    plt.grid(True)
    plt.legend()
    plt.show()
