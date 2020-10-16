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
from numpy import pi
import random


def pierson_moskowitz_spectra(omega, hs, tp):
    """
    Generate Pierson-Moskowitz wave specturm \n
    ref: DNVGL-RP-C205 ver.2018,pp64
    :param omega: numpy.ndarray | Array of frequencies
    :param hs: float  |  significant wave height [m]
    :param tp: float  |  peak wave period [s]
    :return: numpy.ndarray  |     Array of shape omega with wave energy densities
    """
    omega_p = float(2 * pi / tp)
    spectra = 5.0 / 16.0 * pow(hs, 2) * pow(omega_p, 4) * pow(omega, -5) * np.exp(-1.25 * pow(omega / omega_p, -4))
    return spectra


def jonswap_spectra(omega, hs, tp, gamma=3.3, gamma_auto=False):
    """
    Generate JONSWAP spectrum
    The Jonswap wave spectrum is expected to be a reasonable model for:
    3.6 < Tp/sqrt(hs) < 5 \n
    ref: DNVGL-RP-C205 ver.2018,pp64
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
    omega_p = float(2 * pi / tp)
    # pm_spectra = 5.0 / 16.0 * pow(hs, 2) * pow(omega_p, 4) * pow(omega, -5) * np.exp(-1.25 * pow(omega / omega_p, -4))
    pm_spectra = pierson_moskowitz_spectra(omega, hs, tp)
    # JONSWAP
    a_gamma = 1 - 0.287 * np.log(gamma)
    sigma = np.ones(omega.shape) * sigma_low
    sigma[omega > omega_p] = sigma_high
    spectra = a_gamma * pm_spectra * pow(gamma, np.exp(-0.5 * pow((omega - omega_p) / (sigma * omega_p), 2)))
    return spectra


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    # # Plot1  wave spectras
    time_max = 3600*2  # [s]
    dt = 0.1
    time_frame = np.arange(0, 3600, dt)

    fre_max = 3
    d_fre = 2 * pi / time_max
    fre_range = np.arange(d_fre, fre_max, d_fre)
    # fre_range=np.arange(0.001,3,0.001)
    # fre_range=np.linspace(d_fre,fre_max,10000)
    plt.figure(figsize=(6.3, 4.0))
    Hs = 4
    Tp = 8
    plt.plot(fre_range, jonswap_spectra(fre_range, Hs, Tp, gamma=5), label="jonswap_gama=5")
    plt.plot(fre_range, jonswap_spectra(fre_range, Hs, Tp, gamma=2), label="jonswap_gama=2")
    plt.plot(fre_range, jonswap_spectra(fre_range, Hs, 8.4, gamma=1), label="jonswap_gama=1")
    plt.plot(fre_range, pierson_moskowitz_spectra(fre_range, Hs, Tp), label="PM")
    plt.plot(fre_range, jonswap_spectra(fre_range, Hs, Tp, gamma=3.3), label="jonswap_gama=3.3")
    plt.xlabel("omega (rad)")
    plt.ylabel("S(omega)")
    plt.xlim(0, 3)
    plt.ylim(0, 6)
    plt.grid(True)
    plt.legend()
    plt.savefig('./figures/waveSpectra.png', dpi=600)
    # plt.show()

    # plot 2
    import Airywave as aw
    g1 = 3
    g2 = 1
    gs = gridspec.GridSpec(g1, g2)           # Create 1x2 sub plots

    xi_range = np.sqrt(2*d_fre*jonswap_spectra(fre_range, 4, 8.4, gamma=1.0))
    yita_com = np.zeros((len(xi_range), len(time_frame)))
    yita_com2 = np.zeros((len(xi_range), len(time_frame)))

    space_slice = np.ones((1000, 3))
    x_axis = []
    for posi in range(1000):
        space_slice[posi] = [posi/2, 0, 0]
        x_axis.append(posi/2)
    wave_elevation_with_x = np.zeros((len(xi_range), len(x_axis)))

    list_of_waves = []
    for index, each in enumerate(list(zip(xi_range, fre_range))):
        wave_period = 2*pi/each[1]
        list_of_waves.append(aw.Airywave(2 * each[0], wave_period, 60, 0, random.uniform(0, 360)))

    for index, each in enumerate(list(zip(xi_range, fre_range))):
        yita_com[index, :] = list_of_waves[index].get_elevations_with_time(
            np.array([0, 0, 0]), time_frame)
        yita_com2[index, :] = list_of_waves[index].get_elevations_with_time(
            np.array([100, 0, 0]), time_frame)
        # yita_com[index,:]= each[0]*np.cos(each[1]*time_frame-random.uniform(0,2*pi))
        # yita_com2[index,:]= each[0]*np.cos(each[1]*time_frame-random.uniform(0,2*pi))
    yita = np.sum(yita_com, axis=0)
    yita2 = np.sum(yita_com2, axis=0)
    print("The maximum elevation is "+str(max(yita)))
    print("The minimum elevation is "+str(min(yita)))
    print("The standard deviation is "+str(np.std(yita)))

    print("The maximum elevation is "+str(max(yita2)))
    print("The minimum elevation is "+str(min(yita2)))
    print("The standard deviation is "+str(np.std(yita2)))

    t_ser = np.arange(0, 3600, 1).tolist()
    for t in t_ser:
        for index, each in enumerate(list(zip(xi_range, fre_range))):
            wave_elevation_with_x[index, :] = list_of_waves[index].get_elevation_at_nodes(space_slice, t)
        wave_elevation = np.sum(wave_elevation_with_x, axis=0)

        plt.figure(figsize=(16, 9))
        ax = plt.subplot(gs[0, 0])
        plt.title("surface elvation along x axis when time = "+str(t)+"s")
        plt.plot(x_axis, wave_elevation)
        plt.xlabel("X (m)")
        plt.ylabel("surface elevation (m)")
        plt.xlim(0, 500)
        plt.ylim(-5, 5)

        ax = plt.subplot(gs[1, 0])
        plt.title("surface elvation with time at position x=0,y=0")
        plt.plot(time_frame, yita)
        plt.xlabel("time (s)")
        plt.ylabel("surface elevation (m)")
        plt.xlim(0, 3600)
        plt.ylim(-5, 5)
        plt.tight_layout()
        # plt.show()

        ax = plt.subplot(gs[2, 0])
        plt.title("surface elvation with time at position x=100,y=0")
        plt.plot(time_frame, yita2)
        plt.xlabel("time (s)")
        plt.ylabel("surface elevation (m)")
        plt.xlim(0, 3600)
        plt.ylim(-5, 5)
        plt.tight_layout()
        # plt.show()

        print("time is "+str(t))
        plt.savefig('./figures/waveElevations_'+str(t)+'.png', dpi=300)
