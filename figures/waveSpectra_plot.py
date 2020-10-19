import numpy as np
import matplotlib.pyplot as plt
import sys

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

try:
    # The insertion index should be 1 because index 0 is this file
    # sys.path.insert(1, '/absolute/path/to/folder/a')  # the type of path is string
    # sys.path.insert(1, 'E:\Hui_Win\Documents\GitHub\OpenAqua')  # the type of path is string
    sys.path.insert(1, '/home/hui/PycharmProjects/OpenAqua')  # the type of path is string
    # because the system path already have the absolute path to folder a
    # so it can recognize file_a.py while searching 
    from scr.enviromentModules.wave_spectrum import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

# # Plot1  wave spectras
time_max = 3600 * 2  # [s]
dt = 0.1
time_frame = np.arange(0, 3600, dt)
fre_max = 3
d_fre = 2 * np.pi / time_max
fre_range = np.arange(d_fre, fre_max, d_fre)
plt.figure(figsize=(6.3, 4.0))
Hs = 4
Tp = 8
plt.plot(fre_range, jonswap_spectra(fre_range, Hs, Tp, gamma=5),   label="JONSWAP_gama=5")
plt.plot(fre_range, jonswap_spectra(fre_range, Hs, Tp, gamma=2),   label="JONSWAP_gama=2")
plt.plot(fre_range, jonswap_spectra(fre_range, Hs, 8.4, gamma=1),  label="JONSWAP_gama=1")
plt.plot(fre_range, jonswap_spectra(fre_range, Hs, Tp, gamma=3.3), label="JONSWAP_gama=3.3")
plt.plot(fre_range, pierson_moskowitz_spectra(fre_range, Hs, Tp), label="PM")
plt.xlabel("Omega (rad)")
plt.ylabel("$S$")
plt.xlim(0, 3)
plt.ylim(0, 6)
plt.grid(True)
plt.legend()
# print(sys.path[0])
plt.savefig('./png/waveSpectra.png', dpi=600)
plt.show()
