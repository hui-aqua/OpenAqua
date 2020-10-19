"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Plot figure(s)
please email: hui.cheng@uis.no \n
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

import matplotlib
# del matplotlib.font_manager.weight_dict['roman']
# matplotlib.font_manager._rebuild()

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
    from scr.enviromentModules.Airywave import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

g1 = 1
g2 = 2
gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots

water_d = [10, 20, 30, 40, 50, 60, 80, 100, 1000]
waves_length = []
waves_phasevelocity = []
wave_period = np.linspace(1, 20, 50)
for item in water_d:
    waves_length.append([Airywave(1, i, item, 0).wave_Length for i in wave_period])
    waves_phasevelocity.append([Airywave(1, i, item, 0).wave_phase_velocity for i in wave_period])

plt.figure(figsize=(6.3, 3.5))

ax = plt.subplot(gs[0, 0])
for item in waves_length:
    plt.plot(wave_period, item, label="Depth " + str(water_d[waves_length.index(item)]))
plt.xlabel("Wave period (s)")
plt.ylabel("Wave length (m)")
plt.xlim(0, 20)
plt.ylim(0, 700)
plt.grid(True)

ax = plt.subplot(gs[0, 1])
for item in waves_phasevelocity:
    plt.plot(wave_period, item, label="Depth " + str(water_d[waves_phasevelocity.index(item)]))
plt.xlabel("Wave period (s)")
plt.ylabel("Phase velocity (m/s)")
plt.xlim(0, 20)
plt.ylim(0, 35)
plt.grid(True)
plt.legend(frameon=False, loc='center', bbox_to_anchor=(1.4, 0.5))

plt.tight_layout()

plt.savefig('./png/wave_period.png', dpi=600)
plt.show()
## ref to Figure3-3 on page 47 from DNV GL-RP205 Ver. 2008