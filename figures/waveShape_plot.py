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

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

## validation 2 shows the wave elevation according to time and space
g1 = 2
g2 = 1
gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots

water_d = [10, 20, 30, 40, 50, 60, 80, 100, 1000]

time_slice = np.linspace(0, 100, 1000)
velocities_with_time = []

space_slice = np.ones((1000, 3))
x_axis = []
for position in range(1000):
    space_slice[position] = [position / 2, 0, 0]
    x_axis.append(position / 2)
elevation_with_time = []
wave_height = 1.5
wave_period = 10

for item in water_d:
    velocities_with_time.append(
        [Airywave(wave_height, wave_period, item, 0).get_elevation(np.array([0, 0, 0]), i) for i in time_slice])
    elevation_with_time.append(Airywave(wave_height, wave_period, item, 0).get_elevation_at_nodes(space_slice, 0))
plt.figure()

ax = plt.subplot(gs[0, 0])
plt.title("Wave elevation with time at x=0m, y=0m")
for item in water_d:
    plt.plot(time_slice, velocities_with_time[water_d.index(item)], label="Depth " + str(item))
plt.xlabel("Time (s)")
plt.ylabel("Wave elevation (m)")
plt.xlim(0, 100)
plt.ylim(-3, 3)
plt.grid(True)
plt.legend(frameon=False, loc='center', bbox_to_anchor=(1.2, 0.5))

ax = plt.subplot(gs[1, 0])
plt.title("Wave elevation with space when T=0")
for item in water_d:
    plt.plot(x_axis, elevation_with_time[water_d.index(item)], label="Depth " + str(item))
plt.xlabel("X (m)")
plt.ylabel("Wave elevation (m)")
plt.xlim(0, 500)
plt.ylim(-3, 3)
plt.grid(True)

plt.tight_layout()
plt.savefig('./png/wave_shape.png', dpi=600)
plt.show()