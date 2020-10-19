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
# NOTE remember to remove the last line in airywave get velocity/acceleration at nodes
try:
    # The insertion index should be 1 because index 0 is this file
    # sys.path.insert(1, '/absolute/path/to/folder/a')  # the type of path is string
    # sys.path.insert(1, 'E:\Hui_Win\Documents\GitHub\OpenAqua')  # the type of path is string
    sys.path.insert(1, '/home/hui/PycharmProjects/OpenAqua')  # the type of path is string
    # because the system path already have the absolute path to folder a
    # so it can recognize file_a.py while searching 
    from scr.enviromentModules.irregularwaves import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

plt.rcParams['image.cmap'] = 'summer'

sea_state = irregular_sea(5, 8, 3, 60, 0)
time_max = 3600  # [s]
dt = 1
time_frame = np.arange(0, time_max, dt)
print(sea_state)

g1 = 3
g2 = 1
gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots
elevations = sea_state.get_elevations_with_time(np.array([0, 0, 0]), time_frame)
position = np.zeros((600, 3))
for i in range(20):
    for j in range(30):
        position[i + 20 * j] = [5 * i, 0, 5.0 - j]
position.flags.writeable = False
position2 = np.zeros((100, 3))
for i in range(100):
    position2[i] = [i, 0, 0]
plt.figure(figsize=(6.7, 6.5))
# for t in time_frame:
t=0
plt.clf()

ax = plt.subplot(gs[0, 0])
plt.title("Elevation with time at position x=0, y=0")
plt.plot(time_frame, elevations, color='b', linewidth=0.5)
plt.xlabel("time (s)")
plt.ylabel("surface elevation (m)")
plt.xlim(0, 3600)
plt.ylim(-5, 5)

ax = plt.subplot(gs[1, 0])
plt.title("Elevation and velocity at " + str(t) + " s")
velocity = sea_state.get_velocity_at_nodes(position, t)

plt.text(65, -35, "Maximum velocity " + str(round(max(np.linalg.norm(velocity, axis=1)), 2)) + " m/s")
plt.quiver(position[:, 0],
           position[:, 2],
           velocity[:, 0],
           velocity[:, 2],
           units='xy',
           scale=0.5,
           )
plt.plot([i for i in range(100)], sea_state.get_elevation_at_nodes(position2, t), color='b')
plt.plot([0, 100], [0, 0], color='r', linewidth=0.5)
plt.xlabel("X (m)")
plt.ylabel("elevation (m)")
plt.xlim(0, 100)
plt.ylim(-26, 6)


ax = plt.subplot(gs[2, 0])
plt.title("Elevation and acceleration at " + str(t) + " s")
acceleration=sea_state.get_acceleration_at_nodes(position,t)
plt.text(65, -35, "Maximum acceleration " + str(round(max(np.linalg.norm(acceleration, axis=1)), 2)) + " m/s$^{2}$")
plt.quiver(position[:, 0],
           position[:, 2],
           acceleration[:, 0],
           acceleration[:, 2],
           units='xy',
           scale=0.5,
           )
plt.plot([i for i in range(100)], sea_state.get_elevation_at_nodes(position2, t), color='b')
plt.plot([0, 100], [0, 0], color='r', linewidth=0.5)
plt.xlabel("X (m)")
plt.ylabel("elevation (m)")
plt.xlim(0, 100)
plt.ylim(-26, 6)
    
# print("velocity at "+str(position[360])+" is "+ str(velocity[360]))
plt.tight_layout()
plt.savefig('./png/irregular_waves_' + str(t) + '.png', dpi=600)
# plt.savefig('./png/timeframes/irrugularwaves_' + str(t) + '.png', dpi=600)
plt.show()
