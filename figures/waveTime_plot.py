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

# validation 2 shows the wave velocity & acceleration according to time

g1 = 2
g2 = 3
gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots
wave1=Airywave(5,10,600,0)

time_max = 60  # [s]
dt = 0.2
time_frame = np.arange(0, time_max, dt)

positions = np.array([[0,0,2],[0,0,0],[0,0,-2]])


velocities=np.zeros((len(positions),len(time_frame),3))
acceleration=np.zeros((len(positions),len(time_frame),3))
eta=wave1.get_elevations_with_time(positions[0],time_frame)
for i in range(len(positions)):
    velocities[i] = wave1.get_velocity_with_time(positions[i],time_frame)
    acceleration[i] = wave1.get_acceleration_with_time(positions[i],time_frame)

plt.figure()
for i in range(3):
    ax = plt.subplot(gs[0, i])
    plt.title("Wave elevation with time at x=0m, y=0m and velocity u"+str(i))
    for j in range(len(positions)):
        plt.plot(time_frame, velocities[j,:,i], label="position_z="+str(positions[j,2]))
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")
    plt.xlim(0, 60)
    plt.ylim(-3, 3)
    plt.grid(True)
plt.legend(frameon=False, loc='center', bbox_to_anchor=(1.2, 0.5))

for i in range(3):
    ax = plt.subplot(gs[1, i])
    plt.title("Wave elevation with time at x=0m, y=0m and acceleration u"+str(i))
    for j in range(len(positions)):
        plt.plot(time_frame, acceleration[j,:,i], label="position_z="+str(positions[j,2]))
    plt.xlabel("Time (s)")
    plt.ylabel("acceleration (m/s$^{2}$)")
    plt.xlim(0, 60)
    plt.ylim(-3, 3)
    plt.grid(True)
    # plt.legend(frameon=False, loc='center', bbox_to_anchor=(1.2, 0.5))

plt.tight_layout()
plt.savefig('./png/wave_time.png', dpi=600)
plt.show()