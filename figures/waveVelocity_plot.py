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
    from scr.Airywave import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

g1 = 2
g2 = 1
gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots

wave1 = Airywave(6, 8, 600, 0, 0)

print(wave1)
x_list = np.linspace(0, 90, 10)
z_list = np.linspace(5, -60, 30)
x_axis = np.linspace(0, 500, 500)

yita_list = []
for x in x_axis:
    yita_list.append(wave1.get_elevation(np.array([x, 0, 0]), 0))

posi = []
for x in x_list:
    for z in z_list:
        posi.append([x, 0, z])
posi = np.array(posi)
velo = wave1.get_velocity_at_nodes(posi, 0)
acce = wave1.get_acceleration_at_nodes(posi, 0)

print("maximum velocity mag is " + str(np.max(np.linalg.norm(velo, axis=1))))

plt.figure(figsize=(6.3, 4.5))

ax = plt.subplot(gs[0, 0])
plt.title("Velocity")
plt.plot(x_axis, yita_list, color="b")
plt.quiver(posi[:, 0], posi[:, 2], velo[:, 0], velo[:, 2],units='xy')
plt.text(65, -80, "Maximum velocity " + str(round(max(np.linalg.norm(velo, axis=1)), 2)) + " m/s")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.xlim(-10, 100)
plt.ylim(-60, 10)
plt.grid(True)

ax = plt.subplot(gs[1, 0])
plt.title("Acceleration")
plt.plot(x_axis, yita_list, color="b")
plt.quiver(posi[:, 0], posi[:, 2], acce[:, 0], acce[:, 2],units='xy')
plt.text(65, -80, "Maximum acceleration " + str(round(max(np.linalg.norm(acce, axis=1)), 2)) + " m/s$^{2}$")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.xlim(-10, 100)
plt.ylim(-60, 10)
plt.grid(True)

plt.tight_layout()
plt.savefig('./png/wave_velocity.png', dpi=600)
plt.show()
