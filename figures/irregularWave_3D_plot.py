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
import sys

try:
    # The insertion index should be 1 because index 0 is this file
    # sys.path.insert(1, '/absolute/path/to/folder/a')  # the type of path is string
    # sys.path.insert(1, 'E:\Hui_Win\Documents\GitHub\OpenAqua')  # the type of path is string
    sys.path.insert(1, '/home/hui/PycharmProjects/OpenAqua')  # the type of path is string
    # because the system path already have the absolute path to folder a
    # so it can recognize file_a.py while searching 
    from scr.irregularwaves import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

hs = 4
tp = 8
r = 3
direction = 45
sea_state = irregular_sea(hs, tp, r, 60, direction)
print(sea_state)

position = np.zeros((20 * 35, 3))
for i in range(20):
    for j in range(35):
        position[i + 20 * j] = [5 * i, 5 * i, 5 - j]

t_inst = 0
velocity = sea_state.get_velocity_at_nodes(position, t_inst)
fig = plt.figure()
ax = fig.gca(projection='3d')

velocity_mag = np.linalg.norm(velocity, axis=1)
normal_velo = (velocity_mag.ravel() - velocity_mag.min()) / velocity_mag.ptp()
# Repeat for each body line and two head lines
c = np.concatenate((normal_velo, np.repeat(normal_velo, 2)))
# Colormap
c = plt.cm.summer(c)
q = ax.quiver(position[:, 0],
              position[:, 1],
              position[:, 2],
              velocity[:, 0],
              velocity[:, 1],
              velocity[:, 2],
              colors=c,
              normalize=False,
              )
colorbar_size = fig.add_axes([0.1, 0.1, 0.03, 0.5])  # left bottom length hight
fig.colorbar(q, colorbar_size, label="Velocity (m/s)")

for i in range(25):
    position2 = np.zeros((50, 3))
    for j in range(50):
        position2[j] = [i * 4, j * 2, 0]
    eta = sea_state.get_elevation_at_nodes(position2, t_inst)
    print()
    ax.plot(position2[:, 0], position2[:, 1], eta, color="b", linewidth=0.5)

ax.set_title("JONSWAP sea condition $Hs$=" + str(hs) + "m, $Tp$=" + str(tp) + "s, $r$=" + str(r) + ", \n"
             "water depth=60m, wave direction=" + str(direction) + " $\degree$, Time=" + str(t_inst) + "s")
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_zlim(-30, 5)
plt.savefig('./png/irregular_waves_3D_' + str(t_inst) + '.png', dpi=300)
plt.show()
