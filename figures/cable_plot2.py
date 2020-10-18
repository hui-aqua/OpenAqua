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
    from scr.cable import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

gravity = np.array([0, 0, -9.81])  # [m/s2]
# define nodes
nodes = np.zeros((20, 3))
for i in range(20):
    nodes[i] = [i * 0.05, 0, 0]
# define connection
elements = []
for i in range(20 - 1):
    elements.append([i, i + 1])

structure = []
for each in elements:
    structure.append(cable(nodes[each[0]], nodes[each[1]], 0.05 * 1.4142, 0.01, 1025, 2e9, 15e6))

# TODO solve the mass motion equations
dt = 0.0001  # [s]
t_end = 2.3
position = nodes.copy()
displacement = np.zeros((len(nodes), 3))
acceleration = np.zeros((len(nodes), 3))
velocity = np.zeros((len(nodes), 3))
mass = np.ones((len(nodes), 1)) * structure[0].area * structure[0].l_0 * structure[0].row

force_external = np.zeros((len(nodes), 3))
for i in range(len(nodes) - 1):
    force_external[i] = structure[i].area * structure[i].l_0 * structure[i].row * gravity

t_inst = 0
for i in range(200):
    t_inst += dt
    force_internal = np.zeros((len(nodes), 3))
    for index, item in enumerate(structure):
        position += displacement
        # print(item.tension_mag)
        force_internal[elements[index][0]] = \
            item.map_tensions(position[elements[index][0]], position[elements[index][1]])[0]
        force_internal[elements[index][1]] = \
            item.map_tensions(position[elements[index][0]], position[elements[index][1]])[1]

    acceleration = (force_external + force_internal) / mass
    velocity = acceleration * dt
    displacement = velocity * dt + 0.5 * acceleration * pow(dt, 2)
    print(max(np.linalg.norm(displacement, axis=1)))
    # boundary condition
    if max(np.linalg.norm(displacement, axis=1)) > 0.05*0.1:
        dt/=2
    elif max(np.linalg.norm(displacement, axis=1))<1e-5:
        dt*=2
    else:
        pass
    velocity[0] = 0
    displacement[0] = 0
    velocity[-1] = 0
    displacement[-1] = 0
    print("t_inst is " + str(t_inst))
    print("dt is " + str(dt))


import matplotlib.pyplot as plt

for each in elements:
    plt.plot([0, 1], [0, 0], color='r')
    plt.plot([0, 1], [-1, -1], color='r')
    plt.plot([0, 0], [0, -1], color='r')
    plt.plot([1, 1], [0, -1], color='r')
    plt.plot(position[each][:, 0], position[each][:, 2])
plt.scatter(position[:, 0], position[:, 2], color='k')
plt.show()
