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
    from scr.structuralModules.cable import *
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
nodes = np.zeros((15, 3))
for i in range(15):
    nodes[i] = [i * 0.1, 0, 0]
# define connection
elements = []
for i in range(15 - 1):
    elements.append([i, i + 1])

structure = []
for each in elements:
    structure.append(cable(each[0], each[1], 0.1, 0.01, 1025, 2e6, 15e6))


dt = 1e-4   # [s]
t_end = 0.1
position = nodes.copy()
displacement = np.zeros((len(nodes), 3))
velocity = np.zeros((len(nodes), 3))
total_force = np.zeros((len(nodes), 3))
mass = np.ones((len(nodes), 1)) * structure[0].area * structure[0].l_0 * structure[0].row

force_external = np.zeros((len(nodes), 3))
for i in range(len(nodes) - 1):
    force_external[i] = structure[i].area * structure[i].l_0 * structure[i].row * gravity

# boundary condition on the end
force_external[-1] = [0, 0, -1]

t_inst = 0
for i in range(int(t_end / dt)):
    t_inst += dt
    force_internal = np.zeros((len(nodes), 3))
    displacement = velocity * dt
    for index, item in enumerate(structure):
        position += displacement
        print(item.tension_mag)
        force_internal[elements[index][0]] = \
            item.map_tensions(position[elements[index][0]], position[elements[index][1]])[0]
        force_internal[elements[index][1]] = \
            item.map_tensions(position[elements[index][0]], position[elements[index][1]])[1]
    total_force=force_external + force_internal
    velocity = (total_force) / mass * dt
    # velocity=np.where(abs(velocity)<0.1/dt,velocity,velocity*0.01)

    # boundary condition
    velocity[0] = 0
    displacement[0] = 0

    print("t_inst is " + str(t_inst))
    # print(position)
    # plt.clf()

for each in elements:
    plt.plot([0, 1], [0, 0], color='r')
    plt.plot([0, 1], [-1, -1], color='r')
    plt.plot([0, 0], [0, -1], color='r')
    plt.plot([1, 1], [0, -1], color='r')
    plt.plot(position[each][:, 0], position[each][:, 2])
plt.axis('equal')
plt.scatter(position[:, 0], position[:, 2], color='k')
plt.show()
