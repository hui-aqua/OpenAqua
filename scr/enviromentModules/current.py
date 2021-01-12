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



def linear_interpolation_current(z):
    if z > np.max(np.array(U['depth'])):
        u_final = U['velocity'][0]
    elif z < np.min(np.array(U['depth'])):
        u_final = U['velocity'][-1]
    else:
        index = np.argmin(np.abs(np.array(U['depth']) - z))
        u1 = U['velocity'][index]
        z1 = U['depth'][index]
        if z > z1:
            z2 = U['depth'][index - 1]
            u2 = U['velocity'][index - 1]
        else:
            z2 = U['depth'][index + 1]
            u2 = U['velocity'][index + 1]

        dz = abs(float(z1 - z2))
        u_final = np.array(u2) * float(abs(z - z1) / dz) + np.array(u1) * abs(z - z2) / dz
    return u_final


def net2netWake(position):
    # get r_net, flow velocity reduction factor for net2net wake effect
    pass


def cage2cageWake():
    # get r_cage, flow velocity reduction factor for cage2cage wake effect
    pass


def get_velocity_at_nodes(list_of_point):
    node_velocity = np.zeros((len(list_of_point), 3))
    for index, position in enumerate(list_of_point.tolist()):
        node_velocity[index] = linear_interpolation_current(position[2])

    print('\n')
    return node_velocity


def current_reduce_simple(factor, node_position, velocity_fluid):
    """
    modified the current velocity
    :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
    :param velocity_fluid: velocity_fluid: np.array[n,3] | Unit [m/s]| fluid velocity at the coordinates of nodes, n is the number of nodes
    :return: modified_velocity_fluid: np.array[n,3] | Unit [m/s]| fluid velocity at the coordinates of nodes, n is the number of nodes
    """
    origin = np.mean(node_position, axis=0)
    u_mean = np.mean(velocity_fluid, axis=0)
    velocity_modified = np.zeros(velocity_fluid.shape)

    for index in range(len(node_position)):
        if np.dot(u_mean, node_position[index] - origin) > 0:
            velocity_modified[index] = velocity_fluid[index]
        else:
            velocity_modified[index] = velocity_fluid[index] * factor
    return velocity_modified


if __name__ == "__main__":
    U = {'depth': [0, -10, -30, -60],
         'velocity': [[0.5, 0, 0],
                      [0.6, 0, 0],
                      [0.3, 0, 0],
                      [0.0, 0, 0]]
         }

    points = np.array([[100, 0, 1],
                       [150, 0, -5],
                       [150, 0, -10],
                       [150, 0, -15],
                       [100, 0, -31],
                       [150, 0, -59],
                       [111, 5, -80]])
    u = get_velocity_at_nodes(points)
    print(u)
    print(U.keys())
    # k=current_reduce_simple(0.8,np.array([[1,0,0],[-1,0,0]]),np.array([[1,0,0],[1,0,0]]))
    # print(k)
