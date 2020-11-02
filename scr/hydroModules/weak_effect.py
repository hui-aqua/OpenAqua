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

def current_reduce_simple(factor,node_position,velocity_fluid):
    """
    modified the current velocity
    :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
    :param velocity_fluid: velocity_fluid: np.array[n,3] | Unit [m/s]| fluid velocity at the coordinates of nodes, n is the number of nodes
    :return: modified_velocity_fluid: np.array[n,3] | Unit [m/s]| fluid velocity at the coordinates of nodes, n is the number of nodes
    """
    origin = np.mean(node_position,axis=0)
    u_mean=np.mean(velocity_fluid,axis=0)
    velocity_modified=np.zeros(velocity_fluid.shape)

    for index in range(len(node_position)):
        if np.dot(u_mean,node_position[index]-origin)>0:
            velocity_modified[index]=velocity_fluid[index]
        else:
            velocity_modified[index] = velocity_fluid[index]*factor
    return velocity_modified

if __name__ == "__main__":
    k=current_reduce_simple(0.8,np.array([[1,0,0],[-1,0,0]]),np.array([[1,0,0],[1,0,0]]))
    print(k)
