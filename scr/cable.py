"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
class for cable element
"""
import numpy as np
from numpy import pi


class cable:
    """using a cable theory\n
    Ref is not defined yet
    """
    def __init__(self, point1, point2, initial_length, cross_section_diameter, density, Young_module, breaking_strength):
        """ Initialize a cable element
        Args:
            point1 (np.array[1,3]): [description]
            point2 (np.array[1,3]): [description]
            initial_length (float): [description]
            cross_section_diameter (float): [description]
            density (float): material density [kg/m3]
            Young_module (float): [description]
            breaking_strength (float):
        """
        self.l_0 = initial_length
        self.p1 = point1
        self.p2 = point2
        self.length = np.linalg.norm(point1 - point2)
        self.area = 0.25 * pi * pow(cross_section_diameter, 2)
        self.row = density
        self.elasticity = Young_module
        self.bs = breaking_strength
        self.tension_mag = self.cal_tension(point1, point2)[0]

    def cal_tension(self, position1, position2):
        element_unit_vector = (position2 - position1) / np.linalg.norm(position2 - position1)
        # epsilon is always > 0, because cable can not be compressed
        epsilon = max((np.linalg.norm(position2 - position1) - self.l_0) / self.l_0, 0)
        tension_mag = self.elasticity * epsilon * self.area
        tension_vector = element_unit_vector * tension_mag
        self.tension_mag = tension_mag
        return tension_mag, tension_vector

    def map_tensions(self,position1, position2):
        tension_mag, tension_vector=self.cal_tension(position1, position2)
        force1 = -tension_vector
        force2 = tension_vector
        return force1, force2

if __name__ == "__main__":
    gravity = 9.81  # [m/s2]
    # define nodes
    nodes=np.zeros((20,3))
    for i in range(20):
        nodes[i]=[i*0.05,0,-i*0.05]
    # define connection
    elements=[]
    for i in range(20-1):
        elements.append([i,i+1])


    structure=[]
    for each in elements:
        structure.append(cable(nodes[each[0]],nodes[each[1]], 0.05, 0.01, 1025, 2e6, 15e6))
    #TODO solve the mass motion equations

    print(structure[0].tension_mag)
    print(structure[0].cal_tension(np.array([0, 0, 0]), np.array([0, 0, 1.1])))
    print(structure[0].map_tensions(np.array([0, 0, 0]), np.array([0, 0, 1.1])))



    import matplotlib.pyplot as plt


    for each in elements:
        plt.plot(nodes[each][:,0], nodes[each][:,2])
    plt.scatter(nodes[:,0],nodes[:,2],color='k')
    plt.show()