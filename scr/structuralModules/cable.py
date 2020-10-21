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
    def __init__(self, point1_index, point2_index, initial_length, cross_section_diameter, density, Young_module,
                 breaking_strength):
        """ Initialize a cable element
        Args:
            point1 (np.array[1,3]): [description]
            point2 (np.array[1,3]): [description]
            initial_length (float)| uint [m]: the length of the cable that do not have any tension
            cross_section_diameter (float): [description]
            density (float): material density [kg/m3]
            Young_module (float): [description]
            breaking_strength (float):
        """
        self.l_0 = float(initial_length)
        self.p1 = point1_index
        self.p2 = point2_index
        self.length = float(initial_length)
        self.area = 0.25 * pi * pow(cross_section_diameter, 2)
        self.row = float(density)
        self.elasticity = Young_module
        self.bs = breaking_strength
        self.tension_mag = 0
        self.mass=initial_length*self.area*density
        self.element_unit_vector = np.array([0,0,0])

    def cal_property(self, position1, position2):
        """
        A private function, should be invoke by others
        :param position1:
        :param position2:
        :return: No return value, Only change properties of element
        """
        self.element_unit_vector = (position2 - position1) / np.linalg.norm(position2 - position1)
        # epsilon is always > 0, because cable can not be compressed
        epsilon = max((np.linalg.norm(position2 - position1) - self.l_0) / self.l_0, 0)
        self.length=np.linalg.norm(position2 - position1)
        self.tension_mag = self.elasticity * epsilon * self.area
        # break
        if self.tension_mag >self.bs:
            self.elasticity=0
            self.tension_mag=0



    def map_tensions(self, position1, position2):
        self.cal_property(position1, position2)
        force1 = self.element_unit_vector*self.tension_mag
        return force1, -force1


if __name__ == "__main__":
    pass
