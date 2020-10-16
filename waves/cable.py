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
    """usig a cabel theory\n
    Ref is not defined yet
    """
    def __init__(self, point1,point2,initialLenght, corss_section_diameter, Young_module, breaking_strength):
        """ Initialize a cable element
        Args:
            point1 (np.array[1,3]): [description]
            point2 (np.array[1,3]): [description]
            initialLenght (float): [description]
            Young_module (float): [description]
            corss_section_diameter (float): [description]
        """
        self.p1 = point1
        self.p2 = point2
        self.length=np.linalg.norm(point1-point2)
        self.area=0.25*pi*pow(corss_section_diameter,2)
        self.elasticity= Young_module
        self.bs=breaking_strength
        