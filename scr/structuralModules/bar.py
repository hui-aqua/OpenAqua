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
from scr.structuralModules.cable import *

class bar(cable):
    def cal_tension(self, position1, position2):
        element_unit_vector = (position2 - position1) / np.linalg.norm(position2 - position1)
        # epsilon is always > 0, because cable can not be compressed
        epsilon = (np.linalg.norm(position2 - position1) - self.l_0) / self.l_0
        tension_mag = self.elasticity * epsilon * self.area
        tension_vector = element_unit_vector * tension_mag
        self.tension_mag = tension_mag
        return tension_mag, tension_vector