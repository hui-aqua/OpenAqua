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
import Airywave as aw

class Stocks2wave(aw.Airywave):
    def __str__(self):
        s0 = 'The environment is stokes second order wave theory condition and the specific parameters are:\n'
        s1 = 'water Depth= ' + str(self.waterDepth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wavePeriod) + ' s\n'
        s3 = 'wave Length= ' + str(self.waveLength) + ' m\n'
        s4 = 'wave Height= ' + str(self.waveHeight) + ' m\n'
        S = s0 + s1 + s2 + s3 + s4
        return S
#TODO modified the original function by adding additional nonlinear terms.
if __name__ == "__main__":
    pass