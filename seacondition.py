"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
A module can be used to provide the wave velocity
In order to use this module, we recommend ``import seacondition as sc`` in the beginning of your code.

    .. note:: In this library, I assume the wave is coming from -X, and long X-axis, Z is the gravity direction,Z=0 is the water level, and the water is below z=0.

reference: 2000linearwavetheory_NTNU.pdf
"""
import numpy as np
# from scr.model4aster.waves import Airywave
import waves.Airywave as Airywave
# from scr.model4aster.waves import stokeswave
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    yita=2
    w_d=60
    z=np.arange(-w_d,0,1)
    zs=np.arange(-w_d,yita,1)
    z2=zs-yita/(1+yita/w_d)
    plt.plot(z,z,label='z')
    plt.plot(zs,z2,label='z2')
    plt.legend()
    
    plt.show()
