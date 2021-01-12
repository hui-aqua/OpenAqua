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
from scr.enviromentModules import ocean
import numpy as np
# importing the module
import ast

# reading the data from the file
with open('dataddd.py') as f:
    data = f.read()

print("Data type before reconstruction : ", type(data))

# reconstructing the data as a dictionary
meshinfo = ast.literal_eval(data)

print("Data type after reconstruction : ", type(meshinfo))
print(meshinfo['multi_ring'])
print(meshinfo['multi_net'])
U = {'depth': [0, -10, -30, -60],
     'velocity': [[0.5, 0, 0],
                  [0.6, 0, 0],
                  [0.3, 0, 0],
                  [0.0, 0, 0]]
     }


rings=meshinfo['multi_ring']
net=meshinfo['multi_net']
point=np.array(meshinfo['Nodes'])
current=ocean.current(U,rings,net,len(point),0.2)

print('number of nodes is')
print(len(point))
print('number of nodes on one cage nets  '+str(len(net['cage_0'])))
u=current.get_velocity_at_nodes(point)
print(u)