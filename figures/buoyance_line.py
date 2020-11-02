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
import matplotlib.gridspec as gridspec
import sys

try:
    # The insertion index should be 1 because index 0 is this file
    # sys.path.insert(1, '/absolute/path/to/folder/a')  # the type of path is string
    # sys.path.insert(1, 'E:\Hui_Win\Documents\GitHub\OpenAqua')  # the type of path is string
    sys.path.insert(1, '/home/hui/PycharmProjects/OpenAqua')  # the type of path is string
    # because the system path already have the absolute path to folder a
    # so it can recognize file_a.py while searching 
    from scr.hydroModules.one_dimensional import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"


# translation
elevation = np.array([[0,0,0],[1,0,0]])
line1=MorisonModel('M4',[[1,2]],0.2,0.01,0.01)
num_p=1000
z=np.linspace(-2,2,num_p)
# z = np.linspace(-3, 3, num_p)
posi=np.zeros((num_p,2,3))
buoy=np.zeros((num_p,3))
for i in range(num_p):
    posi[i]=np.array([[0,0,z[i]],[1,0,z[i]-0.1]])
    buoy[i] = line1.cal_buoy_force(posi[i],elevation)
# print(buoy)
print(pi*0.25*0.01*0.01)
plt.figure()
plt.plot(z,buoy)
plt.show()
#
# rotate

num_p=1000
x = np.cos(np.linspace(0,2*pi,num_p))*0.5
z = np.sin(np.linspace(0, 2 * pi, num_p)) * 0.5
posi=np.zeros((num_p,2,3))
buoy=np.zeros((num_p,3))
for i in range(num_p):
    elevation = np.array([[x[i], 0, 0], [-x[i], 0, 0]])
    posi[i]=np.array([[x[i],0,z[i]+0.002],[-x[i],0,-z[i]+0.002]])
    buoy[i] = line1.cal_buoy_force(posi[i],elevation)
# print(buoy)
plt.figure()
plt.plot(np.linspace(0,360,num_p),buoy)
plt.show()

