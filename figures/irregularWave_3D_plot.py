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
    from scr.wave_spectrum import *
except (ModuleNotFoundError, ImportError) as e:
    print("{} fileure".format(type(e)))
else:
    print("Import succeeded")

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.weight'] = 'regular'
plt.rcParams["mathtext.default"] = "it"
plt.rcParams["mathtext.fontset"] = "stix"

fig = plt.figure()
ax = fig.gca(projection='3d')
# time_frame=np.arange(1000,4600,dt)
# yita=sea_state.get_elevations_with_time([0,0,0],time_frame)
# print(np.std(yita))
# plt.plot(time_frame,yita)
x_axis=np.array(np.arange(-10,10,1.1)).tolist()
position=np.zeros((10*len(x_axis),3))
for i in range(10):
    for index, item in enumerate(x_axis):
        position[i*len(x_axis)+index]=[6*i,6*i,-float(item)/2]
# print(position.shape)
# print(position)
velocity=sea_state.get_velocity_at_nodes(position, 0)
velocity_mag=[]
for each in velocity:
    velocity_mag.append(np.linalg.norm(each))
print("max velocity is " +str(max(velocity_mag)) + " m/s")
print("and the position is "+str(position[velocity_mag.index(max(velocity_mag))]))
# Flatten and normalize
velocity_mag=np.array(velocity_mag)
normal_velo = (velocity_mag.ravel() - velocity_mag.min()) / velocity_mag.ptp()
# print(normal_velo)
# Repeat for each body line and two head lines
c = np.concatenate((normal_velo, np.repeat(normal_velo, 2)))
# Colormap
c = plt.cm.summer(c)
q=ax.quiver(position[:,0],
            position[:,1],
            position[:,2],
            velocity[:,0],
            velocity[:,1],
            velocity[:,2],
            colors=c,
            normalize=False,
            )
# q.set_array(np.linspace(0,2,10))
fig.colorbar(q)

for i in range(60):
    position=np.zeros((len(x_axis),3))
    for index, item in enumerate(x_axis):
        position[index]=[item,i,0]
    yita=sea_state.get_elevation_at_nodes(position, 0)
    ax.plot(x_axis,[i]*len(x_axis),yita,color="b")

# plt.plot(x_axis,yita)
ax.set_title("JONSWAP sea condition Hs=4m, Tp=8s r=3, t="+str(0)+"s")
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(0, 60)
ax.set_ylim(0, 60)
ax.set_zlim(-30, 5)
plt.savefig('./figures/waves.png', dpi=300)
plt.show()