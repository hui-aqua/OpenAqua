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
    gravity = np.array([0,0,-9.81])  # [m/s2]
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
        structure.append(cable(nodes[each[0]],nodes[each[1]], 0.05*1.4142, 0.01, 1025, 2e6, 15e6))
    
    #TODO solve the mass motion equations
    dt=0.1 #[s]
    t_end=1
    
    displacement=np.zeros((len(nodes),3))
    acceleration=np.zeros((len(nodes),3))
    velocity=np.zeros((len(nodes),3))
    mass=np.ones((len(nodes),1))*structure[0].area*structure[0].l_0*structure[0].row
    
    force_external=np.zeros((len(nodes),3))
    for i in range(len(nodes)-1):
        force_external[i]=structure[i].area*structure[i].l_0*structure[i].row*gravity
    # boundary condition on the end        
    force_external[-1]=[0,0,-1]
    
    t_inst=0
    for i in range(int(t_end/dt)):
        t_inst+=dt
        force_internal=np.zeros((len(nodes),3))
        for index, item in enumerate(structure):
            position=nodes.copy()+displacement
            force_internal[elements[index][0]]= item.map_tensions(position[elements[index][0]],position[elements[index][1]])[0]
            force_internal[elements[index][1]]= item.map_tensions(position[elements[index][0]],position[elements[index][1]])[1]

        acceleration=(force_external+force_internal)/mass
        velocity+=acceleration*dt
        displacement=velocity*dt+0.2*acceleration*pow(dt,2)
        # boundary condition
        velocity[0]=0
        displacement[0]=0
        position+=displacement
        print("t_inst is " +str(t_inst))
        print(position)




    import matplotlib.pyplot as plt


    for each in elements:
        plt.plot(position[each][:,0], position[each][:,2])
    plt.scatter(position[:,0],position[:,2],color='k')
    plt.show()