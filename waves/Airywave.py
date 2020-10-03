"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
wave direction is x+
"""
import numpy as np
from numpy import pi


class Airywave:
    """
    Using Airy wave theory      \n
    Ref. DNV GL-RP205 Ver. 2008:P45
    Linear wave theory:
    
    """
    def __init__(self, waveHeight=1.0, wavePeriod=10, waterDepth=60.0, direction=0):
        """
        :param waveHeight: [float] Unit: [m]. wave height.
        :param wavePeriod: [float] Unit: [s]. wave period.
        :param waterDepth: [float] Unit: [m]. wave depth.
        :param direction： [float] Unit: [degree]. direction of propagation, measured from the positive x-axis.
        """
        self.gravity = 9.81
        self.wave_Height = waveHeight
        self.wave_Period = wavePeriod
        self.water_Depth = waterDepth
        self.beta = pi*direction/180.0

        ## calculation
        # angular frequency 
        self.omega=2*pi/wavePeriod
        # wave length
        alpha=[1,0.666,0.445,-0.105,0.272]
        omega_ba=4.0*pow(pi,2)*waterDepth/self.gravity/pow(wavePeriod,2)
        f_omega=0.0
        for i,a in enumerate(alpha):
            f_omega+=a*pow(omega_ba,i)
        self.wave_Length=self.wave_Period*pow(self.gravity*self.water_Depth,0.5)*pow(f_omega/(1+omega_ba*f_omega),0.5)
        # wave number 
        self.wave_k=2*pi/self.wave_Length
        # phase velocity
        self.c=pow(self.gravity/self.wave_k*np.tanh(self.wave_k*self.water_Depth),0.5)
        
        ### for easy calculation       
        self.pi_h_t = pi * waveHeight / wavePeriod
        self.pi_h_l = pi * waveHeight / self.wave_Length
        self.pi_h_t_2 = 2 * waveHeight * pow(pi / wavePeriod, 2)

    def __str__(self):
        s0 = 'The environment is airy wave (deep water) wave condition and the specific parameters are:\n'
        s1 = 'water Depth= ' + str(self.water_Depth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wave_Period) + ' s\n'
        s3 = 'wave Length= ' + str(self.wave_Length) + ' m\n'
        s4 = 'wave Height= ' + str(self.wave_Height) + ' m\n'
        s4 = 'wave phase velocity= ' + str(self.c) + ' m/s\n'
        S = s0 + s1 + s2 + s3 + s4
        return S

    def get_surface(self, positions, time):
        """
        A private function. \n
        :param positions: [np.array].shape=(1,3) coordinates Unit: [m]. The position of the point which you want to know the wave surface elevation.
        :param time: [float] Unit: [s].
        :return: [float] Unit: [m]. The sea surface level in Z direction. At the targeted position.
        """
        zeta = self.wave_k * (positions[0]*np.cos(self.beta)+positions[1]*np.sin(self.beta)) -self.omega* time
        yita = float(self.wave_Height / 2 * np.cos(zeta))
        return yita

    def get_velocity(self, positions, time):
        """
        :param positions: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave velocity
        :param time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return:  [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the velocity at the targeted point.
        """
        zeta = self.wave_k * (positions[0]*np.cos(self.beta)+positions[1]*np.sin(self.beta)) -self.omega* time
        yita = self.get_surface(positions, time)
        if positions[2] < yita:
            horizonvelocity = self.pi_h_t * np.cosh(self.wave_k * (positions[2] + self.water_Depth)) * np.cos(zeta) / np.sinh(self.wave_k * self.water_Depth)
            vericalvelocity = self.pi_h_t * np.sinh(self.wave_k * (positions[2] + self.water_Depth)) * np.sin(zeta) / np.sinh(self.wave_k * self.water_Depth)
        else:
            horizonvelocity = 0.0
            vericalvelocity = 0.0
            
        velo = np.array([0.0, 0.0, 0.0])
        velo[0] = horizonvelocity*np.cos(self.beta)
        velo[1] = horizonvelocity*np.sin(self.beta)
        velo[2] = vericalvelocity
        return velo

    def get_acceleration(self, positions, time):
        """
        :param posi: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave acceleration.
        :param time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return: [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the acceleration at the targeted point.
        """

        zeta = self.wave_k * (positions[0]*np.cos(self.beta)+positions[1]*np.sin(self.beta)) -self.omega* time
        yita = self.get_surface(positions, time)
        if positions[2] < yita:
            horizontalacceleration = self.pi_h_t_2 * np.cosh(self.wave_k * (positions[2] + self.water_Depth)) * np.sin(zeta) / np.sinh(self.wave_k * self.water_Depth)
            vericalaccelateration = -self.pi_h_t_2 * np.sinh(self.wave_k * (positions[2] + self.water_Depth)) * np.cos(zeta) / np.sinh(self.wave_k * self.water_Depth)
        else:
            horizontalacceleration = 0.0
            vericalaccelateration = 0.0
        acce = np.array([0.0, 0.0, 0.0])
        acce[0] = horizontalacceleration*np.cos(self.beta)
        acce[1] = horizontalacceleration*np.sin(self.beta)
        acce[2] = vericalaccelateration
        return acce


    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        list_of_velocity = []
        for point in list_of_point:
            velocity_at_point = self.get_velocity(point, global_time)
            list_of_velocity.append(velocity_at_point)
        return list_of_velocity

    def get_velocity_at_elements(self, position_nodes, elements, global_time):
        """
        :param position_nodes: a numpy list of position \n
        :param elements: a python list of element \n
        :param global_time: time [s] \n
        :return: Get a numpy array of velocity at a list of elements \n
        """
        velocity_list = []
        for element in elements:
            element_center = np.array([0, 0, 0])
            for node in element:
                element_center += position_nodes[node] / len(element)
            velocity_on_element = self.get_velocity(element_center, global_time)
            velocity_list.append(velocity_on_element)
        return np.array(velocity_list)
    
    
    def get_elevation_at_elements(self, position_nodes,global_time):
        """
        :param position_nodes: a numpy list of position \n
        :param elements: a python list of element \n
        :param global_time: time [s] \n
        :return: Get a numpy array of velocity at a list of elements \n
        """
        velocity_list = []
        for element in elements:
            element_center = np.array([0, 0, 0])
            for node in element:
                element_center += position_nodes[node] / len(element)
            velocity_on_element = self.get_velocity(element_center, global_time)
            velocity_list.append(velocity_on_element)
        return np.array(velocity_list)
    
    
    
if __name__ == "__main__":
    wave1=Airywave(2,20,10,0)
    