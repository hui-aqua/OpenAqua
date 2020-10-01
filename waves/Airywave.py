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
    """
    def __init__(self, waveHeight=1.0, waveLength=25.0, waterDepth=60.0):
        """
        :param waveHeight: [float] Unit: [m]. wave height.
        :param waveLength: [float] Unit: [m]. wavelength.
        :param waterDepth: [float] Unit: [m]. water depth.
        """
        self.gravity = 9.81
        self.waveHeight = waveHeight
        self.waveLength = waveLength
        self.waterDepth = waterDepth
        self.phase_velocity = np.sqrt(
            self.gravity * waveLength * np.tanh(2 * pi * waterDepth / waveLength) / (2.0 * pi))
        self.wavePeriod = waveLength / self.phase_velocity
        self.wave_k = 2 * pi / waveLength
        self.pi_h_t = pi * waveHeight / self.wavePeriod
        self.pi_h_l = pi * waveHeight / waveLength
        self.pi_h_t_2 = 2 * waveHeight * pow(pi / self.wavePeriod, 2)

    def __str__(self):
        s0 = 'The environment is airy wave (deep water) wave condition and the specific parameters are:\n'
        s1 = 'water Depth= ' + str(self.waterDepth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wavePeriod) + ' s\n'
        s3 = 'wave Length= ' + str(self.waveLength) + ' m\n'
        s4 = 'wave Height= ' + str(self.waveHeight) + ' m\n'
        S = s0 + s1 + s2 + s3 + s4
        return S

    def get_velocity(self, posi, time):
        """
        :param posi: [np.array].shape=(1,3) or a [list] of coordinates Unit: [m]. The position of the point which you want to know the wave velocity
        :param time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return:  [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the velocity at the targeted point.
        """
        zeta = self.wave_k * posi[0] - 2 * pi / self.wavePeriod * time
        yita = self.get_surface(posi, time)
        if posi[2] < yita:
            horizonvelocity = self.pi_h_t * np.cosh(self.wave_k * (posi[2] + self.waterDepth)) * np.cos(zeta) / np.sinh(
                self.wave_k * self.waterDepth)
            vericalvelocity = self.pi_h_t * np.sinh(self.wave_k * (posi[2] + self.waterDepth)) * np.sin(zeta) / np.sinh(
                self.wave_k * self.waterDepth)

            if self.waterDepth > self.waveLength * 0.5:  # if it is deep water
                horizonvelocity = self.pi_h_t * np.exp(self.wave_k * posi[2]) * np.cos(zeta)
                vericalvelocity = self.pi_h_t * np.exp(self.wave_k * posi[2]) * np.sin(zeta)
        else:
            horizonvelocity = 0.0
            vericalvelocity = 0.0
        velo = np.array([0.0, 0.0, 0.0])
        velo[0] = horizonvelocity
        velo[1] = 0.0
        velo[2] = vericalvelocity
        return velo

    def get_acceleration(self, posi, time):
        """
        :param posi: [np.array].shape=(1,3) or a [list] of coordinates Unit: [m]. The position of the point which you want to know the wave acceleration.
        :param time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return: [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the acceleration at the targeted point.
        """

        zeta = self.wave_k * posi[0] - 2 * pi / self.wavePeriod * time
        yita = self.get_surface(posi, time)
        if posi[2] < yita:
            horizontalacceleration = self.pi_h_t_2 * np.cosh(self.wave_k * (posi[2] + self.waterDepth)) * np.sin(
                zeta) / np.sinh(self.wave_k * self.waterDepth)
            vericalaccelateration = -self.pi_h_t_2 * np.sinh(self.wave_k * (posi[2] + self.waterDepth)) * np.cos(
                zeta) / np.sinh(self.wave_k * self.waterDepth)
            if self.waterDepth > self.waveLength * 0.5:  # if it is deep water
                horizontalacceleration = self.pi_h_t_2 * np.exp(self.wave_k * posi[2]) * np.sin(zeta)
                vericalaccelateration = -self.pi_h_t_2 * np.exp(self.wave_k * posi[2]) * np.cos(zeta)
        else:
            horizontalacceleration = 0.0
            vericalaccelateration = 0.0
        acce = np.array([0.0, 0.0, 0.0])
        acce[0] = horizontalacceleration
        acce[1] = 0.0
        acce[2] = vericalaccelateration
        return acce

    def get_surface(self, positions, time):
        """
        A class private function. \n
        :param positions: [np.array].shape=(1,3) or a [list] of coordinates Unit: [m]. The position of the point which you want to know the wave velocity or acceleration.
        :param time: [float] Unit: [s].
        :return: [float] Unit: [m]. The sea surface level in Z direction. At the targeted position.
        """
        zeta = self.wave_k * positions[0] - 2 * pi / self.wavePeriod * time
        yita = self.waveHeight / 2 * np.cos(zeta)
        return yita

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
