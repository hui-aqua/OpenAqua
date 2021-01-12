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
from scr.enviromentModules import Airywave as wave
from scr.enviromentModules import wave_spectrum as wsp

import numpy as np
from numpy import pi


class irregular_sea:
    """
    The default wave spectra is JONSWAP
    ----------
    *NOTE: The maximum applied simulation time should be less than 3h.
    """

    def __init__(self, significant_wave_height, peak_period, gamma, water_depth, wave_direction):
        """
        Parameters
        ------------
        significant_wave_height: significant wave height | float | Unit [m]
        peak_period: peak period | float | Unit [s]
        gamma: gamma | float | Unit [-]
        water_depth: water depth of the sea, assume flat sea floor. A position number | float | Unit [m]
        wave_direction: direction of wave propagation. | float | Unit [degree]
        """
        self.hs = significant_wave_height
        self.tp = peak_period
        self.water_depth = water_depth
        self.list_of_waves = []

        time_max = 3600 * 3  # 3h, We assume the simulations will not exceed 3h.
        fre_max = 3  # we assume the highest eigen frequency of studied structure is below 3 Hz.
        d_fre = 2 * np.pi / time_max  # get the resolution for frequency.
        fre_range = np.arange(d_fre, fre_max, d_fre)
        design_wave_spectra = wsp.jonswap_spectra(fre_range, significant_wave_height, peak_period, gamma)
        list_xi = np.sqrt(2 * d_fre * design_wave_spectra)
        for index, item in enumerate(list_xi):
            wave_period = 2 * np.pi / fre_range[index]
            self.list_of_waves.append(
                wave.Airywave(item * 2, wave_period, water_depth, wave_direction, np.random.uniform(0, 360)))

    def __str__(self):
        """ Print the information of the present object. """
        s0 = 'The environment is irregular waves condition and the specific parameters are:\n'
        s1 = 'significant wave height = ' + str(self.hs) + ' m\n'
        s2 = 'peak period= ' + str(self.tp) + ' s\n'
        s3 = 'Number of wave components is ' + str(len(self.list_of_waves))
        return s0 + s1 + s2 + s3

    ## elevation

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) | Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations \n
        :return: Get a list of elevations at one position with a time sequence \n
        """
        wave_elevations = np.zeros((len(self.list_of_waves), len(time_list)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index] = each_wave.get_elevation(position, time_list)
        return np.sum(wave_elevations, axis=0)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        wave_elevations = np.zeros((len(self.list_of_waves), len(list_of_point)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index] = each_wave.get_elevation_at_nodes(list_of_point, global_time)
        return np.sum(wave_elevations, axis=0)

    ## velocity
    def get_velocity_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) | Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations \n
        :return: Get a list of elevations at one position with a time sequence \n
        """
        waves_velocities = np.zeros((len(self.list_of_waves), len(time_list), 3))
        for index, each_wave in enumerate(self.list_of_waves):
            waves_velocities[index] = each_wave.get_velocity_with_time(position, time_list, irregularwaves=True)
        return np.sum(waves_velocities, axis=0)

    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_velocity_set = np.zeros((len(self.list_of_waves), len(list_of_point), 3))
        list_elevation = self.get_elevation_at_nodes(list_of_point, global_time)
        points = list_of_point.copy()
        # wheeler stretching method
        points[:, 2] = (list_of_point[:, 2] - list_elevation) * self.water_depth / (self.water_depth + list_elevation)
        # linear stretching method
        # points[:, 2][points[:, 2] >= 0] = 0
        for index, each_wave in enumerate(self.list_of_waves):
            node_velocity_set[index] = each_wave.get_velocity_at_nodes(points, global_time, irregularwaves=True)
        velocities = np.sum(node_velocity_set, axis=0)
        # ensure the velocity is zero above the wave elevation
        velocities[list_of_point[:, 2] > list_elevation] = 0
        return velocities

    def get_acceleration_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of node positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_acceleration_set = np.zeros((len(self.list_of_waves), len(list_of_point), 3))
        list_elevation = self.get_elevation_at_nodes(list_of_point, global_time)
        points = list_of_point.copy()
        # wheeler stretching method
        points[:, 2] = (list_of_point[:, 2] - list_elevation) * self.water_depth / (self.water_depth + list_elevation)
        # linear stretching method
        # points[:, 2][points[:, 2] >= 0] = 0
        for index, each_wave in enumerate(self.list_of_waves):
            node_acceleration_set[index] = each_wave.get_acceleration_at_nodes(points, global_time, irregularwaves=True)
        accelerations = np.sum(node_acceleration_set, axis=0)
        accelerations[list_of_point[:, 2] > list_elevation] = 0
        return accelerations


class current:
    def __init__(self, current_profile: dict, nodes_on_top_rings: dict, nodes_on_net: dict, number_of_nodes: int,
                 Sn: float):
        """
        Parameters
        ----------------
        :param current_profile: dictionary.
        :param nodes_on_top_rings: dictionary.
        :param nodes_on_net: dictionary.
        :param number_of_nodes: int.
        """
        self.cd= 0.33 * Sn + 6.54 * pow(Sn, 2) - 4.88 * pow(Sn, 3)
        self.current_profile = current_profile
        self.ring_nodes = nodes_on_top_rings
        self.netting_nodes = nodes_on_net
        self.number_of_cages = len(nodes_on_top_rings.keys())
        self.number_of_nodes_per_ring = len(nodes_on_top_rings['cage_0'])
        self.number_of_nodes_per_netting = len(nodes_on_net['cage_0'])
        self.net2net_factors = [1.0] * number_of_nodes
        self.cage2cage_factors = [1.0] * number_of_nodes
        # The first current vector define the direction for cage-to-cage wake effect direction
        current_top = np.array(current_profile['velocity'][0])
        self.current_vector = current_top / np.linalg.norm(current_top)  # unit vector
        # bounding box for cage-to-cage wake effect
        self.boxes = {}
        for each in range(self.number_of_cages):
            self.boxes['cage_' + str(each)] = [0, 0, 0, 0, 0]  # x_0,y_0,z_min, z_max, dy

    def linear_interpolation_current(self, z):
        if z > np.max(np.array(self.current_profile['depth'])):
            u_final = self.current_profile['velocity'][0]
        elif z < np.min(np.array(self.current_profile['depth'])):
            u_final = self.current_profile['velocity'][-1]
        else:
            index = np.argmin(np.abs(np.array(self.current_profile['depth']) - z))
            u1 = self.current_profile['velocity'][index]
            z1 = self.current_profile['depth'][index]
            if z > z1:
                z2 = self.current_profile['depth'][index - 1]
                u2 = self.current_profile['velocity'][index - 1]
            else:
                z2 = self.current_profile['depth'][index + 1]
                u2 = self.current_profile['velocity'][index + 1]
            dz = abs(float(z1 - z2))
            u_final = np.array(u2) * float(abs(z - z1) / dz) + np.array(u1) * abs(z - z2) / dz
        return u_final

    def update_bounding_boxes(self, node_position):
        """
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :return: Only update the variable without return value.
        """
        for i in range(self.number_of_cages):
            top_ring_center = node_position[self.ring_nodes['cage_' + str(i)]].mean(axis=0)
            z_min = node_position[self.netting_nodes['cage_' + str(i)]][:, 2].min()
            z_max = node_position[self.netting_nodes['cage_' + str(i)]][:, 2].max()
            dy = node_position[self.netting_nodes['cage_' + str(i)]][:, 1].ptp()
            # x_0,y_0,z_min, z_max, dy
            self.boxes['cage_' + str(i)] = [top_ring_center[0], top_ring_center[1], z_min, z_max, dy]
    def initial_wake_factor(self,node_position):
        """

        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :return:
        """
        self.get_net2net_wake(node_position)

    def get_net2net_wake(self, node_position):
        for i in range(self.number_of_cages):
            origin = np.array(self.boxes['cage_' + str(i)][:2])
            for each in self.netting_nodes['cage_' + str(i)]:
                vector_to_origin = node_position[each][:2] - origin
                if np.dot(vector_to_origin, self.current_vector[:2]) >= 1.0:
                    self.net2net_factors[each] = 1 - 0.46 * self.cd
                else:
                    self.net2net_factors[each] = 1.0




    def get_cage2cage_wake(self, node_position):

        def factor(x,y):
            """
            :param x: local x
            :param y: local y
            :return:
            """

            a0 = 0.02402
            a1 = 0.04824
            a2 = 0.002303
            a3 = -0.01288
            a4 = 0.0006093
            a5 = 0.005875
            a6 = -0.001147
            a7 = -0.002984
            w = 2.692

            val = a0 + a1 * np.cos(y * w) + a2 * np.cos(2 * x * w) + a3 * np.cos(3 * x * w) + a4 * np.cos(
            4 * x * w) + a5 * np.cos(
            5 * x * w) + a6 * np.cos(6 * x * w) + a7 * np.cos(7 * x * w)

            return  val * 5
        # outer loop is by each nodes
        for i in range(self.number_of_cages):
            origin = np.array(self.boxes['cage_' + str(i)][:2])
            dx = self.boxes['cage_' + str(i)][-1]*25  # wake length is 25 times of fish cage diameter
            dy = self.boxes['cage_' + str(i)][-1]*2 # wake width is 2 times of fish cage diameter
            dz = self.boxes['cage_' + str(i)][-2]-self.boxes[-3]


            for each in self.netting_nodes['cage_' + str(i)]:
                vector_to_origin = node_position[each][:2] - origin
                local_x=np.dot(vector_to_origin, self.current_vector[:2])
                local_y = pow((pow(np.linalg.norm(vector_to_origin),2)-pow(local_x,2)),0.5)
                local_z=node_position[each][-1]
                if ( 1.5*dy/4<=local_x<=dx  and local_y <=dy and local_z):


                    reducx = val * m_averageSn / 0.25;
                    self.net2net_factors[each] = 1 - 0.46 * self.cd
                else:
                    self.net2net_factors[each] = 1.0

        pass



    def update_wake(self, node_position):
        """
        update the cage to cage wake factor and net to net wake factor
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :return: Only update the variable without return value.
        """
        self.update_bounding_boxes(node_position)
        self.get_cage2cage_wake(node_position)
        print(self.boxes)

    def get_velocity_at_nodes(self, node_position):
        node_velocity = np.zeros((len(node_position), 3))
        self.update_wake(node_position)
        for index, position in enumerate(node_position.tolist()):
            node_velocity[index] = self.linear_interpolation_current(position[2]) * self.net2net_factors[index] * \
                                   self.cage2cage_factors[index]

        return node_velocity


if __name__ == "__main__":
    pass
    k = np.array([[0, 0, 0], [1.2, 1.2, 1.2], [0, 1, 2], [10, 12, 1]])
    node = [0, 1, 3]
    print(k[node][:, 1].ptp())
