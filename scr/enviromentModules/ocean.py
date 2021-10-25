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
    def __init__(self, current_profile):
        """
        Parameters
        ----------------
        :param current_profile: dictionary or list
        * If it is a dictionary, this means the flow velocity is varied at different water depth,
        and it has a format as follows:
        current_profile={'depth': [0, -10, -30, -60],
                         'velocity': [[0.5, 0, 0],
                                      [0.6, 0, 0],
                                      [0.3, 0, 0],
                                      [0.0, 0, 0]]
                        }

        * If it is a list, this means the flow velocity is uniform along water depth, and it has a format as:
        current_profile=[0.5,0,0]
        """
        self.current_profile = current_profile
        self.flow_velocity = []

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

    def update_velocity_at_nodes(self, node_position):
        if type(self.current_profile) == type([]):
            self.flow_velocity = np.ones_like(node_position) * self.current_profile
        elif type(self.current_profile) == type({}):
            self.flow_velocity = np.zeros_like(node_position)
            for index, position in enumerate(node_position.tolist()):
                self.flow_velocity[index] = self.linear_interpolation_current(position[2])

    def get_velocity_at_nodes(self, node_position):
        self.update_velocity_at_nodes(node_position)
        return self.flow_velocity


class wakeEffect:
    def __init__(self, current_direction: list, net2netWake=False, cage2cageWake=False):
        """

        :param current_direction:  [ux,uy,uz] Unit: m/s, give the direction of current direction
        :param net2netWake: True or False, trigger for net-to-net wake effect
        :param cage2cageWake: True or False, trigger for cage-to-cage wake effect
        """
        self.direction = np.array(current_direction) / np.linalg.norm(np.array(current_direction))
        self.wake1 = net2netWake
        self.wake2 = cage2cageWake

    def initial_wake(self, node_position: list, fish_cages: dict, Sn: float):
        """

        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :param fish_cages:
        :param Sn:
        :return:
        """

        def list2set(list_of_list):
            outer_list = []

            def lists_to_list(list_of_list):
                for el in list_of_list:
                    if type(el) == list:
                        lists_to_list(el)
                    else:
                        outer_list.append(el)
                return set(outer_list)

            return lists_to_list(list_of_list)

        self.Sn = Sn
        self.cd = 0.33 * Sn + 6.54 * pow(Sn, 2) - 4.88 * pow(Sn, 3)
        self.ru1 = [1.0] * len(node_position)
        self.ru2 = [1.0] * len(node_position)
        self.number_of_cages = len(fish_cages.keys())

        self.cage_tri_elem = {}
        self.cage_nodes = {}
        self.boxes = {}
        self.width={}
        for key in fish_cages.keys():
            self.cage_tri_elem[key] = set(fish_cages[key])
            self.cage_nodes[key] = list2set(fish_cages[key])
            # x_0,y_0,z_min, z_max, D
            self.boxes[key] = [0, 0, 0, 0]
            self.width[key] = 0.5 * (node_position[self.cage_nodes[key]][:, 1].ptp() +
                                     node_position[self.cage_nodes[key]][:,0].ptp()
                                    )

        if self.wake1:
            self.update_ru_net2net(node_position)
        else:
            pass

    def update_ru_net2net(self, node_position):
        """
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :return: Only update the variable without return value.
        """
        for i in range(self.number_of_cages):
            origin = np.array(self.boxes['cage_' + str(i)][:2])
            for each in self.cage_tri_elem['cage_' + str(i)]:
                vector_from_origin = node_position[each][:2] - origin
                if np.dot(vector_from_origin, self.direction[:2]) >= 0.5:
                    self.ru1[each] = 1 - 0.46 * self.cd
                else:
                    self.ru1[each] = 1.0

    def update_boxes(self, node_position):
        for key in self.boxes.keys():
            nodes = node_position[self.cage_nodes[key]]
            top_ring_center = nodes.mean(axis=0)
            z_min = nodes[:, 2].min()
            z_max = nodes[:, 2].max()
            # x_0,y_0,z_min, z_max,
            self.boxes[key] = [top_ring_center[0], top_ring_center[1], z_min, z_max]

    def update_ru_cage2cage(self, node_position):
        """
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :return: Only update the variable without return value.
        """

        def factor(x, y):
            """
            :param x: local x  [m]
            :param y: local y  [m]
            :return: cage to cage flow reduction factor
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
            # ref: experiments
            val = a0 + a1 * np.cos(y * w) + a2 * np.cos(2 * y * w) + a3 * np.cos(3 * y * w) + a4 * np.cos(
                4 * y * w) + a5 * np.cos(5 * y * w) + a6 * np.cos(6 * y * w) + a7 * np.cos(7 * y * w)

            return val * 5 * self.Sn / 0.25 * pow(np.exp(-(x - 1.5) / 25.0), 0.5)

        self.update_boxes(node_position)

        # outer loop is by each nodes
        for index in range(len(node_position)):
            for key in self.boxes.keys():
                origin = np.array(self.boxes[key][:2])
                diameter = self.width[key]
                if index not in self.cage_nodes[key]:
                    vector_from_origin = node_position[index][:2] - origin
                    dz = self.boxes[key][3] - self.boxes[key][2]  # max-min
                    local_x = (vector_from_origin @ self.direction[:2]) / diameter
                    local_y = pow(pow(np.linalg.norm(vector_from_origin) / diameter, 2) - pow(local_x, 2), 0.5)
                    local_z = node_position[index][-1] - self.boxes[key][2]
                    if (1 <= local_x <= 25 and local_y <= 1 and 0 < local_z < dz):
                        self.ru2[index] *= factor(local_x, local_y)
                    else:
                        self.ru2[index] *= 1.0

    def get_reduction_factor(self, node_position):
        """
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :return: Only update the variable without return value.
        """
        self.update_ru_cage2cage(node_position)
        return np.array(self.ru1) * np.array(self.ru2)

    def update_u(self, flow_velocity, node_position):
        """

        :param flow_velocity:
        :param node_position:
        :return:
        """
        self.update_ru_cage2cage(node_position)
        u = np.zeros_like(node_position)
        for i in range(len(flow_velocity)):
            u[i] = flow_velocity[i] * self.ru1[i] * self.ru2[i]
        return u


if __name__ == "__main__":
    U = {'depth': [0, -10, -30, -60],
         'velocity': [[0.5, 0, 0],
                      [0.6, 0, 0],
                      [0.3, 0, 0],
                      [0.0, 0, 0]]
         }
    env = current(U)
    uc = env.get_velocity_at_nodes(np.array([[0, 0, 0],
                                             [0.5, 10, 10],
                                             [0, 0, -12],
                                             [0, 0, -32],
                                             [0, 10, -52],
                                             [0, 10, -72]]))
    print(uc)


    def list2set(list_of_list):
        outer_list = []

        def lists_to_list(list_of_list):
            for el in list_of_list:
                if type(el) == list:
                    lists_to_list(el)
                else:
                    outer_list.append(el)
            return set(outer_list)

        return lists_to_list(list_of_list)


    ee = [[2, 3], [2, 0]]
    ee2 = [55, [23, 2], [2, 5, 18], [158, 5, 6, 4]]
    t = list(list2set(ee))
    t2 = list2set(ee2)
    print(uc[t])
