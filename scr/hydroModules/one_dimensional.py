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
import sys
from numpy import pi

np.set_printoptions(threshold=sys.maxsize)
row_air = 1.225  # [kg/m3]   air density
row_water = 1025.0  # [kg/m3]   sea water density
kinematic_viscosity = 1.004e-6  # when the water temperature is 20 degree.
dynamic_viscosity = 1.002e-3  # when the water temperature is 20 degree.
gravity = 9.81


class MorisonModel:
    """
    For Morison hydrodynamic models, the forces on netting are calculated based on individual twines.
    The twines are taken as cylindrical elements. In practice, the force is usually decomposed into two components:
    normal drag force F_n and tangential drag force F_t (Cheng et al., 2020)
    """

    def __init__(self, model_index, hydro_element, solidity, dw0=0.0, dwh=0.0):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'M1', 'M2', 'M3'.
        :param hydro_element: [list] Unit: [-]. A python list to indicate how the lines are connected.
        :param solidity: [float] Unit: [-]. The solidity of netting.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
        :param dwh: [float] Unit: [m]. The hydrodynamic diameter of the numerical net twines. It is used for the force calculation (reference area)
        """
        self.model_index = str(model_index)
        self.line_elements = (np.array(hydro_element) - 1).tolist()
        self.dw0 = dw0  # used for the hydrodynamic coefficients
        self.dwh = dwh  # used for the force calculation (reference area)
        self.dws = np.sqrt(dwh / dw0) * dw0
        self.sn = solidity
        self.FEtime = 0
        self.hydro_dynamic_forces = np.zeros((len(hydro_element), 3))
        self.hydro_static_forces = np.zeros((len(hydro_element), 3))
        self.hydro_total_forces = np.zeros((len(hydro_element), 3))

    def __str__(self):
        """Print information of the present object."""
        s0 = "Morison model\n"
        s1 = "The model index is " + str(self.model_index) + "\n"
        s2 = "In total, there are " + str(len(self.line_elements)) + " hydrodynamic line elements. \n"
        s3 = "The total force on the netting are \nFx=" + str(sum(self.hydro_total_forces[:, 0])) + "N\n" + "Fy=" + str(
            sum(self.hydro_total_forces[:, 2])) + "N\n" + "Fz=" + str(sum(self.hydro_total_forces[:, 2])) + "N\n"
        return s0 + s1 + s2 + s3

    def output_hydro_element(self):
        """
        :return: [list] Unit: [-]. A list of indexes of the elements in the wake region.
        """
        return self.line_elements

    def hydro_coefficients(self, current_velocity):
        """
        :param current_velocity: [np.array].shape=(1,3) Unit: [m/s]. Current velocity [ux,uy,uz] in cartesian coordinate system.
        :return: normal and tangential drag force coefficients. [float] Unit: [-].
        """
        drag_normal = 0
        drag_tangent = 0
        if self.model_index not in 'M1,M2,M3,M4,M5':
            print("The selected hydrodynamic model is not included in the present program")
            exit()
        elif self.model_index == 'M1':  # Bessonneau 1998
            drag_normal = 1.2
            drag_tangent = 0.1
        elif self.model_index == 'M2':  # Wan 2002
            drag_normal = 1.3
            drag_tangent = 0.0
        elif self.model_index == 'M3':  # Takagi 2004
            reynolds_number = row_water * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity
            if reynolds_number < 200:
                drag_normal = pow(10, 0.7) * pow(reynolds_number, -0.3)
            else:
                drag_normal = 1.2
            drag_tangent = 0.1
        elif self.model_index == 'M4':  # choo 1971
            reynolds_number = row_water * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity
            drag_tangent = np.pi * dynamic_viscosity * (
                    0.55 * np.sqrt(reynolds_number) + 0.084 * pow(reynolds_number, 2.0 / 3.0))
            s = -0.07721565 + np.log(8.0 / reynolds_number)
            if 0 < reynolds_number < 1:
                drag_normal = 8 * np.pi * (1 - 0.87 * pow(s, -2)) / (s * reynolds_number)
            elif reynolds_number < 30:
                drag_normal = 1.45 + 8.55 * pow(reynolds_number, -0.9)
            elif reynolds_number < 2.33e5:
                drag_normal = 1.1 + 4 * pow(reynolds_number, -0.5)
            elif reynolds_number < 4.92e5:
                drag_normal = (-3.41e-6) * (reynolds_number - 5.78e5)
            elif reynolds_number < 1e7:
                drag_normal = 0.401 * (1 - np.exp(-reynolds_number / 5.99 * 1e5))
            else:
                print("Reynold number=" + str(reynolds_number) + ", and it exceeds the range.")
                exit()
        elif self.model_index == 'M5':  # cifuentes 2017
            reynolds_number = row_water * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity
            drag_normal = -3.2891e-5 * pow(reynolds_number * self.sn * self.sn,
                                           2) + 0.00068 * reynolds_number * self.sn * self.sn + 1.4253
        return drag_normal, drag_tangent

    def force_on_element(self, node_position, velocity_fluid, velocity_structure=np.zeros((99999, 3))):
        """
        calculate hydrodynamic forces on line-type structure.
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :param velocity_fluid: np.array[n,3] | Unit [m/s]|  fluid velocity at the coordinates of nodes, n is the number of nodes
        :param velocity_structure: np.array[n,3] | Unit [m/s]| structure velocity of nodes, n is the number of nodes
        :return: np.array[m,3] | Unit [N]| hydrodynamic forces on elements, m is the number of element
        """
        num_line = len(self.line_elements)

        force_on_element = np.zeros((num_line, 3))  # force on line element, initial as zeros
        for index in range(num_line):
            element_length = np.linalg.norm(
                node_position[int(self.line_elements[index][0])] - node_position[int(self.line_elements[index][1])])
            element_direction = (node_position[int(self.line_elements[index][0])] - node_position[
                int(self.line_elements[index][1])]) / element_length
            element_u = (velocity_fluid[int(self.line_elements[index][0])] + velocity_fluid[
                int(self.line_elements[index][1])]) / 2.0
            element_v = (velocity_structure[int(self.line_elements[index][0])] + velocity_structure[
                int(self.line_elements[index][1])]) / 2.0
            relative_velocity = element_u - element_v
            drag_n, drag_t = self.hydro_coefficients(relative_velocity)
            ft = 0.5 * row_water * self.dwh * (element_length - self.dwh) * drag_t * np.dot(element_direction,
                                                                                            relative_velocity) * element_direction * np.linalg.norm(
                np.dot(element_direction, relative_velocity))
            fn = 0.5 * row_water * self.dwh * (element_length - self.dwh) * drag_n * (
                        relative_velocity - np.dot(element_direction,
                                                   relative_velocity) * element_direction) * np.linalg.norm(
                (relative_velocity - np.dot(element_direction, relative_velocity) * element_direction))
            force_on_element[index] = ft + fn
        self.hydro_dynamic_forces = np.array(force_on_element)
        return np.array(force_on_element)

    def cal_buoy_force(self, node_position, elevation):
        """
        calculate hydrodynamic forces on line-type structure.
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :param elevation: np.array[n,1] |Unit [m]| elevation of the sea, n is the number of nodes
        :return: np.array[m,3] |Unit [N]| buoyancy force on elements, m is the number of elements
        """
        num_line = len(self.line_elements)
        force_on_element = np.zeros((num_line, 3))  # force on line element, initial as zeros
        for index in range(num_line):
            element_vector = node_position[int(self.line_elements[index][0])] - node_position[
                int(self.line_elements[index][1])]
            element_length = np.linalg.norm(element_vector)
            element_volume = 0.25 * pi * pow(self.dws, 2) * element_length

            ae_vector = elevation[int(self.line_elements[index][0])] - node_position[int(self.line_elements[index][0])]
            be_vector = elevation[int(self.line_elements[index][1])] - node_position[int(self.line_elements[index][1])]
            list_z = [ae_vector[2], be_vector[2]]
            sign_set = set(np.sign(list_z))

            if sign_set == {0}:
                force_on_element[index] = [0, 0, element_volume * gravity * 0.5 * (row_air + row_water)]
            elif sign_set == {1} or sign_set == {0, 1}:
                force_on_element[index] = [0, 0, element_volume * gravity * 0.5 * (row_water + row_water)]
            elif sign_set == {-1} or sign_set == {0, -1}:
                force_on_element[index] = [0, 0, element_volume * gravity * 0.5 * (row_air + row_air)]
            else:
                # partly submerged in water + -
                vertical_distance = abs(np.array(list_z)).sum()
                if vertical_distance > self.dws:
                    ratio = max(np.array(list_z)) / vertical_distance
                    force_on_element[index] = [0, 0,
                                               element_volume * gravity * (row_water * ratio + row_air * (1 - ratio))]
                # we assume it is horizontal to the water level
                else:
                    print("The vertical distance is too small. thus, the line is token as horizontal")
                    alpha = np.arccos(np.mean(list_z) / (0.5 * self.dws))
                    section_area = 0.25 * pi * pow(self.dws, 2)
                    area_fg = 0.25 * pow(self.dws, 2) * alpha - 0.5 * np.mean(list_z) * self.dws * np.sin(alpha)
                    submerged_area = 0.5 * section_area + np.sign(np.mean(list_z)) * (0.5 * section_area - area_fg)
                    submerged_volume = submerged_area * element_length
                    force_on_element[index] = [0, 0, submerged_volume * gravity * row_water]

        self.hydro_static_forces = np.array(force_on_element)
        return np.array(force_on_element)

    def distribute_force(self, number_of_node):
        """
        Transfer the forces on line element to their corresponding nodes.\n
        :return: [np.array].shape=(N,3) Unit [N]. The hydrodynamic forces on all N nodes
        """
        force_on_nodes = np.zeros((number_of_node, 3))  # force on nodes, initial as zeros
        self.hydro_total_forces = self.hydro_static_forces + self.hydro_dynamic_forces
        for index, line in enumerate(self.line_elements):
            force_on_nodes[line[0]] += (self.hydro_total_forces[index]) / 2
            force_on_nodes[line[1]] += (self.hydro_total_forces[index]) / 2
        return force_on_nodes


if __name__ == "__main__":
    pass
