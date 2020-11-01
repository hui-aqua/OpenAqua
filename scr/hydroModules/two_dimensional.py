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


class screenModel:
    """
    For Screen hydrodynamic models, the forces on netting are calculated based on individual a panel section of netting.
    The twines and knots in the net panel are considered as an integrated structure. In this module, the net panel is defined by
    three nodes because three (non-collinear) points can determine a unique plane in Euclidean geometry.
    In practice, the force is usually decomposed into two components: drag force F_D and lift force F_L (Cheng et al., 2020).
    """

    def __init__(self, model_index, hydro_element, solidity, dw0=0.0, dwh=0.0):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'S1', 'S2', 'S3'.
        :param hydro_element: [[list]] Unit: [-]. A python list to indicate how the net panel are connected. e.g.:[[p1,p2,p3][p2,p3,p4,p5]...]. If the input net panel contains 4 nodes, it will automaticly decomposed to 3 node net panel.
        :param solidity: [float] Unit: [-]. The solidity of netting.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
        :param dwh: [float] Unit: [m]. The hydrodynamic diameter of the numerical net twines. It is used for the force calculation (reference area)
        """
        self.modelIndex = str(model_index)
        self.triangular_elements = []
        self.dw0 = dw0
        self.dwh = dwh
        self.dws = np.sqrt(dwh / dw0) * dw0
        self.sn = solidity
        self.FEtime = 0
        converted_index = hydro_element.copy()
        # converted_index = list(hydro_element) #if is python2
        for index, item in enumerate(converted_index):
            for j in range(len(item)):
                converted_index[index][j] -= 1
        for panel in converted_index:  # loop based on the hydrodynamic elements
            if len([int(k) for k in set(panel)]) <= 3:  # the hydrodynamic element is a triangle
                self.triangular_elements.append(
                    [k for k in set([int(k) for k in set(panel)])])  # a list of the node sequence
            else:
                for i in range(len(panel)):
                    nodes = [int(k) for k in panel]  # get the list of nodes [p1,p2,p3,p4]
                    nodes.pop(i)  # delete the i node to make the square to a triangle
                    self.triangular_elements.append(nodes)  # delete the i node to make the square to a triangle
        self.hydro_dynamic_forces = np.zeros((len(self.triangular_elements), 3))
        self.hydro_static_forces = np.zeros((len(self.triangular_elements), 3))
        self.hydro_total_forces = np.zeros((len(self.triangular_elements), 3))

    def __str__(self):
        """Print information of the present object."""
        s0 = "Screen model"
        s1 = "The model index is " + str(self.modelIndex) + "\n"
        s2 = "In total, there are " + str(len(self.triangular_elements)) + " hydrodynamic triangular elements. \n"
        s3 = "The total force on the netting are \nFx=" + str(
            sum(self.hydro_dynamic_forces[:, 0])) + "N\n" + "Fy=" + str(
            sum(self.hydro_dynamic_forces[:, 2])) + "N\n" + "Fz=" + str(sum(self.hydro_dynamic_forces[:, 2])) + "N\n"
        return s0 + s1 + s2 + s3

    def output_hydro_element(self):
        """
        :return: [[list]] of the indexes of points in elements. e.g.:[[1,2,3],[2,3,4]...]
        """
        return self.triangular_elements

    def hydro_coefficients(self, point1, point2, point3, fluid_velocity, knot=False):
        """
        :return:
        :param inflow_angle: [float] Unit [rad]. Definition: the angle between normal vector of a net panel and the flow direction
        :param fluid_velocity: [np.array].shape=(1,3) Unit: [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :param knot: [boolean] knot option. *Default=False*
        :return: drag and lift force coefficients. [float] Unit: [-].
        """
        a1 = point2 - point1
        a2 = point3 - point1
        unit_normal_vector = np.cross(a1, a2) / np.linalg.norm(np.cross(a1, a2))
        if np.dot(unit_normal_vector, fluid_velocity) < 0:
            unit_normal_vector = -unit_normal_vector
        surface_area = 0.5 * np.linalg.norm(np.cross(a1, a2))
        lift_vector = np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) / np.linalg.norm(
            np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) + np.finfo(np.float64).eps)
        # NOTE: avoid dividing by zero
        coin_alpha = np.dot(unit_normal_vector, fluid_velocity) / np.linalg.norm(fluid_velocity)
        inflow_angle = np.arccos(coin_alpha)

        drag_coefficient, lift_coefficient = 0, 0
        if self.modelIndex not in 'S1,S2,S3,S4,S5,S6':
            print("The selected hydrodynamic model is not included in the present program")
            exit()
        elif self.modelIndex == 'S1':  # aarsnes 1990
            drag_coefficient = 0.04 + (-0.04 + self.sn - 1.24 * pow(self.sn, 2) + 13.7 * pow(self.sn, 3)) * np.cos(
                inflow_angle)
            lift_coefficient = (0.57 * self.sn - 3.54 * pow(self.sn, 2) + 10.1 * pow(self.sn, 3)) * np.sin(
                2 * inflow_angle)

        elif self.modelIndex == 'S2':  # Loland 1991
            drag_coefficient = 0.04 + (
                    -0.04 + 0.33 * self.sn + 6.54 * pow(self.sn, 2) - 4.88 * pow(self.sn, 3)) * np.cos(
                inflow_angle)
            lift_coefficient = (-0.05 * self.sn + 2.3 * pow(self.sn, 2) - 1.76 * pow(self.sn, 3)) * np.sin(
                2 * inflow_angle)

        elif self.modelIndex == 'S3':  # Kristiansen 2012
            a1 = 0.9
            a3 = 0.1
            b2 = 1.0
            b4 = 0.1
            reynolds_number = row_water * self.dw0 * np.linalg.norm(fluid_velocity) / dynamic_viscosity / (
                    1 - self.sn)  # Re
            cd_cylinder = -78.46675 + 254.73873 * np.log10(reynolds_number) - 327.8864 * pow(np.log10(reynolds_number),
                                                                                             2) + 223.64577 * pow(
                np.log10(reynolds_number), 3) - 87.92234 * pow(
                np.log10(reynolds_number), 4) + 20.00769 * pow(np.log10(reynolds_number), 5) - 2.44894 * pow(
                np.log10(reynolds_number), 6) + 0.12479 * pow(np.log10(reynolds_number), 7)
            cd_zero = cd_cylinder * (self.sn * (2 - self.sn)) / (2.0 * pow((1 - self.sn), 2))
            cn_pi_4 = 0.5 * cd_cylinder * self.sn / pow(1 - self.sn, 2)
            cl_pi_4 = (0.5 * cd_zero - np.pi * cn_pi_4 / (8 + cn_pi_4)) / np.sqrt(2)
            drag_coefficient = cd_zero * (a1 * np.cos(inflow_angle) + a3 * np.cos(3 * inflow_angle))
            lift_coefficient = cl_pi_4 * (b2 * np.sin(2 * inflow_angle) + b4 * np.sin(4 * inflow_angle))

        elif self.modelIndex == 'S4':  # Fridman 1973
            reynolds_number = np.linalg.norm(fluid_velocity) * self.dw0 * row_water / dynamic_viscosity
            reynolds_star = reynolds_number / (2 * self.sn)
            coe_tangent = 0.1 * pow(reynolds_number, 0.14) * self.sn
            coe_normal = 3 * pow(reynolds_star, -0.07) * self.sn
            drag_coefficient = coe_normal * np.cos(inflow_angle) * pow(np.cos(inflow_angle), 2) + coe_tangent * np.sin(
                inflow_angle) * pow(
                np.sin(inflow_angle), 2)
            lift_coefficient = coe_normal * np.sin(inflow_angle) * pow(np.cos(inflow_angle), 2) + coe_tangent * np.cos(
                inflow_angle) * pow(
                np.sin(inflow_angle), 2)
        elif self.modelIndex == 'S5':  # Lee 2005 # polynomial fitting
            drag_coefficient = 0.556 * pow(inflow_angle, 7) - 1.435 * pow(inflow_angle, 6) - 2.403 * pow(
                inflow_angle, 5) + 11.75 * pow(inflow_angle, 4) - 13.48 * pow(inflow_angle, 3) + 5.079 * pow(
                inflow_angle, 2) - 0.9431 * pow(inflow_angle, 1) + 1.155
            lift_coefficient = -10.22 * pow(inflow_angle, 9) + 69.22 * pow(inflow_angle, 8) - 187.9 * pow(
                inflow_angle, 7) + 257.3 * pow(inflow_angle, 6) - 181.6 * pow(inflow_angle, 5) + 59.14 * pow(
                inflow_angle, 4) - 7.97 * pow(inflow_angle, 3) + 2.103 * pow(
                inflow_angle, 2) + 0.2325 * pow(inflow_angle, 1) + 0.01294

        elif self.modelIndex == 'S6':  # Balash 2009
            reynolds_cylinder = row_water * self.dw0 * np.linalg.norm(fluid_velocity) / dynamic_viscosity / (
                    1 - self.sn) + 0.000001
            cd_cylinder = 1 + 10.0 / (pow(reynolds_cylinder, 2.0 / 3.0))
            drag_coefficient = cd_cylinder * (0.12 - 0.74 * self.sn + 8.03 * pow(self.sn, 2)) * pow(inflow_angle, 3)
            if knot:
                mesh_size = 10 * self.dw0  # assume Sn=0.2
                diameter_knot = 2 * self.dw0  # assume the knot is twice of the diameter of the twine
                reynolds_sphere = row_water * diameter_knot * np.linalg.norm(fluid_velocity) / dynamic_viscosity / (
                        1 - self.sn) + 0.000001
                coe_sphere = 24.0 / reynolds_sphere + 6.0 / (1 + np.sqrt(reynolds_sphere)) + 0.4
                drag_coefficient = (cd_cylinder * 8 * pow(diameter_knot,
                                                          2) + coe_sphere * np.pi * mesh_size * self.dw0) / np.pi * mesh_size * self.dw0 * (
                                           0.12 - 0.74 * self.sn + 8.03 * pow(self.sn, 2)) * pow(inflow_angle, 3)
            else:
                pass
        elif self.modelIndex == 'bi2014':  # from Table 1 in Bi et al., 2014, JFS
            p1 = 0.1873
            p2 = 0.4921
            p3 = 0.04
            drag_coefficient = p3 + p2 * np.cos(inflow_angle) + p1 * pow(np.cos(inflow_angle), 2)
            p1 = -0.169
            p2 = 0.4159
            p3 = 0
            lift_coefficient = p3 + p2 * np.sin(2 * inflow_angle) + p1 * pow(np.sin(2 * inflow_angle), 2)

        return surface_area, drag_coefficient, lift_coefficient, unit_normal_vector, lift_vector

    def force_on_element(self, node_position, velocity_fluid, velocity_structure=np.zeros((99999, 3)), elevation=False):
        """
        calculate hydrodynamic forces on triangular-type structure.
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :param velocity_fluid: np.array[n,3] | Unit [m/s]| fluid velocity at the coordinates of nodes, n is the number of nodes
        :param velocity_structure: np.array[n,3] | Unit [m/s]| structure velocity of nodes, n is the number of nodes
        :return: np.array[m,3] | Unit [N]| hydrodynamic forces on elements, m is the number of element
        """
        num_element = len(self.triangular_elements)
        force_on_element = np.zeros((num_element, 3))  # force on line element, initial as zeros
        for index, panel in enumerate(self.triangular_elements):  # loop based on the hydrodynamic elements
            p1 = node_position[panel[0]]
            p2 = node_position[panel[1]]
            p3 = node_position[panel[2]]
            element_u = (velocity_fluid[panel[0]] + velocity_fluid[panel[1]] + velocity_fluid[panel[2]]) / 3.0
            element_v = (velocity_structure[panel[0]] + velocity_structure[panel[1]] + velocity_structure[
                panel[2]]) / 3.0
            relative_velocity = element_u - element_v
            net_area, c_d, c_l, drag_direction, lift_direction = self.hydro_coefficients(p1, p2, p3, relative_velocity)
            if elevation is not False:
                element_center = (p1 + p2 + p3) / 3.0
                element_elevation = (elevation[int(panel[0])] + elevation[int(panel[1])] + elevation[
                    int(panel[2])]) / 3.0
                if element_center[2] > element_elevation:
                    row = row_air
                else:
                    row = row_water
            fd = 0.5 * row * net_area * c_d * pow(np.linalg.norm(relative_velocity), 2) * drag_direction
            fl = 0.5 * row * net_area * c_l * pow(np.linalg.norm(relative_velocity), 2) * lift_direction
            force_on_element[index] = (fd + fl) / 2.0
        self.hydro_dynamic_forces = np.array(force_on_element)
        return np.array(force_on_element)

    def screen_fsi(self, node_position, velocity_fluid, velocity_structure=np.zeros((99999, 3))):
        """
        calculate hydrodynamic forces on triangular-type structure. Usage for coupling simulation
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :param velocity_fluid: np.array[m,3] | Unit [m/s]| fluid velocity [ux,uy,uz] on net panels centroid, m is the number of element
        :param velocity_structure: np.array[n,3] | Unit [m/s]| structure velocity of nodes, n is the number of nodes
        :return: np.array[m,3] | Unit [N]| hydrodynamic forces on elements, m is the number of element
        """
        num_element = len(self.triangular_elements)
        force_on_element = np.zeros((num_element, 3))  # force on line element, initial as zeros
        if len(velocity_fluid) < len(self.triangular_elements):
            print("position is " + str(node_position))
            print("Velocity is " + str(velocity_structure))
            print("velocity elements is " + str(velocity_fluid))
            exit()
        for index, panel in enumerate(self.triangular_elements):  # loop based on the hydrodynamic elements
            p1 = node_position[panel[0]]
            p2 = node_position[panel[1]]
            p3 = node_position[panel[2]]
            element_u = velocity_fluid[index]
            element_v = (velocity_structure[panel[0]] + velocity_structure[panel[1]] + velocity_structure[
                panel[2]]) / 3.0
            relative_velocity = element_u - element_v
            net_area, c_d, c_l, drag_direction, lift_direction = self.hydro_coefficients(p1, p2, p3, relative_velocity)
            relative_velocity = velocity_fluid[index] * np.sqrt(2.0 / (2.0 - c_d - c_l))
            fd = 0.5 * row_water * net_area * c_d * pow(np.linalg.norm(relative_velocity), 2) * drag_direction
            fl = 0.5 * row_water * net_area * c_l * pow(np.linalg.norm(relative_velocity), 2) * lift_direction
            force_on_element[index] = (fd + fl) / 2.0
        self.hydro_dynamic_forces = np.array(force_on_element)
        return np.array(force_on_element)

    def cal_buoy_force(self, node_position, elevation):
        """
        calculate hydrodynamic forces on line-type structure.
        :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
        :param elevation: np.array[n,1] |Unit [m]| elevation of the sea, n is the number of nodes
        :return: np.array[m,3] |Unit [N]| buoyancy force on elements, m is the number of elements
        """
        num_element = len(self.triangular_elements)
        force_on_element = np.zeros((num_element, 3))  # force on line element, initial as zeros
        for index, panel in enumerate(self.triangular_elements):  # loop based on the hydrodynamic elements
            p1 = node_position[panel[0]]
            p2 = node_position[panel[1]]
            p3 = node_position[panel[2]]
            element_area = self.hydro_coefficients(p1, p2, p3, np.array([1, 0, 0]))[0]  # get the area
            element_volume = element_area * self.sn * self.dws * 0.25 * pi

            ae_vector = elevation[int(self.triangular_elements[index][0])] - node_position[int(self.triangular_elements[index][0])]
            be_vector = elevation[int(self.triangular_elements[index][1])] - node_position[int(self.triangular_elements[index][1])]
            ce_vector = elevation[int(self.triangular_elements[index][2])] - node_position[int(self.triangular_elements[index][2])]
            list_z = [ae_vector[2], be_vector[2], ce_vector[2]]
            sign_set = set(np.sign(list_z))

            if sign_set == {-1} or sign_set == {0, -1}:
                force_on_element[index] = [0, 0, element_volume * gravity * row_air]
            elif sign_set == {1} or sign_set == {0, 1}:
                force_on_element[index] = [0, 0, element_volume * gravity * row_water]
            elif sign_set == {0}:
                force_on_element[index] = [0, 0, element_volume * gravity * 0.5 * (row_water + row_air)]
            elif sign_set == {1, -1}:
                # ++- or --+
                if np.sign(list_z).sum() > 0:
                    # ++-
                    z1 = -min(np.array(list_z))
                    ratio = (z1 / abs(np.array(list_z))).prod()
                    force_on_element[index] = [0, 0,
                                               element_volume * gravity * (row_air * ratio + row_water * (1 - ratio))]
                else:  # --+
                    pass
                    z1 = max(np.array(list_z))
                    ratio = (z1 / abs(np.array(list_z))).prod()
                    force_on_element[index] = [0, 0,
                                               element_volume * gravity * (row_water * ratio + row_air * (1 - ratio))]
            else:
                # 0 + -
                total_height = abs(np.array(list_z)).sum()
                ratio = max(np.array(list_z)) / total_height
                force_on_element[index] = [0, 0, element_volume * gravity * (row_air * (1 - ratio) + row_water * ratio)]
        self.hydro_static_forces = np.array(force_on_element) / 2
        return np.array(force_on_element) / 2

    def screen_buoy_fsi(self, node_position, elevation):
        # todo later
        pass

    def distribute_force(self, number_of_node):
        """
        Transfer the forces on line element to their corresponding nodes.\n
        :return: [np.array].shape=(N,3) Unit [N]. The hydrodynamic forces on all N nodes.
        """
        force_on_nodes = np.zeros((number_of_node, 3))  # force on nodes, initial as zeros
        self.hydro_total_forces = self.hydro_static_forces + self.hydro_dynamic_forces
        for index, panel in enumerate(self.triangular_elements):
            force_on_nodes[panel[0]] += (self.hydro_total_forces[index]) / 3
            force_on_nodes[panel[1]] += (self.hydro_total_forces[index]) / 3
            force_on_nodes[panel[2]] += (self.hydro_total_forces[index]) / 3
        return force_on_nodes


if __name__ == "__main__":
    pass
