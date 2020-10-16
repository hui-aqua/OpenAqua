"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
Modules can be used to calculate the hydrodynamic forces on nets.
In order to use this module, we recommend ``import screenModel as sm`` in the beginning of your code.
Please refer to Cheng et al. (2020) [https://doi.org/10.1016/j.aquaeng.2020.102070] for details.
# In the present program, the nodes' position is stored as numpy array not a python list for easier manipulation.
# Velocity of current and/or wave is also a numpy array.
# Element index is stored as a python list for fast running.

"""

import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

row = 1025  # [kg/m3]   sea water density
kinematic_viscosity = 1.004e-6  # when the water temperature is 20 degree.
dynamic_viscosity = 1.002e-3  # when the water temperature is 20 degree.

class net2netWeak:
    """
    A module that can be use to deal with net to net wake effect.\n

    .. note: Can only apply to a single fish cage.
    """

    def __init__(self, model_index, node_initial_position, hydro_element, current_velocity, origin, dw0, net_solidity):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'factor-0.9', 'loland-0.5', 'hui-1'.
        :param node_initial_position: [np.array].shape=(N,3) Unit: [m]. The initial coordinates of N nodes in cartesian coordinate system.
        :param hydro_element: [list] Unit: [-]. A python list to indicate how the lines are connected.
        :param current_velocity: [np.array].shape=(1,3) or a [list] of three values Unit: [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :param origin: [np.array].shape=(1,3) or a [list] of three values Unit: [m]. The origin [x,y,z] for detecting the elements in the wake region. For a fish cage, the origin is usually sit in the floating collar.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines.
        :param net_solidity: [float] Unit: [-]. The solidity of netting.
        """
        self.positions = np.array(node_initial_position)
        self.elements = index_convert(hydro_element)
        self.flow_direction = np.array(current_velocity) / np.linalg.norm(current_velocity)
        # np.array for incoming velocity to a cage
        self.origin = np.array(origin) - np.array(current_velocity) / np.linalg.norm(current_velocity) * (
                2 * dw0 / net_solidity)
        # The coordinate of origin, it can be cage center [0,0,0]
        # (2 * dw0 / net_solidity) is a half mesh size, it is a safety factor.
        self.Sn = net_solidity
        self.wake_type = str(model_index).split("-")[0]
        self.wake_value = str(model_index).split("-")[1]
        self.wake_element_indexes = self.get_element_in_wake()
        print("\n net2net weak effect is initialized.\n")

    def __str__(self):
        """Print information of the present object."""
        s0 = "The selected wake model is " + str(self.wake_type) + "\n"
        s1 = "The index of the element in the wake region is " + str(self.wake_element_indexes) + "\n"
        S = s0 + s1
        return S

    def is_element_in_wake(self, one_element):
        """
        A private function to tell if an element is in wake region or not.\n
        :param one_element: [List] of three values Unit [-]. A list to contain the index of nodes in one element. The list should contain at least two values.
        :return: True if the element in the wave region.
        """
        element_center = np.array([0.0, 0.0, 0.0])  # must be a float
        for node in one_element:
            element_center += np.array(self.positions[int(node)] / len(one_element))
        vector_point_to_origin = np.array(element_center - self.origin)
        if np.dot(vector_point_to_origin, self.flow_direction) < 0:
            return True
        else:
            return False

    def get_element_in_wake(self):
        """
        A private function to go go through all the elements and find the elements in the wake region.\n
        :return: [List] Unit [-]. A list of "indexes of the elements" in the wake region.
        """
        elements_in_wake = []
        for index,element in enumerate(self.elements):
            if self.is_element_in_wake(element):
                elements_in_wake.append(index)
        return elements_in_wake

    def reduction_factor(self, element_index, alpha=0):
        """
        A public function to return the flow reduction factor according to the index of element (and inflow angle).\n
        :param element_index: [integer] Unit [-]. A index of element
        :param alpha: [float] Unit [rad]. *Default value=0*; The inflow angle of the net panel.
        :return: [float] Unit [-]. *Range [0,1]* Flow reduction factor. If value=1, the flow is no change. If value is 0, the flow is reduced to zero.
        """
        factor = 1.0
        if element_index in self.wake_element_indexes:
            if self.wake_type in ['factor']:
                factor = min(1.0, float(self.wake_value))
            elif self.wake_type in ['loland']:
                factor = float(self.cal_factor1())
            elif self.wake_type in ['hui']:
                factor = float(self.cal_factor2(alpha))
            else:
                print("the selected wake type " + str(self.wake_type) + " is not supported.")
                print()
        else:
            factor = 1.0
        return factor

    def cal_factor1(self):
        """
        A private function calculate the flow reduction factor according to Loland (1991). r=1-0.14*CD (CD is the input value in the ```wakeModel```) \n
        :return: Flow reduction factor: [float] Unit [-].
        """
        return 1 - 0.46 * float(self.wake_value)

    def cal_factor2(self, alf):
        """
        A private function calculate the flow reduction factor according to Hui et al., (2020).\n
        :param alf: [float] Unit [rad]. Definition: the angle between normal vector of a net panel and the flow direction.
        :return: Flow reduction factor: [float] Unit [-].
        """
        alf = np.abs(alf)
        reduction_factor = (np.cos(alf) + 0.05 - 0.38 * self.Sn) / (np.cos(alf) + 0.05)
        return max(0, reduction_factor)


class morisonModel:
    """
    For Morison hydrodynamic models, the forces on netting are calculated based on individual twines.
    The twines are taken as cylindrical elements. In practice, the force is usually decomposed into two componnets:
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
        self.modelIndex = str(model_index)
        self.elements = index_convert(hydro_element)  # the connections of the twines[[p1,p2][p2,p3]]
        self.dwh = dwh  # used for the force calculation (reference area)
        self.dw0 = dw0  # used for the hydrodynamic coefficients
        self.Sn = solidity
        self.force_on_elements = np.zeros((len(self.elements), 3))

    def __str__(self):
        """Print information of the present object."""
        s0="Morison modelm\n"
        s1="The model index is "+ str(self.modelIndex)+ "\n"
        s2="In total, there are "+str(len(self.elements))+" hydrodynamic line elements. \n"
        s3="The total force on the nettings are \nFx="+str(sum(self.force_on_elements[:,0])) +"N\n" +  "Fy="+str(sum(self.force_on_elements[:,2])) +"N\n" +"Fz="+str(sum(self.force_on_elements[:,2])) +"N\n"
        return s0+s1+s2+s3

    def output_hydro_element(self):
        """
        :return: [list] Unit: [-]. A list of indexes of the elements in the wake region.
        """
        return self.elements

    def hydro_coefficients(self, current_velocity):
        """
        :param current_velocity: [np.array].shape=(1,3) Unit: [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :return: normal and tangential drag force coefficients. [float] Unit: [-].
        """
        drag_normal = 0
        drag_tangent = 0
        if self.modelIndex not in 'M1,M2,M3,M4,M5':
            print("The selected hydrodynamic model is not included in the present program")
            exit()
        elif self.modelIndex == 'M1':  # Bessonneau 1998
            drag_normal = 1.2
            drag_tangent = 0.1
        elif self.modelIndex == 'M2':  # Wan 2002
            drag_normal = 1.3
            drag_tangent = 0.0
        elif self.modelIndex == 'M3':  # Takagi 2004
            reynolds_number = row * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity
            if reynolds_number < 200:
                drag_normal = pow(10, 0.7) * pow(reynolds_number, -0.3)
            else:
                drag_normal = 1.2
            drag_tangent = 0.1
        elif self.modelIndex == 'M4':  # choo 1971
            reynolds_number = row * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity
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
        elif self.modelIndex == 'M5':  # cifuentes 2017
            reynolds_number = row * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity
            drag_normal = -3.2891e-5 * pow(reynolds_number * self.Sn * self.Sn,
                                           2) + 0.00068 * reynolds_number * self.Sn * self.Sn + 1.4253
        return drag_normal, drag_tangent

    def force_on_element(self, net_wake, realtime_node_position, current_velocity, wave=False, fe_time=0):
        """
        :param net_wake: A object wake model, net2net wake model. Must create first.
        :param realtime_node_position: [np.array].shape=(N,3) Unit: [m]. The coordinates of N nodes in cartesian coordinate system.
        :param current_velocity: [np.array].shape=(1,3) Unit [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :param wave:  A wake model object. *Default value=False* Must create first.
        :param fe_time: [float] Unit [s]. The time in Code_Aster. *Default value=0* Must give if wave is added.
        :param current_velocity: numpy array ([ux,uy,uz]), Unit [m/s]
        :return: [np.array].shape=(M,3) Unit [N]. The hydrodynamic forces on all M elements. Meanwhile, update the self.force_on_elements
        """
        num_line = len(self.elements)
        hydro_force_on_element = []  # force on line element, initial as zeros
        wave_velocity = np.zeros((num_line, 3))
        if wave:
            for index,line in enumerate(self.elements):
                element_center = (realtime_node_position[int(line[0])] + realtime_node_position[int(line[1])]) / 2.0
                wave_velocity[index] = wave.get_velocity(element_center, fe_time)

        for i in range(num_line):
            b = get_distance(realtime_node_position[int(self.elements[i][0])],
                             realtime_node_position[int(self.elements[i][1])])
            a = get_orientation(realtime_node_position[int(self.elements[i][0])],
                                realtime_node_position[int(self.elements[i][1])])

            velocity = np.array(current_velocity) * net_wake.reduction_factor(i) + wave_velocity[i]
            drag_n, drag_t = self.hydro_coefficients(velocity)
            ft = 0.5 * row * self.dwh * (b - self.dwh) * drag_t * np.dot(a, velocity) * a * np.linalg.norm(
                np.dot(a, velocity))
            fn = 0.5 * row * self.dwh * (b - self.dwh) * drag_n * (velocity - np.dot(a, velocity) * a) * np.linalg.norm(
                (velocity - np.dot(a, velocity) * a))
            hydro_force_on_element.append(ft + fn)
        self.force_on_elements = np.array(hydro_force_on_element)
        return np.array(hydro_force_on_element)

    def distribute_force(self, number_of_node):
        """
        Transfer the forces on line element to their corresponding nodes.\n
        :return: [np.array].shape=(N,3) Unit [N]. The hydrodynamic forces on all N nodes
        """
        force_on_nodes = np.zeros((number_of_node, 3))  # force on nodes, initial as zeros
        for index, line in enumerate(self.elements):
            force_on_nodes[line[0]] += (self.force_on_elements[index]) / 2
            force_on_nodes[line[1]] += (self.force_on_elements[index]) / 2
        return force_on_nodes

class screenModel:

    """
    For Morison hydrodynamic models, the forces on netting are calculated based on individual a panel section of netting.
    The twines and knots in the net panel are considered as an integrated structure. In this module, the net panel is defined by
    three nodes because three (non-colinear) points can determine a unique plane in Euclidean geometry.
    In practice, the force is usually decomposed into two componnets: drag force F_D and lift force F_L (Cheng et al., 2020).
    """

    def __init__(self, model_index, hydro_element, solidity, dw0=0.0,dwh=0.0):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'S1', 'S2', 'S3'.
        :param hydro_element: [[list]] Unit: [-]. A python list to indicate how the net panel are connected. e.g.:[[p1,p2,p3][p2,p3,p4,p5]...]. If the input net panel contains 4 nodes, it will automaticly decomposed to 3 node net panel.
        :param solidity: [float] Unit: [-]. The solidity of netting.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
        :param dwh: [float] Unit: [m]. The hydrodynamic diameter of the numerical net twines. It is used for the force calculation (reference area)
        """
        self.modelIndex = str(model_index)
        self.hydro_element = convert_hydro_element(hydro_element)
        self.dwh = dwh
        self.dw0 = dw0
        self.Sn = solidity
        self.FEtime=0
        self.force_on_elements = np.zeros((len(self.hydro_element), 3))
    def __str__(self):
        """Print information of the present object."""
        s0="Screen model"
        s1="The model index is "+ str(self.modelIndex)+ "\n"
        s2="In total, there are "+str(len(self.hydro_element))+" hydrodynamic triangular elements. \n"
        s3="The total force on the nettings are \nFx="+str(sum(self.force_on_elements[:,0])) +"N\n" + "Fy="+str(sum(self.force_on_elements[:,2])) +"N\n" +"Fz="+str(sum(self.force_on_elements[:,2])) +"N\n"
        return s0+s1+s2+s3

    def output_hydro_element(self):
        """
        :return: [[list]] of the indexes of points in elements. e.g.:[[1,2,3],[2,3,4]...]
        """
        return self.hydro_element

    def hydro_coefficients(self, inflow_angle, current_velocity, knot=False):
        """
        :return:
        :param inflow_angle: [float] Unit [rad]. Definition: the angle between normal vector of a net panel and the flow direction
        :param current_velocity: [np.array].shape=(1,3) Unit: [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :param knot: [boolean] knot option. *Default=False*
        :return: drag and lift force coefficients. [float] Unit: [-].
        """
        drag_coefficient, lift_coefficient = 0, 0
        if self.modelIndex not in 'S1,S2,S3,S4,S5,S6':
            print("The selected hydrodynamic model is not included in the present program")
            exit()
        elif self.modelIndex == 'S1':  # aarsnes 1990
            drag_coefficient = 0.04 + (-0.04 + self.Sn - 1.24 * pow(self.Sn, 2) + 13.7 * pow(self.Sn, 3)) * np.cos(
                inflow_angle)
            lift_coefficient = (0.57 * self.Sn - 3.54 * pow(self.Sn, 2) + 10.1 * pow(self.Sn, 3)) * np.sin(
                2 * inflow_angle)

        elif self.modelIndex == 'S2':  # Loland 1991
            drag_coefficient = 0.04 + (
                    -0.04 + 0.33 * self.Sn + 6.54 * pow(self.Sn, 2) - 4.88 * pow(self.Sn, 3)) * np.cos(
                inflow_angle)
            lift_coefficient = (-0.05 * self.Sn + 2.3 * pow(self.Sn, 2) - 1.76 * pow(self.Sn, 3)) * np.sin(
                2 * inflow_angle)

        elif self.modelIndex == 'S3':  # Kristiansen 2012
            a1 = 0.9
            a3 = 0.1
            b2 = 1.0
            b4 = 0.1
            reynolds_number = row * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity / (
                    1 - self.Sn)  # Re
            cd_cylinder = -78.46675 + 254.73873 * np.log10(reynolds_number) - 327.8864 * pow(np.log10(reynolds_number),
                                                                                             2) + 223.64577 * pow(
                np.log10(reynolds_number), 3) - 87.92234 * pow(
                np.log10(reynolds_number), 4) + 20.00769 * pow(np.log10(reynolds_number), 5) - 2.44894 * pow(
                np.log10(reynolds_number), 6) + 0.12479 * pow(np.log10(reynolds_number), 7)
            cd_zero = cd_cylinder * (self.Sn * (2 - self.Sn)) / (2.0 * pow((1 - self.Sn), 2))
            cn_pi_4 = 0.5 * cd_cylinder * self.Sn / pow(1 - self.Sn, 2)
            cl_pi_4 = (0.5 * cd_zero - np.pi * cn_pi_4 / (8 + cn_pi_4)) / np.sqrt(2)
            drag_coefficient = cd_zero * (a1 * np.cos(inflow_angle) + a3 * np.cos(3 * inflow_angle))
            lift_coefficient = cl_pi_4 * (b2 * np.sin(2 * inflow_angle) + b4 * np.sin(4 * inflow_angle))

        elif self.modelIndex == 'S4':  # Fridman 1973
            reynolds_number = np.linalg.norm(current_velocity) * self.dw0 * row / dynamic_viscosity
            reynolds_star = reynolds_number / (2 * self.Sn)
            coe_tangent = 0.1 * pow(reynolds_number, 0.14) * self.Sn
            coe_normal = 3 * pow(reynolds_star, -0.07) * self.Sn
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
            reynolds_cylinder = row * self.dw0 * np.linalg.norm(current_velocity) / dynamic_viscosity / (
                    1 - self.Sn) + 0.000001
            cd_cylinder = 1 + 10.0 / (pow(reynolds_cylinder, 2.0 / 3.0))
            drag_coefficient = cd_cylinder * (0.12 - 0.74 * self.Sn + 8.03 * pow(self.Sn, 2)) * pow(inflow_angle, 3)
            if knot:
                mesh_size = 10 * self.dw0  # assume Sn=0.2
                diameter_knot = 2 * self.dw0  # assume the knot is twice of the diameter of the twine
                reynolds_sphere = row * diameter_knot * np.linalg.norm(current_velocity) / dynamic_viscosity / (
                        1 - self.Sn) + 0.000001
                coe_sphere = 24.0 / reynolds_sphere + 6.0 / (1 + np.sqrt(reynolds_sphere)) + 0.4
                drag_coefficient = (cd_cylinder * 8 * pow(diameter_knot,
                                                          2) + coe_sphere * np.pi * mesh_size * self.dw0) / np.pi * mesh_size * self.dw0 * (
                                           0.12 - 0.74 * self.Sn + 8.03 * pow(self.Sn, 2)) * pow(inflow_angle, 3)
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

        return drag_coefficient, lift_coefficient

    def force_on_element(self, net_wake, realtime_node_position, current_velocity, velocity_nodes=np.zeros((99999, 3)),
                         wave=False, fe_time=0):
        """
        :param velocity_nodes: [np.array].shape=(1,3) Unit [m/s]. The structural velocity of nodes [ux,uy,uz] in cartesian coordinate system.
        :param net_wake: A object wake model, net2net wake model. Must create first.
        :param realtime_node_position: [np.array].shape=(N,3) Unit: [m]. The coordinates of N nodes in cartesian coordinate system.
        :param current_velocity: [np.array].shape=(1,3) Unit [m/s]. The current velocity [ux,uy,uz] in cartesian coordinate system.
        :param wave:  A wake model object. *Default value=False* Must create first.
        :param fe_time: [float] Unit [s]. The time in Code_Aster. *Default value=0* Must give if wave is added.
        :param current_velocity: numpy array ([ux,uy,uz]), Unit [m/s]
        :return: [np.array].shape=(M,3) Unit [N]. The hydrodynamic forces on all M elements. Meanwhile, update the self.force_on_elements
        """
        self.FEtime=fe_time
        wave_velocity = np.zeros((len(self.hydro_element), 3))
        if wave:
            for index, panel in enumerate(self.hydro_element):
                element_center = (realtime_node_position[int(panel[0])] + realtime_node_position[int(panel[1])] +
                                  realtime_node_position[int(panel[2])]) / 3
                wave_velocity[index] = wave.get_velocity(element_center, fe_time)
        hydro_force_on_element = []  # force on net panel, initial as zeros
        for index, panel in enumerate(self.hydro_element):  # loop based on the hydrodynamic elements
            p1 = realtime_node_position[panel[0]]
            p2 = realtime_node_position[panel[1]]
            p3 = realtime_node_position[panel[2]]
            velocity_structure = (velocity_nodes[panel[0]] + velocity_nodes[panel[1]] + velocity_nodes[panel[2]]) / (
                len(panel))
            alpha, lift_direction, net_area = calculation_on_element(p1, p2, p3, np.array(current_velocity))
            # calculate the inflow angel, normal vector, lift force factor, area of the hydrodynamic element
            velocity = np.array(current_velocity) * net_wake.reduction_factor(index, alpha) + wave_velocity[index]
            # if not in the wake region, the effective velocity is the undisturbed velocity
            while np.linalg.norm(velocity_structure) > np.linalg.norm(velocity):
                velocity_structure *= 0.1
            if np.dot(velocity_structure, velocity) > 0:
                velocity_relative = velocity - velocity_structure
            else:
                velocity_relative = velocity - velocity_structure * 0
            drag_coefficient, lift_coefficient = self.hydro_coefficients(alpha, velocity_relative, knot=False)
            fd = 0.5 * row * net_area * drag_coefficient * np.linalg.norm(velocity_relative) * velocity_relative
            fl = 0.5 * row * net_area * lift_coefficient * pow(np.linalg.norm(velocity_relative), 2) * lift_direction
            hydro_force_on_element.append((fd + fl) / 2.0)
            # print("panel is "+str(panel))
            # print("area of panel is "+str(net_area))
            # print("density of water is "+str(row))
            # print("drag_coefficient is "+str(drag_coefficient))
            # print("lift_coefficient is "+str(lift_coefficient))
            # print("velocity_relative "+str(velocity_relative))
            # print("fd is "+str(fd))
            # print("fl is "+str(fl))
            # print("fh is "+str((fd + fl) / 2.0))

        if np.size(np.array(hydro_force_on_element)) == np.size(self.hydro_element):
            self.force_on_elements = np.array(hydro_force_on_element)
            return np.array(hydro_force_on_element)
        else:
            print("\nError! the size of hydrodynamic force on element is not equal to the number of element."
                  "\nPlease cheack you code.")
            print("\nThe size of element is " + str(len(self.hydro_element)))
            print("\nThe size of hydrodynamic force is " + str(len(np.array(hydro_force_on_element))))
            exit()

    def screen_fsi(self, realtime_node_position, velocity_on_element, velocity_of_nodes=np.zeros((9999, 3))):
        """
        :param realtime_node_position: [np.array].shape=(N,3) Unit: [m]. The coordinates of N nodes in cartesian coordinate system.
        :param velocity_on_element: [np.array].shape=(M,3) Unit [m/s]. The current velocity [ux,uy,uz] on all net panels in cartesian coordinate system.
        :param velocity_of_nodes: [np.array].shape=(N,3) Unit: [m/s]. The velocities of N nodes in cartesian coordinate system.
        :return: update the self.force_on_elements and output the forces on all the net panels.
        """
        # print("The length of U vector is " + str(len(velocity_on_element)))
        hydro_force_on_element = []  # force on net panel, initial as zeros
        # print("velocity elementy is 1" + str(velocity_on_element))
        # print("velocity_of_nodes is " + str(velocity_of_nodes))
        # print("realtime_node_position is " + str(realtime_node_position))
        if len(velocity_on_element) < len(self.hydro_element):
            print("position is " + str(realtime_node_position))
            print("Velocity is " + str(velocity_of_nodes))
            print("velocity elements is " + str(velocity_on_element))
            exit()
        for index, panel in enumerate(self.hydro_element):  # loop based on the hydrodynamic elements
            p1 = realtime_node_position[panel[0]]
            p2 = realtime_node_position[panel[1]]
            p3 = realtime_node_position[panel[2]]
            alpha, lift_direction, surface_area = calculation_on_element(p1, p2, p3, velocity_on_element[index])
            drag_coefficient, lift_coefficient = self.hydro_coefficients(alpha, velocity_on_element[index], knot=False)

            velocity_fluid = velocity_on_element[index] * np.sqrt(2.0 / (2.0 - drag_coefficient - lift_coefficient))
            velocity_structure = (velocity_of_nodes[panel[0]] + velocity_of_nodes[panel[1]] + velocity_of_nodes[
                panel[2]]) / (len(panel))
            while np.linalg.norm(velocity_structure) > np.linalg.norm(velocity_fluid):
                velocity_structure *= 0.1
            if np.dot(velocity_structure, velocity_fluid) > 0:
                velocity_relative = velocity_fluid - velocity_structure
            else:
                velocity_relative = velocity_fluid - velocity_structure * 0

            fd = 0.5 * row * surface_area * drag_coefficient * np.linalg.norm(np.array(velocity_relative)) * np.array(
                velocity_relative)
            fl = 0.5 * row * surface_area * lift_coefficient * pow(np.linalg.norm(velocity_relative),
                                                                   2) * lift_direction
            hydro_force_on_element.append((fd + fl) / 2.0)
        if np.size(np.array(hydro_force_on_element)) == np.size(self.hydro_element):
            self.force_on_elements = np.array(hydro_force_on_element)
            return np.array(hydro_force_on_element)
        else:
            print("Error!, the size of hydrodynamic force on element is not equal to the number of element."
                  "Please cheack you code.")
            exit()

    # def distribute_velocity(self, current_velocity, wave_velocity=0):
    #     """
    #     :param current_velocity:
    #     :param wave_velocity:
    #     :return:
    #     """
    #     velocity_on_element = []  # velocity on net panel, initial as zeros
    #     if len(current_velocity) < 4:
    #         velocity_on_element = np.ones((len(self.hydro_element), 3)) * current_velocity
    #     elif len(current_velocity) == len(self.hydro_element):
    #         velocity_on_element = np.array(current_velocity)
    #
    #     if not wave_velocity == 0:
    #         for panel in self.hydro_element:  # loop based on the hydrodynamic elements
    #             velocity_on_element[self.hydro_element.index(panel)] = +wave_velocity[self.hydro_element.index(panel)]
    #     if len(velocity_on_element) == len(self.hydro_element):
    #         return velocity_on_element
    #     else:
    #         print("\nError! the size of velocity on element is not equal to the number of element."
    #               "\nPlease cheack you code.")
    #         print("\nThe size of element is " + str(len(self.hydro_element)))
    #         print("\nThe size of hydrodynamic force is " + str(len(velocity_on_element)))
    #         exit()

    def distribute_force(self, number_of_node, force_increasing_factor=1):
        """
        Transfer the forces on line element to their corresponding nodes.\n
        :return: [np.array].shape=(N,3) Unit [N]. The hydrodynamic forces on all N nodes.
        """
        forces_on_nodes = np.zeros((number_of_node, 3))  # force on nodes, initial as zeros
        for index, panel in enumerate(self.hydro_element):
            forces_on_nodes[panel[0]] += force_increasing_factor * (self.force_on_elements[index]) / 3
            forces_on_nodes[panel[1]] += force_increasing_factor * (self.force_on_elements[index]) / 3
            forces_on_nodes[panel[2]] += force_increasing_factor * (self.force_on_elements[index]) / 3
        return forces_on_nodes



def calculation_on_element(point1, point2, point3, velocity):
    """
    Moudule private function.\n
    :param point1: point1 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :param point2: point2 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :param point3: point3 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :param velocity: [np.array].shape=(1,3) Unit: [m/s]. Flow velocity
    :return: inflow angle [float][rad], lift vector[np.array].shape=(1,3) [-], area of net panel [float][m^2]
    """
    # because the mesh construction, the first two node cannot have same index
    a1 = get_orientation(point1, point2)
    a2 = get_orientation(point1, point3)
    ba1 = get_distance(point1, point2)
    ba2 = get_distance(point1, point3)
    normal_vector = np.cross(a1, a2) / np.linalg.norm(np.cross(a1, a2))
    if np.dot(normal_vector, velocity) < 0:
        normal_vector = -normal_vector
    surface_area = 0.5 * np.linalg.norm(np.cross(a1 * ba1, a2 * ba2))
    lift_vector = np.cross(np.cross(velocity, normal_vector), velocity) / \
                  np.linalg.norm(np.cross(np.cross(velocity, normal_vector), velocity) + 0.000000001)

    coin_alpha = abs(np.dot(normal_vector, velocity) / np.linalg.norm(velocity))
    alpha = np.arccos(coin_alpha)
    return alpha, lift_vector, surface_area


def convert_hydro_element(elements):
    """
    :param elements: [[list]] Unit: [-]. A list of indexes of elements which can contain 4 or 3 nodes.
    :return: [[list]] Unit: [-]. A list of indexes of elements only contain 3 nodes.
    """
    hydro_elements = []
    for panel in index_convert(elements):  # loop based on the hydrodynamic elements
        if len([int(k) for k in set(panel)]) <= 3:  # the hydrodynamic element is a triangle
            hydro_elements.append([k for k in set([int(k) for k in set(panel)])])  # a list of the node sequence
        else:
            for i in range(len(panel)):
                nodes = [int(k) for k in panel]  # get the list of nodes [p1,p2,p3,p4]
                nodes.pop(i)  # delete the i node to make the square to a triangle
                hydro_elements.append(nodes)  # delete the i node to make the square to a triangle
    return hydro_elements


def get_distance(p1, p2):
    """
    Module private function.\n
    :param p1: point1 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :param p2: point2 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :return: The distance between two points.  [float] Unit [m].
    """

    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    return np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)


def get_orientation(p1, p2):
    """
    Module private function.\n
    :param p1: point1 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :param p2: point2 [np.array].shape=(1,3) or a [list] of coordinates Unit: [m].
    :return: The unit vector from p1 to p2.
    """
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    p = np.array([dx, dy, dz])
    return p / np.linalg.norm(p)


# two function used by code_aster
def get_position_aster(table_aster):
    """
    Module public function.\n
    :param table_aster: A table from Code_Aster by ``POST_RELEVE_T`` command with NOM_CHAM=('DEPL')
    :return:  [np.array].shape=(N,3) Unit: [m]. A numpy array of all the nodes positions.
    """
    content = table_aster.EXTR_TABLE()
    original_x = content.values()['COOR_X']
    original_y = content.values()['COOR_Y']
    original_z = content.values()['COOR_Z']
    delta_x = content.values()['DX']
    delta_y = content.values()['DY']
    delta_z = content.values()['DZ']
    position = np.array([original_x, original_y, original_z]) + np.array([delta_x, delta_y, delta_z])
    return np.transpose(position)


def get_velocity_aster(table_aster):  # to get the velocity
    """
    Module public function.\n
    :param table_aster: A table from Code_Aster by ``POST_RELEVE_T`` command with NOM_CHAM=('VITE')
    :return:  [np.array].shape=(N,3) Unit: [m/s]. A numpy array of all the nodes velocities.
    """
    content = table_aster.EXTR_TABLE()
    velocity_x = content.values()['DX']
    velocity_y = content.values()['DY']
    velocity_z = content.values()['DZ']
    velocity = np.array([velocity_x, velocity_y, velocity_z])
    return np.transpose(velocity)

def index_convert(input_connection):
    out=input_connection
    for i, con in enumerate(out):
        for j in range(len(con)):
            out[i][j]-=1
    return out



if __name__ == "__main__":
    # test code here
    pass
