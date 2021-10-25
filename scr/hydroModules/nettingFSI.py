"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
In order to use this module, we recommend ``import nettingFSI as fsi`` in the beginning of your code.
To use this module, one should be import this into the input file for calculations.

"""
import os
import time
import sys
import numpy as np
import json

np.set_printoptions(threshold=sys.maxsize)


def screen_fsi(self,
               node_position,
               velocity_fluid,
               velocity_structure=np.zeros((99999, 3))):
    """
    calculate hydrodynamic forces on triangular-type structure. Usage for coupling simulation
    :param node_position: np.array[n,3] | Unit [m]| coordinates of nodes, n is the number of nodes
    :param velocity_fluid: np.array[m,3] | Unit [m/s]| fluid velocity [ux,uy,uz] on net panels centroid, m is the number of element
    :param velocity_structure: np.array[n,3] | Unit [m/s]| structure velocity of nodes, n is the number of nodes
    :return: np.array[m,3] | Unit [N]| hydrodynamic forces on elements, m is the number of element
    """
    num_line = len(self.line_elements)
    force_on_element = np.zeros((num_line, 3), dtype=float)
    # loop based on the hydrodynamic elements
    if len(velocity_fluid) < len(self.hydro_element):
        print("position is " + str(node_position))
        print("Velocity is " + str(velocity_structure))
        print("velocity elements is " + str(velocity_fluid))
        exit()
    # loop based on the hydrodynamic elements
    for index, element in enumerate(self.hydro_element):
        p1 = node_position[element[0]]
        p2 = node_position[element[1]]
        p3 = node_position[element[2]]
        element_u = velocity_fluid[index]
        element_v = np.array([0.0] * 3)
        for each in element:
            element_v += velocity_fluid[each] / float(len(element))
        relative_velocity = element_u - element_v
        net_area, c_d, c_l, drag_direction, lift_direction = self.hydro_coefficients(
            p1, p2, p3, relative_velocity)
        relative_velocity = velocity_fluid[index] * \
                            np.sqrt(2.0 / (2.0 - c_d - c_l))
        fd = 0.5 * row_water * net_area * c_d * \
             pow(np.linalg.norm(relative_velocity), 2) * drag_direction
        fl = 0.5 * row_water * net_area * c_l * \
             pow(np.linalg.norm(relative_velocity), 2) * lift_direction
        force_on_element[index] = (fd + fl) / 2.0
    self.hydro_dynamic_forces = np.array(force_on_element)
    return np.array(force_on_element)








# here we assume Code_aster is much faster than OpenFoam, thus OpenFOAM do not need to wait .

def start_flag(cwd, flag):
    """
    :param cwd: work path
    :param flag: file name [string]
    :return:  create a empty file to tell openfoam the file starts writing
    """
    if os.path.isfile(os.path.join(cwd,str(flag))):
        print("path is :"+str(os.path.join(cwd,str(flag))))
        os.remove(os.path.join(cwd,str(flag)))
    else:
        pass


def finish_flag(cwd, flag):
    """
    :param cwd: work path
    :param flag: file name [string]
    :return:  create a empty file to tell openfoam the file finishes writing
    """
    # print("in finish_flag, the cwd is "+cwd)
    os.mknod(os.path.join(cwd,str(flag)))


def write_position(position, cwd):
    """
    :param position: A numpy array of all the nodes' position
    :param cwd: work path,
    :return: write the nodes' positions to "constant" folder as a file named "posi"
    """
    # print(cwd)
    start_flag(cwd, "position_flag")
    head_file = ["// Input for the nets in openfoam.",
                 "// Author: Hui Cheng",
                 "// Email: hui.cheng@uis.no",
                 "FoamFile",
                 "{",
                 "    version     2.0;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      posi;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "numOfPoint   " + str(len(position)) + ";"]
    with open(os.path.join(cwd, 'posi.tmp'), 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(position):
            output_file.write(
                "p" + str(index) + "\t( " + str(round(item[0],5)) + "\t" + str(round(item[1],5)) + "\t" + str(round(item[2],5)) + " );\n")
        output_file.write("}\n")
        output_file.write("// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //\n")
    output_file.close()

    # os.remove(os.path.join(cwd, 'posi'))
    try:
        os.rename(os.path.join(cwd, 'posi.tmp'), os.path.join(cwd, 'posi'))
    except FileExistsError:
        os.remove(os.path.join(cwd, 'posi'))
        os.rename(os.path.join(cwd, 'posi.tmp'), os.path.join(cwd, 'posi'))
    finish_flag(cwd, "position_flag")


def write_element(hydro_element, cwd):
    """
    :param hydro_element: A numpy array of all the net panel element
    :param cwd: work path,
    :return: write the net panels to "constant" folder as a file named "surf"
    """
    head_file = ["// Input for the nets in openfoam.",
                 "// Author: Hui Cheng",
                 "// Email: hui.cheng@uis.no",
                 "FoamFile",
                 "{",
                 "    version     2.0;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      surc;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "numOfSurf   " + str(len(hydro_element)) + ";"]
    with open(os.path.join(cwd,'../constant','surf') , 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(hydro_element):
            output_file.write(
                "e" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")
        output_file.write("// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //\n")
    output_file.close()


def write_fh(hydro_force, time_fe, cwd):
    """
    :param time_fe: simulation in FEM
    :param hydro_force: A numpy array of all the hydrodynamic forces on net panel
    :param cwd: work path,
    :return: write the hydrodynamic forces to "constant" folder as a file named "Fh" and save the total hydrodynamic forces on netting to "forceOnNetting.txt"
    """
    print("Here>>>>>>>>>>>>>Fh>>>  " + str(time_fe))
    start_flag(cwd, "fh_flag")

    head_file = ["// Input for the nets in openfoam.",
                 "// Author: Hui Cheng",
                 "// Email: hui.cheng@uis.no",
                 "FoamFile",
                 "{",
                 "    version     2.0;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      Fh;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "timeInFE   " + str(round(time_fe, 3)) + ";",
                 "numOfFh    " + str(len(hydro_force)) + ";"]
    with open(os.path.join(cwd, 'Fh.tmp'), 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")

        for index, item in enumerate(hydro_force):

            output_file.write(
                "fh" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")
        output_file.write("// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //\n")
    output_file.close()

    try:
        os.rename(os.path.join(cwd, 'Fh.tmp'), os.path.join(cwd, 'Fh'))
    except FileExistsError:
        os.remove(os.path.join(cwd, 'Fh'))
        os.rename(os.path.join(cwd, 'Fh.tmp'), os.path.join(cwd, 'Fh'))

    # write total force
    total_force = [round(time_fe, 3),
                   round(hydro_force[:, 0].sum(),3),
                   round(hydro_force[:, 1].sum(),3),
                   round(hydro_force[:, 2].sum(),3)]
    with open(os.path.join(cwd,'forceOnNetting.txt'), 'a+') as output_file:
        output_file.write(str(total_force) + '\n')
    output_file.close()

    finish_flag(cwd, "fh_flag")

velocity_dict = {'time_record': ['0']}


def get_velocity(cwd, length_velocity, time_aster):
    """
    :param cwd: working path for code aster
    :param length_velocity: the length of the list of element(or vlocity)
    :param time_aster: the time in Code_Aster
    :return: a numpy array of velocity on elements
    """
    print("Here>>>>>>>>>>>>>velo>>>  " + str(time_aster))

    velocity_file=os.path.join(cwd,'../velocity_on_elements.txt')
    velocity_dict, time_foam = read_velocity(velocity_file, length_velocity)
    time_foam = velocity_dict['time_record'][-1]
    # print(velocity_dict)
    # print(length_velocity)

    # exit()
    while float(time_aster) > float(time_foam) or str(time_foam) == str(0):
        time.sleep(0.1)
        velocity_dict, time_foam = read_velocity(velocity_file, length_velocity)
        # print(">>>>>>>>>>>>>>>.rean")
    else:
        with open(cwd + "/velocity_on_element.txt", "w") as f:
            json.dump(velocity_dict,f)
        f.close()
        return np.array(velocity_dict["velocities_at_" + str(time_foam)])


def read_velocity(file, length_velocity):
    """
    :param cwd_foam_root: work path
    :param length_velocity: The length of velocity file, normally equal to the length of elements
    :return: return the velocity dictionary and time in openfoam
    """
    while not os.path.exists(file):
        # print("path= "+ str(os.path.join(cwd_foam_root,"velocity_on_elements.txt")))
        print("Wait for velocity from OpenFoam......")
        time.sleep(0.1)

    else:
        with open(file, "r") as f:
            lines = f.readlines()
        f.close()
        # print(lines)
        time_foam = str(0)
        for line in lines[-length_velocity * 5:]:
            if "The velocities at" in line:

                time_foam = line.split(" ")[3]
                # print(time_foam)
            # else:
                # if len(line.split()) == 6:
                start_line = lines.index(line) + 3
                end_line = lines.index(line) + 3 + length_velocity
                if (end_line < len(lines)) and (str(time_foam) not in velocity_dict['time_record']):
                    velocity_dict['time_record'].append(str(time_foam))
                    velocity_dict["velocities_at_" + str(time_foam)] = []
                    for item in range(length_velocity):
                        try:
                            velocity_dict["velocities_at_" + str(time_foam)].append(
                                [float(lines[item + start_line].split()[0][1:]),
                                 float(lines[item + start_line].split()[1]),
                                 float(lines[item + start_line].split()[2][:-1])])
                        except:
                            velocity_dict["velocities_at_" + str(time_foam)].append(
                                [float(lines[start_line].split()[0][1:]),
                                 float(lines[start_line].split()[1]),
                                 float(lines[start_line].split()[2][:-1])])
    return velocity_dict, velocity_dict['time_record'][-1]

if __name__ == "__main__":
    pass
