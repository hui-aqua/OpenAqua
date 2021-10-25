'''
/--------------------------------\
|    University of Stavanger     |
|           Hui Cheng            |
\--------------------------------/
Any questions about this code, please email: hui.cheng@uis.no
The center of the floating collar is (0,0,0)
Fish cage is along the Z- direction
Z=0 is the free surface
Z<0 is the water zone
Z>0 is the air zone
The sinkers are attached to the floating collar
'''


import os
import sys
import json
import numpy as np
from numpy import pi

# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter template start
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# default value
parameters = {
    "CageShape":
        {
            "origin": [[0, 0, 0]],
            "cageCircumference": 100.0,
            "bottom_top_ratio": 0.8,
            "cageHeight": 10.0,
            "cageConeHeight": 3,
            "element_length": 2.5
        },
    "Frame":
        {
            "topringDensity": 958,
            "topring_sec_dia": 0.35,
            "topring_thickness": 0.0185,
            "rope_depth_height_ratio": 2,
            "rope_sec_dia": 0.05,
            "weight_tip": 981,
            "weight_per_metter": 40.001
        },
}

# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter template finish
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


floater_center = parameters['CageShape']['origin']
cr_top = float(parameters['CageShape']['cageCircumference'] / pi / 2.0)
cr_bottom = cr_top * float(parameters['CageShape']['bottom_top_ratio'])
dr = cr_top - cr_bottom  #
cage_height = float(parameters['CageShape']['cageHeight'])
cage_cone_height = float(parameters['CageShape']['cageConeHeight'])

NT = int(parameters['CageShape']['cageCircumference'] / parameters['CageShape'][
    'element_length'])  # Number of the nodes in circumference
NN = int(cage_height / parameters['CageShape'][
    'element_length'])  # number of section in the height, thus, the nodes should be NN+1
BN = int(cr_bottom / parameters['CageShape']['element_length'] / np.sqrt(
    2))  # number of section along the cone, thus, the nodes should be NN+1

# below is the final information for outputing
# con_pipe
# con_rope
# con_netting
# sur_netting
# point


# generate the point coordinates matrix for cylindrical netting
# procedure for generating point_netting:
# 1. Generate the node on one cage based on rotating the nodes on one side using transcope matrix.
# 2. Generate the node on (all) netting through translation.
point_one_cage = []
# Step 1.
for i in range(0, NT):
    for j in range(NN):
        point_one_cage.append(
            [(cr_top - dr * j / float(NN)) * np.cos(i * 2 * pi / float(NT)),
             (cr_top - dr * j / float(NN)) * np.sin(i * 2 * pi / float(NT)),
             - j * cage_height / float(NN)])
    for j in range(BN):
        point_one_cage.append(
            [cr_bottom * ((BN - j) / BN) * np.cos(i * 2 * pi / float(NT)),
             cr_bottom * ((BN - j) / BN) * np.sin(i * 2 * pi / float(NT)),
             - cage_height - j * (cage_cone_height) / float(BN)])
point_one_cage.append([0, 0, -cage_cone_height - cage_height])  # the last point should be at the cone tip

number_of_point_one_netting = len(point_one_cage)
# print(number_of_point_one_netting)

# generate the point coordinates matrix for mooring part
r_top = cr_top + 1
NN_frame = NN
BN_frame = NN
NT_frame = NT
r_bottom = cr_bottom + 1

for i in range(0, NT_frame):
    for j in range(0, NN_frame):
        point_one_cage.append(
            [(r_top - dr * j / float(NN_frame)) * np.cos(i * 2 * pi / float(NT_frame)),
             (r_top - dr * j / float(NN_frame)) * np.sin(i * 2 * pi / float(NT_frame)),
             -cage_height * j / float(NN_frame)])
    for j in range(BN_frame):
        point_one_cage.append(
            [r_bottom * ((BN_frame - j) / BN_frame) * np.cos(i * 2 * pi / float(NT)),
             r_bottom * ((BN_frame - j) / BN_frame) * np.sin(i * 2 * pi / float(NT)),
             - cage_height - j * (float(parameters['Frame']['rope_depth_height_ratio']) * cage_height) / float(
                 BN_frame)])
point_one_cage.append([0, 0, -float(parameters['Frame'][
                                        'rope_depth_height_ratio']) * cage_height - cage_height])  # the last point should be at the cone tip

number_of_point_one_frame = len(point_one_cage) - number_of_point_one_netting
number_of_point_one_cage = len(point_one_cage)
# print(number_of_point_one_cage)
# Step 2
point = []
for center in floater_center:
    point += (np.array(point_one_cage) + np.ones((len(point_one_cage), 3)) * np.array(center)).tolist()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Generate connection and sur
# similar with the nodes, first is for one side, then for one cage, finally for all the cages.
# 1.one cage

con_in_one_netting = []
sur_in_one_netting = []
sur_in_one_netting_tri = []
# horizontal con for netting
for i in range(BN + NN):
    con_in_one_netting.append([1 + i, 1 + i + (BN + NN) * (NT - 1)])
    for j in range(NT - 1):
        con_in_one_netting.append([1 + i + j * (BN + NN), 1 + i + (j + 1) * (BN + NN)])

# vertical con for netting
for i in range(BN + NN - 1):
    for j in range(NT):
        con_in_one_netting.append([1 + i + j * (BN + NN), 1 + 1 + i + j * (BN + NN)])
        if j < NT - 1:
            sur_in_one_netting.append([1 + i + j * (BN + NN), 1 + 1 + i + j * (BN + NN),
                                       1 + i + (j + 1) * (BN + NN), 1 + 1 + i + (j + 1) * (BN + NN)])
        else:
            pass
            sur_in_one_netting.append([1 + i + j * (BN + NN), 1 + 1 + i + j * (BN + NN),
                                       1 + i, 1 + 1 + i])
        # print(sur_in_one_cage)                #test
# bottom part, triangular element
for j in range(NT):
    con_in_one_netting.append([number_of_point_one_netting,
                               (j + 1) * (BN + NN)])
    if j < NT - 1:
        sur_in_one_netting_tri.append([number_of_point_one_netting,
                                       (j + 1) * (BN + NN),
                                       (j + 1 + 1) * (BN + NN)
                                       ])
    else:
        sur_in_one_netting_tri.append([number_of_point_one_netting,
                                       BN + NN,
                                       NT * (BN + NN)
                                       ])

# generate the connection for frame
# print(len(point_netting))             #test
# print(number_of_point_one_frame)           #test
# con_in_one_frame = []
con_in_one_frame_top = []
con_in_one_frame_line = []

for j in range(NT_frame):
    # top
    con_in_one_frame_top.append([1 + j * (BN + NN), 1 + j * (NN_frame + BN_frame) + number_of_point_one_netting])
    # circular
    if j < NT_frame - 1:
        con_in_one_frame_top.append([1 + j * (NN_frame + BN_frame) + number_of_point_one_netting,
                                     1 + (j + 1) * (NN_frame + BN_frame) + number_of_point_one_netting])
    else:
        con_in_one_frame_top.append(
            [1 + j * (NN_frame + BN_frame) + number_of_point_one_netting, 1 + number_of_point_one_netting])
        # vertical lines
    for i in range(NN_frame + BN_frame - 1):
        con_in_one_frame_line.append([1 + i + j * (NN_frame + BN_frame) + number_of_point_one_netting,
                                      1 + 1 + i + j * (NN_frame + BN_frame) + number_of_point_one_netting])
    con_in_one_frame_line.append(
        [NN_frame + BN_frame + j * (NN_frame + BN_frame) + number_of_point_one_netting, number_of_point_one_cage])
    # TODO
    for num in range(int(NN / 2)):
        con_in_one_frame_line.append([1 + NN + j * (BN + NN) - num * 2,
                                      1 + NN_frame + j * (NN_frame + BN_frame) + number_of_point_one_netting - 2 * num])

# Step 2
sur_netting_tri = []
con_one_cage = con_in_one_netting + con_in_one_frame_line + con_in_one_frame_top
num_of_con_in_one_cage = len(con_one_cage)

con_all = []
con_pipe_top = []
con_rope = []
con_netting = []  # lines on netting
sur_netting = []

for num in range(len(floater_center)):
    sur_netting += (np.array(sur_in_one_netting) + num * number_of_point_one_cage).tolist()
    sur_netting_tri += (np.array(sur_in_one_netting_tri) + num * number_of_point_one_cage).tolist()
    con_all += (np.array(con_one_cage) + num * number_of_point_one_cage).tolist()
    con_pipe_top += (np.array(con_in_one_frame_top) + num * number_of_point_one_cage).tolist()
    con_rope += (np.array(con_in_one_frame_line) + num * number_of_point_one_cage).tolist()
    con_netting += (np.array(con_in_one_netting) + num * number_of_point_one_cage).tolist()
sur_netting += sur_netting_tri



if __name__ == '__main__':
    cwd = os.getcwd()
    elements = {}

    for index, item in enumerate(sur_netting):
        for j in range(len(item)):
            sur_netting[index][j] -= 1

    elements['netting'] = sur_netting
    elements['pipes_top'] = (np.array(con_pipe_top, dtype=int)-1).tolist()
    elements['rope'] = (np.array(con_rope, dtype=int)-1).tolist()
    with open(os.path.join(cwd, 'element.json'), "w") as json_file:
        json.dump(elements, json_file)
    json_file.close()

    point_position={}
    point_position['initial']=point
    with open(os.path.join(cwd, 'Nodes_position.json'), "w") as json_file:
        json.dump(point_position, json_file)
    json_file.close()

    print("Number of node is " + str(len(point)))
    print("Number of netting element is " + str(len(sur_netting)))
    print("Number of rope element is " + str(len(con_rope)))
    print("Number of pipe element is " + str(len(con_pipe_top)))