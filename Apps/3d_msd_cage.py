import os
import json
import numpy as np
import scipy.integrate as inte


# ref to CCW-C100-H10-W40-U6

# ------------------------------------------------------------------
row_water=1025.0 #kg/m3
dynamic_viscosity = 1.002e-3  # when the water temperature is 20 degree.

def drag_coefficient(relative_velocity,inflow_angle):
    dw0=2.85e-3
    sn=0.2056
    a1 = 0.9
    a3 = 0.1
    b2 = 1.0
    b4 = 0.1
    # print(row_water * dw0 * np.linalg.norm(relative_velocity))
    reynolds_number = row_water * dw0 * np.linalg.norm(relative_velocity) / dynamic_viscosity / (1 - sn)  # Re

    cd_cylinder = -78.46675 + 254.73873 * np.log10(reynolds_number) - 327.8864 * pow(np.log10(reynolds_number),2) + 223.64577 * pow(
        np.log10(reynolds_number), 3) - 87.92234 * pow(
        np.log10(reynolds_number), 4) + 20.00769 * pow(np.log10(reynolds_number), 5) - 2.44894 * pow(
        np.log10(reynolds_number), 6) + 0.12479 * pow(np.log10(reynolds_number), 7)
    cd_zero = cd_cylinder * (sn * (2 - sn)) / (2.0 * pow((1 - sn), 2))
    cn_pi_4 = 0.5 * cd_cylinder * sn / pow(1 - sn, 2)
    cl_pi_4 = (0.5 * cd_zero - np.pi * cn_pi_4 / (8 + cn_pi_4)) / np.sqrt(2)
    drag_coefficient = cd_zero * (a1 * np.cos(inflow_angle) + a3 * np.cos(3 * inflow_angle))
    lift_coefficient = cl_pi_4 * (b2 * np.sin(2 * inflow_angle) + b4 * np.sin(4 * inflow_angle))
    return drag_coefficient,lift_coefficient

def drag_on_line(p1, p2, u1, u2):
    d = 0.5
    l = np.linalg.norm(np.array(p1) - np.array(p2))
    A = d * l

    u = np.array([1.0, 0.0, 0.0])
    v_m = (u1 + u2) / 2
    ur = np.array(u) - np.array(v_m)
    dx=abs(p1[0]-p2[0])
    dz = abs(p1[2] - p2[2])
    inflow_angle=np.arctan(dx/dz)
    # print(ur)
    cd,cl=drag_coefficient(ur,inflow_angle)
    drag_direction=np.array([1,0,0])
    lift_direction = np.array([0, 0, 1])

    return 0.5 * row_water * A* (cd  * pow(np.linalg.norm(ur), 2)*drag_direction+cl  * pow(np.linalg.norm(ur), 2)*lift_direction)


def tension(p1, p2):
    k = 1e6  # spring constant, N/m
    l0 = 0.1
    l = np.linalg.norm(np.array(p1) - np.array(p2))
    # print('extension is m '+str(l-l0))
    # print('tension is '+str(max(k * (l - l0), 0)))
    return max(k * (l - l0), 0)


def damping():

    return 0


def cal_ddx(point_position, point_velocity, element):
    m = 0.5  # mass, Kg
    c = 300  # damping coefficient
    gravity = np.array([0, 0, -9.81])
    force_on_node = np.zeros((len(point_position), 3))
    T_list=[]
    for index, line in enumerate(element):
        p0 = point_position[element[index][0]]
        p1 = point_position[element[index][1]]
        v0 = point_velocity[element[index][0]]
        v1 = point_velocity[element[index][1]]
        ei = (p1 - p0) / np.linalg.norm(np.array(p1) - np.array(p0))


        T = tension(p0, p1)  # scalar
        T_list.append(T)
        Fd = drag_on_line(p0, p1, v0, v1)

        force_on_node[line[0]] += + T * ei + m * gravity +Fd / 2
        force_on_node[line[1]] +=  - T * ei + m * gravity +Fd / 2

    return force_on_node / m -c*point_velocity


def MassSpringDamper(t, state):

    p = state.flatten()[:int(len(state) / 2)]
    v = state.flatten()[int(len(state) / 2):]
    p_p = np.reshape(p, (-1, 3))
    p_v = np.reshape(v, (-1, 3))

    ddx = cal_ddx(p_p, p_v, line)
    # print(type(v))
    # print(np.array(ddx).flatten().tolist())
    v[:3]=0
    # print(v[:3])
    return v.tolist() + np.array(ddx).flatten().tolist()  # [velocity,acceleration]


def get_data(odesolver, number_point):
    # collect data
    t_values = []
    x_values = []
    dx_values = []
    for i in range(500):
        # get solution step state
        odesolver.step()
        t_values.append(odesolver.t)
        x_values.append([odesolver.y[i * 3:(i + 1) * 3].tolist() for i in range(number_point)])
        dx_values.append(odesolver.y[number_point * 3:])
        # break loop after modeling is finished
        if odesolver.status == 'finished':
            break
    return t_values, x_values, dx_values


# Define the initial conditions and integration boundaries:


# ------------------------------------------------------------------


point = {}
point['p'] = [[0.000, 0.000, 0.00],
              [0.000, 0.000, -0.1],
              [0.000, 0.000, -0.2],
              [0.000, 0.000, -0.3],
              [0.000, 0.000, -0.4],
              [0.000, 0.000, -0.5],
              [0.000, 0.000, -0.6],
              [0.000, 0.000, -0.7],
              [0.000, 0.000, -0.8],
              [0.000, 0.000, -0.9]]

point['v'] = [[0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0]]
line = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5],[5,6],[6,7],[7,8],[8,9]]

initial_conditions = np.array(point['p']).flatten().tolist() + np.array(
    point['v']).flatten().tolist()  # [displacement, velocity]
#

t_end=5.0
num_step=100
solution=inte.solve_ivp(MassSpringDamper, t_span=[0, t_end], y0=initial_conditions,
                        method="BDF", t_eval=np.linspace(0,t_end,num_step))

Nodes_position = {}
cwd = os.getcwd()
for i in range(num_step):
    Nodes_position[i] = [solution.y[j * 3:(j + 1) * 3,i].tolist() for j in range(10)]
    # print('t='+str(i)+' s '+str(x[i]))
# #
with open(os.path.join(cwd, 'odeSolver/Nodes_position.json'), "w") as json_file:
    json.dump(Nodes_position, json_file)
json_file.close()

if __name__ == "__main__":
    pass
