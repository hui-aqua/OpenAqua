import os
import json
import numpy as np
import scipy.integrate as inte
import matplotlib.pyplot as plt
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter22.08-Summary-and-Problems.html
# https://ifcuriousthenlearn.com/blog/2015/06/09/mechanical-vibrations-with-python/

# Model the calling function [X_dot]:
# ------------------------------------------------------------------
def MassSpringDamper(t,state):
    '''
    k=spring constant, Newtons per metre
    m=mass, Kilograms
    c=damping coefficient, Newton*second / meter

    for a mass,spring
        xdd = ((-k*x)/m) + g
    for a mass, spring, damper
        xdd = -k*x/m -c*xd-g
    for a mass, spring, damper with forcing function
        xdd = -k*x/m -c*xd-g + cos(4*t-pi/4)
    '''

    k = 1e6  # spring constant, kN/m
    m = 100  # mass, Kg
    c = 30  # damping coefficient
    # unpack the state vector
    x,y,z, dx,dy,dz = state  # displacement,x and velocity x'
    # print(type(x))
    # print(type(dx))
    g = [0,0,-9.8]  # metres per second**2
    # compute acceleration xdd = x''
    def drag(v,u):
        cd= 1.0
        A=2.0
        row=1025.0
        # print('u is '+str(u))
        # print('v is ' + str(v))
        ur=np.array(u)-np.array(v)
        # print('Ur is '+str(ur))
        return 0.5*row*cd*A*ur*np.linalg.norm(ur)
    # exit()

    # print(drag([dx, dy, dz], [1, 0, 0]))
    ddx = -k * x / m - c * dx + g[0] + drag([dx.item(0), dy.item(0), dz.item(0)], [1, 0, 0])[0] / m
    ddy = -k * y / m - c * dy + g[1] + drag([dx.item(0), dy.item(0), dz.item(0)], [1, 0, 0])[1] / m
    ddz = -k * z / m - c * dz + g[2] + drag([dx.item(0), dy.item(0), dz.item(0)], [1, 0, 0])[2] / m

    return [dx,dy,dz, ddx,ddy,ddz]

def get_data(odesolver):
    # collect data
    t_values = []
    x_values = []
    dx_values = []
    for i in range(500):
        # get solution step state
        odesolver.step()
        t_values.append(odesolver.t)
        x_values.append(odesolver.y[0:3])
        dx_values.append(odesolver.y[3:])
        # break loop after modeling is finished
        if odesolver.status == 'finished':
            break
    return t_values, x_values,dx_values
# Define the initial conditions and integration boundaries:


# ------------------------------------------------------------------
time_step = 0.01
t_upper = 10.5
t_lower = 0

initial_conditions = [0, 0, 0, 0,0,0]  # [displacement(t=0) = 0, velocity(t=0) = 0]

solution1 = inte.RK45(fun=MassSpringDamper, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)
# solution2 = inte.Radau(fun=MassSpringDamper, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)
# solution3 = inte.LSODA(fun=MassSpringDamper, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)

# plt.figure()

t,x,dx=get_data(solution1)
Nodes_position = {}
cwd=os.getcwd()
for i in range(len(t)):
    Nodes_position[round(t[i],3)]=x[i].tolist()

with open(os.path.join(cwd,'Nodes_position.json'), "w") as json_file:
   json.dump(Nodes_position,json_file)
json_file.close()
# plt.show()

