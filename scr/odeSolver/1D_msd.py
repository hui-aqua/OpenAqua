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
    x, xd = state  # displacement,x and velocity x'
    g = 9.8  # metres per second**2
    # compute acceleration xdd = x''
    def drag(xd,u):
        cd= 1.0
        A=2.0
        row=1025.0
        return 0.5*row*cd*A*pow((u-xd),2)*np.sign(u-xd)


    xdd = -k * x / m - c * xd - g + drag(xd,-1.0)/m
    return [xd, xdd]

def get_data(odesolver):
    # collect data
    t_values = []
    x_values = []
    dx_values = []
    for i in range(100):
        # get solution step state
        odesolver.step()
        t_values.append(odesolver.t)
        x_values.append(odesolver.y[0])
        dx_values.append(odesolver.y[1])
        # break loop after modeling is finished
        if odesolver.status == 'finished':
            break
    return t_values, x_values,dx_values
# Define the initial conditions and integration boundaries:


# ------------------------------------------------------------------
time_step = 0.01
t_upper = 1.5
t_lower = 0

initial_conditions = np.array([0, 0])  # [displacement(t=0) = 0, velocity(t=0) = 0]

solution1 = inte.RK45(fun=MassSpringDamper, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)
solution2 = inte.Radau(fun=MassSpringDamper, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)
solution3 = inte.LSODA(fun=MassSpringDamper, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)

plt.figure()

t,x,dx=get_data(solution1)
plt.plot(t,x)
t,x,dx=get_data(solution2)
plt.plot(t,x)
t,x,dx=get_data(solution3)
plt.plot(t,x)

plt.show()

