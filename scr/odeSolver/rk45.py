import numpy as np
import scipy.integrate as inte
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/63953924/issue-understanding-scipy-integrate-rk45-requirements
# https://stackoverflow.com/questions/58999268/rk45-ode-solver-python3

# Model the calling function [X_dot]:
# ------------------------------------------------------------------
def model(t, X):
    # Define the mass, stiffness and damping [m, c, k]:
    m = 2
    c = 10
    k = 1500

    # Define the piecewise forcing function [F]:
    if (t >= 0 and t < 0.1):
        F = 200 * t
    if (t >= 0.1 and t < 0.25):
        F = 20
    else:
        F = 0

    E_matrix = np.array([[0,         1],
                         [(-k / m), (-c / m)]
                         ])
    Q_matrix = np.array([[0],[ F / m]])

    return np.matmul(E_matrix,X) + Q_matrix

def get_data(odesolver):
    # collect data
    t_values = []
    y_values = []
    for i in range(1000):
        # get solution step state
        odesolver.step()
        t_values.append(odesolver.t)
        y_values.append(odesolver.y[0])
        # break loop after modeling is finished
        if odesolver.status == 'finished':
            break
    return t_values, y_values
# Define the initial conditions and integration boundaries:
# ------------------------------------------------------------------
time_step = 0.01
t_upper = 10.5
t_lower = 0

initial_conditions = np.array([0, 0])  # [displacement(t=0) = 0, velocity(t=0) = 0]

solution = inte.RK45(fun=model, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)
solution2 = inte.Radau(fun=model, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)
solution3 = inte.LSODA(fun=model, t0=t_lower, y0=initial_conditions, t_bound=t_upper, vectorized=True)

plt.figure()
x,y=get_data(solution)
plt.plot(x,y)

x,y=get_data(solution2)
plt.plot(x,y)

x,y=get_data(solution3)
plt.plot(x,y)


plt.show()

