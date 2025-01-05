import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

data_lax = np.loadtxt('numerical_solution_lax.txt')
time_lax = data_lax[:, 0]
position_lax = data_lax[:, 1]
u_values_lax = data_lax[:, 2]

data_ana = np.loadtxt('transport_solution.txt')
time_ana = data_ana[:, 0]
position_ana = data_ana[:, 1]
u_values_ana = data_ana[:, 2]

def plot(t0, filename):

    u_lax = u_values_lax[time_lax == t0]
    pos_lax = position_lax[time_lax == t0]

    u_ana = u_values_ana[time_ana == t0]
    pos_ana = position_ana[time_ana == t0]

    plt.figure(figsize=(10, 6))
    plt.plot(pos_ana, u_ana, label='Analytical Solution', linestyle='-', color='red')
    plt.plot(pos_lax, u_lax, label='Lax-Friedrichs Approximation (parallized)', linestyle='--', color='blue')

    plt.title(f'Solution of the Transport Equation at t={t0:.2f}')
    plt.xlabel('Position (x)')
    plt.ylabel('u(x,t)')
    plt.grid()
    plt.legend()


    plt.savefig(filename)
    plt.show()

def rmse_for_method(t0, numerical, t_num, analytical):
    a = analytical[t_num == t0]
    b = numerical[t_num == t0]
    return np.sqrt(mean_squared_error(a, b))

plot(1.26, 'plot_t1.26.png')
plot(2.46, 'plot_t2.46.png')
plot(4.98, 'plot_t4.98.png')

t_ = 2.46
print("RMSE for time",t_)
print(rmse_for_method(t0 = t_, numerical = u_values_lax, t_num= time_lax,
                analytical = u_values_ana))