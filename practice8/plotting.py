import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

# Numpy is so op
data_analytical = np.loadtxt('analytical_solution.txt')
time_analytical = data_analytical[:, 0]
position_analytical = data_analytical[:, 1]
u_values_analytical = data_analytical[:, 2]

data_lax = np.loadtxt('numerical_solution_lax.txt')
time_lax = data_lax[:, 0]
position_lax = data_lax[:, 1]
u_values_lax = data_lax[:, 2]

data_back = np.loadtxt('numerical_solution_back.txt')
time_back = data_lax[:, 0]
position_back = data_lax[:, 1]
u_values_back = data_lax[:, 2]

def plot(t0, filename):
    u_analytical = u_values_analytical[time_analytical == t0]
    pos_analytical = position_analytical[time_analytical == t0]

    u_lax = u_values_lax[time_lax == t0]
    pos_lax = position_lax[time_lax == t0]

    u_back = u_values_back[time_lax == t0]
    pos_back = position_back[time_lax == t0]

    plt.figure(figsize=(10, 6))
    plt.plot(pos_analytical, u_analytical, label='Analytical Solution', linestyle='-', color='red')
    plt.plot(pos_lax, u_lax, label='Lax-Wendroff Approximation (parallized)', linestyle='-.', color='blue')
    plt.plot(pos_back, u_back, label='Backwards Approximation (parallized)', linestyle=':', color='green', linewidth=3)

    plt.title(f'Solution of the Transport Equation at t={t0:.2f}')
    plt.xlabel('Position (x)')
    plt.ylabel('u(x,t)')
    plt.grid()
    plt.legend()


    plt.savefig(filename)
    plt.show()

def rmse_for_method(t0, numerical, t_num, analytical, t_an):
    a = analytical[t_an == t0]
    b = numerical[t_num == t0]
    return np.sqrt(mean_squared_error(a, b))

plot(0.9, 'plot_t0.9.png')
plot(0.5, 'plot_t0.5.png')
plot(0.1, 'plot_t0.1.png')

# Calculating RMSE
rmse_lax = rmse_for_method(0.9, u_values_lax, time_lax, u_values_analytical, time_analytical)
rmse_back = rmse_for_method(0.9, u_values_back, time_back, u_values_analytical, time_analytical)
print(f"RMSE for Lax-Wendroff at t=0.9: {rmse_lax:.4f}")
print(f"RMSE for Backward Approximation at t=0.9: {rmse_back:.4f}")

rmse_lax = rmse_for_method(0.5, u_values_lax, time_lax, u_values_analytical, time_analytical)
rmse_back = rmse_for_method(0.5, u_values_back, time_back, u_values_analytical, time_analytical)
print(f"RMSE for Lax-Wendroff at t=0.5: {rmse_lax:.4f}")
print(f"RMSE for Backward Approximation at t=0.5: {rmse_back:.4f}")

rmse_lax = rmse_for_method(0.1, u_values_lax, time_lax, u_values_analytical, time_analytical)
rmse_back = rmse_for_method(0.1, u_values_back, time_back, u_values_analytical, time_analytical)
print(f"RMSE for Lax-Wendroff at t=0.1: {rmse_lax:.4f}")
print(f"RMSE for Backward Approximation at t=0.1: {rmse_back:.4f}")
