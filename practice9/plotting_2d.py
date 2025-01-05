import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

def plot_transport_solution(data_filename, t0, output_filename):
    data = np.loadtxt(data_filename)

    time_values = np.unique(data[:, 0])
    position_x = np.unique(data[:, 1])
    position_y = np.unique(data[:, 2])

    u_values = np.zeros((len(time_values), len(position_y), len(position_x)))

    for row in data:
        t, x, y, u = row
        t_index = np.where(time_values == t)[0][0]
        x_index = np.where(position_x == x)[0][0]
        y_index = np.where(position_y == y)[0][0]
        u_values[t_index, y_index, x_index] = u

    t_index = np.argmin(np.abs(time_values - t0))
    u_at_t0 = u_values[t_index]

    # print(u_at_t0)


    X, Y = np.meshgrid(position_x, position_y)

    # Plotting
    plt.figure(figsize=(10, 6))
    contour = plt.contourf(X, Y, u_at_t0, levels=50, cmap='viridis')
    plt.colorbar(contour)
    
    plt.title(f'Solution of the Transport Equation at t={t0:.2f}')
    plt.xlabel('Position X')
    plt.ylabel('Position Y')
    
    plt.savefig(output_filename)
    plt.show()

def calculate_rmse(numerical_filename, analytical_filename, t0):
    num_data = np.loadtxt(numerical_filename)
    ana_data = np.loadtxt(analytical_filename)

    num_data_t0 = num_data[np.isclose(num_data[:, 0], t0)]
    ana_data_t0 = ana_data[np.isclose(ana_data[:, 0], t0)]

    rmse = np.sqrt(mean_squared_error(num_data_t0[:, 3], ana_data_t0[:, 3]))
    print(f"RMSE at t={t0}: {rmse}")
    return rmse

# print(calculate_rmse('transport_solution_2D_parallel.txt', 'transport_analytical_2D.txt', 0.9))

# plot_transport_solution('transport_solution_2D.txt', 0.1, 'plot_t0.1.png')
# plot_transport_solution('transport_solution_2D.txt', 0.5, 'plot_t0.5.png')
# plot_transport_solution('transport_solution_2D.txt', 0.9, 'plot_t0.9.png')


plot_transport_solution('transport_solution_2D_parallel.txt', 0.1, 'plot_t0.1_parallel.png')
plot_transport_solution('transport_solution_2D_parallel.txt', 0.5, 'plot_t0.5_parallel.png')
plot_transport_solution('transport_solution_2D_parallel.txt', 0.9, 'plot_t0.9_parallel.png')

# plot_transport_solution('transport_analytical_2D.txt', 0.1, 'plot_t0.1_analytical.png')
# plot_transport_solution('transport_analytical_2D.txt', 0.5, 'plot_t0.5_analytical.png')
# plot_transport_solution('transport_analytical_2D.txt', 0.9, 'plot_t0.9_analytical.png')