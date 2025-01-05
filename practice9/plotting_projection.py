import numpy as np
import matplotlib.pyplot as plt

def plot_transport_solution(data_filename, t0, filename_x, filename_y):
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

    y_fixed_index = len(position_y) // 2
    x_fixed_index = len(position_x) // 2

    u_along_x = u_at_t0[y_fixed_index, :]
    u_along_y = u_at_t0[:, x_fixed_index]

    # Plotting u along the x-axis
    plt.figure(figsize=(10, 6))
    plt.plot(position_x, u_along_x, label=f'u(x, y={position_y[y_fixed_index]:.2f}, t={t0})', color='blue')
    plt.title(f'Change of u along x-axis at t={t0}')
    plt.xlabel('Position X')
    plt.ylabel('u(x, y, t)')
    plt.grid()
    plt.legend()
    
    plt.savefig(filename_x)
    plt.show()

    # Plotting u along the y-axis
    plt.figure(figsize=(10, 6))
    plt.plot(position_y, u_along_y, label=f'u(x={position_x[x_fixed_index]:.2f}, y, t={t0})', color='green')
    plt.title(f'Change of u along y-axis at t={t0}')
    plt.xlabel('Position Y')
    plt.ylabel('u(x, y, t)')
    plt.grid()
    plt.legend()
    
    plt.savefig(filename_y)
    plt.show()

# plot_transport_solution('transport_solution_2D.txt', 0.9, 'plot_x_t0.9.png', 'plot_y_t0.9.png')
# plot_transport_solution('transport_solution_2D.txt', 0.5, 'plot_x_t0.5.png', 'plot_y_t0.5.png')
# plot_transport_solution('transport_solution_2D.txt', 0.1, 'plot_x_t0.1.png', 'plot_y_t0.1.png')

# plot_transport_solution('transport_solution_2D_parallel.txt', 0.9, 'plot_x_t0.9_parallel.png', 'plot_y_t0.9_parallel.png')
# plot_transport_solution('transport_solution_2D_parallel.txt', 0.5, 'plot_x_t0.5_parallel.png', 'plot_y_t0.5_parallel.png')
# plot_transport_solution('transport_solution_2D_parallel.txt', 0.1, 'plot_x_t0.1_parallel.png', 'plot_y_t0.1_parallel.png')

# plot_transport_solution('transport_analytical_2D.txt', 0.9, 'plot_x_t0.9_analytical.png', 'plot_y_t0.9_analytical.png')
# plot_transport_solution('transport_analytical_2D.txt', 0.5, 'plot_x_t0.5_analytical.png', 'plot_y_t0.5_analytical.png')
# plot_transport_solution('transport_analytical_2D.txt', 0.1, 'plot_x_t0.1_analytical.png', 'plot_y_t0.1_analytical.png')