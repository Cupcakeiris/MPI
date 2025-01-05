import numpy as np
import matplotlib.pyplot as plt
import glob

data_arrays = []

# Use glob to find all text files matching the pattern
for filename in glob.glob("laplace_*.txt"):
    data = np.loadtxt(filename)
    data_arrays.append(data)

# Concatenate all arrays into a single array
combined_data = np.concatenate(data_arrays, axis=0) 

x, y, u = combined_data[:, 0], combined_data[:, 1], combined_data[:, 2]

Nx = len(np.unique(x))
Ny = len(np.unique(y))

X = x.reshape((Nx, Ny))
Y = y.reshape((Nx, Ny))
U = u.reshape((Nx, Ny))

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, U, cmap="viridis")

fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x, y)")
plt.savefig("laplace_solution.png", dpi=300)
plt.show()