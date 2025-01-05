import numpy as np
import matplotlib.pyplot as plt

# data = np.loadtxt("parallel_burgers_solution.txt")
data = np.loadtxt("burgers_solution_norm.txt")
x = data[:, 0]
y = data[:, 1]
U = data[:, 2]
V = data[:, 3]

Nx = len(np.unique(x))
Ny = len(np.unique(y))
x = x.reshape(Ny, Nx)
y = y.reshape(Ny, Nx)
U = U.reshape(Ny, Nx)
V = V.reshape(Ny, Nx)

velocity_norm = np.sqrt(U**2 + V**2)

fig, ax = plt.subplots(1, 1, figsize=(6, 6))

c = ax.contourf(x, y, velocity_norm, levels=50, cmap="viridis")
ax.set_xlabel("x", fontsize=12)
ax.set_ylabel("y", fontsize=12)
fig.colorbar(c, ax=ax, label="Velocity Norm")

plt.tight_layout()
# plt.savefig('parallel plot.png')
plt.show()