import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("parallel_burgers_solution.txt")

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


fig, axs = plt.subplots(2, 1, figsize=(8, 8))

c1 = axs[0].contourf(x, y, U, levels=50, cmap="viridis")
axs[0].set_title("Velocity Component U", fontsize=16)
axs[0].set_xlabel("x", fontsize=12)
axs[0].set_ylabel("y", fontsize=12)
fig.colorbar(c1, ax=axs[0], label="U Velocity")


c2 = axs[1].contourf(x, y, V, levels=50, cmap="viridis")
axs[1].set_title("Velocity Component V", fontsize=16)
axs[1].set_xlabel("x", fontsize=12)
axs[1].set_ylabel("y", fontsize=12)
fig.colorbar(c2, ax=axs[1], label="V Velocity")

plt.tight_layout()
plt.show()