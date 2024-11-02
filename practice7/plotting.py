import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np

# Read data
with open("analytical_x.txt", "r") as f:
    analytical_x = [float(line.strip()) for line in f]
with open("newP.txt", "r") as f:
    newP = [float(line.strip()) for line in f]
with open("newP_parallelized.txt", "r") as f:
    newP_parallelized = [float(line.strip()) for line in f]

Nx = len(analytical_x)
x = [i / (Nx - 1) for i in range(Nx)]

# Plot the analytical solution
plt.plot(x, analytical_x, label="Analytical solution")
plt.plot(x, newP, label="Explicit method", linestyle='--')
plt.plot(x, newP_parallelized, label="Parallelized", linestyle='-.')
plt.xlabel("x")
plt.ylabel("P(x)")
plt.title("Analytical and numerical solutions")
plt.grid(True)
plt.legend()
plt.show()

print("RMSE:",np.sqrt(np.mean(np.array(analytical_x) - np.array(newP))**2))
print("RMSE (Parallel):",np.sqrt(np.mean(np.array(analytical_x) - np.array(newP_parallelized))**2))