import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np

with open("solution_backward.txt", "r") as f:
    solution_backward = [float(line.strip()) for line in f]
with open("solution_central.txt", "r") as f:
    solution_central = [float(line.strip()) for line in f]
with open("analytical.txt", "r") as f:
    analytical = [float(line.strip()) for line in f]

Nx = len(solution_backward)
x = [i / (Nx - 1) for i in range(Nx)]

plt.plot(x, analytical, label="Analytical")
plt.plot(x, solution_backward, label="Solution backward", linestyle = '-.')
# plt.plot(x, solution_central, label="Solution central", linestyle=':')
plt.xlabel("x")
plt.ylabel("P(x)")
plt.title("Analytical and numerical solutions")
plt.grid(True)
plt.legend()
plt.show()