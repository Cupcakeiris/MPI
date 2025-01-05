import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('execution_times.txt')
size = data[:, 0]
time = data[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(size, time, marker='o', linestyle='-', color='b')
plt.title("Execution times of MPI Program")
plt.xlabel("Num of procs")
plt.ylabel("Execution time / s")
plt.grid(True)
plt.show()

plt.savefig('execution_times.png')