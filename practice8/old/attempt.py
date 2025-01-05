import numpy as np
import matplotlib.pyplot as plt


L = 1   
T = 1  
Nx = 100
Nt = 1000
dx = L / (Nx - 1)
dt = T / Nt

x = np.linspace(0, L, Nx)
t = np.linspace(0, T, Nt)

u = np.cos(np.pi * x)

def boundary_condition(t):
    return np.exp(-t)

def analytical_solution(x, t):
    return np.cos(np.pi * (x - 2 * t))

u_central = np.zeros((Nt, Nx))
u_backward = np.zeros((Nt, Nx))
u_central[0, :] = u
u_backward[0, :] = u

# Central Difference Method
for n in range(0, Nt - 1):
    for i in range(1, Nx - 1):
        u_central[n + 1, i] = u_central[n, i] - dt / (2 * dx) * (u_central[n, i + 1] - u_central[n, i - 1])
    u_central[n + 1, 0] = boundary_condition(t[n + 1])
    u_central[n + 1, -1] = u_central[n + 1, -2]

# Backward Difference Method
for n in range(0, Nt - 1):
    for i in range(1, Nx):
        u_backward[n + 1, i] = u_backward[n, i] - 2 * dt / dx * (u_backward[n, i] - u_backward[n, i - 1])
    u_backward[n + 1, 0] = boundary_condition(t[n + 1])
    u_backward[n + 1, -1] = u_backward[n + 1, -2]

t0 = 0.99
n0 = int(t0 / dt) # To get index of time

plt.figure(figsize=(10, 5))
plt.plot(x, u_central[n0, :], label='Central Difference Method', linestyle=':')
plt.plot(x, u_backward[n0, :], label='Backward Difference Method', linestyle='--')
plt.plot(x, analytical_solution(x, t0), label='Analytical Solution')
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title(f'Numerical Solution of Transport Equation at t = {t0}')
plt.legend()
plt.show()
