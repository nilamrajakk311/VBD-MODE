import numpy as np
import matplotlib.pyplot as plt

# Parameters 
beta = 24
delta_E = 0.6
mu_E = 0.15
delta_J = 0.08
mu_J = 0.05
alpha = 0.003
omega = 0.5
mu_A = 0.1

# Time settings
t0 = 0
t_end = 365      
h = 0.01

# Initial condition y0 = [E0, J0, A0]
y0 = np.array([10.0, 0.0, 0.0])

# Number of steps: integers 1 to t_end/h
n_steps = int(round((t_end - t0) / h))

# Vector T with pattern: t0 + h * (i - 1), i = 1..n_steps+1
T = t0 + h * np.arange(0, n_steps + 1)

def f(t, y):  # Algorithm 2
    E = y[0]
    J = y[1]
    A = y[2]
    
    dE = beta * A - delta_E * E - mu_E * E
    dJ = delta_E * E - delta_J * J - alpha * J**2 - mu_J * J
    dA = omega * delta_J * J - mu_A * A
    
    return np.array([dE, dJ, dA])

def RK4(f, y_n, t_n, h):  # Algorithm 1 
    k1 = f(t_n, y_n)
    k2 = f(t_n + 0.5*h, y_n + 0.5*h*k1)
    k3 = f(t_n + 0.5*h, y_n + 0.5*h*k2)
    k4 = f(t_n + h,     y_n + h*k3)

    y_next = y_n + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    return y_next

# Save y0 together with RK4 outputs in Y
Y = [y0]
y_current = y0

# Loop over integers 1 to t_end/h (i = 1..n_steps)
for i in range(1, n_steps + 1):
    t_current = T[i - 1]       # first iteration uses first element of T
    y_next = RK4(f, y_current, t_current, h)
    Y.append(y_next)
    y_current = y_next

Y = np.array(Y)

# Plotting
plt.figure(figsize=(12, 8))

# Subplot 1 — E(t)
plt.subplot(3, 1, 1)
plt.plot(T, Y[:, 0])
plt.xlabel("T")
plt.ylabel("E")

# Subplot 2 — J(t)
plt.subplot(3, 1, 2)
plt.plot(T, Y[:, 1])
plt.xlabel("T")
plt.ylabel("J")

# Subplot 3 — A(t)
plt.subplot(3, 1, 3)
plt.plot(T, Y[:, 2])
plt.xlabel("T")
plt.ylabel("A")

plt.tight_layout()
# Save the figure
plt.savefig("VBD_simulation_plot.png", dpi=300)
plt.show()
