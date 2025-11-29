# VBD-MODE
Interview Exercises

This repository contains my solution to the VBD-MODE interview exercise.  
The task is to numerically solve a 3-dimensional ODE system using the 4th-order Runge–Kutta (RK4) method and then visualize the results.

The implementation is written in Python and follows the specification as requested.

---

## 1. Problem Description

We consider the system

\[
\begin{aligned}
\frac{dE}{dt} &= \beta A - \delta_E E - \mu_E E \\
\frac{dJ}{dt} &= \delta_E E - \delta_J J - \alpha J^2 - \mu_J J \\
\frac{dA}{dt} &= \omega \, \delta_J J - \mu_A A
\end{aligned}
\]

with state vector

\[
y(t) = [E(t), J(t), A(t)]^\top.
\]

### Parameters

The parameters are fixed as:

- \(\beta = 24\)
- \(\delta_E = 0.6\)
- \(\mu_E = 0.15\)
- \(\delta_J = 0.08\)
- \(\mu_J = 0.05\)
- \(\alpha = 0.003\)
- \(\omega = 0.5\)
- \(\mu_A = 0.1\)

### Initial Condition and Time Grid

- Initial condition:  
  \[
  y_0 = [E(0), J(0), A(0)] = [10,\; 0,\; 0]
  \]

- Time interval and step size:
  - \(t_0 = 0\)
  - \(t_{\text{end}} = 365\)
  - \(h = 0.01\)

The time vector **T** is constructed as

\[
T_i = t_0 + h \cdot (i - 1), \quad i = 1, 2, \dots, N
\]

so that the first element is \(t_0\) and the last element is \(t_{\text{end}}\).

---

## 2. Numerical Method

### Algorithm 2 – RHS Function \(f\)

The function \(f(t_n, y_n)\) returns the time derivative:

```text
E = y[0]
J = y[1]
A = y[2]

dE = β * A - δ_E * E - μ_E * E
dJ = δ_E * E - δ_J * J - α * J^2 - μ_J * J
dA = ω * δ_J * J - μ_A * A

dy = [dE, dJ, dA]

