import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Van der Pol parameters
mu = 2.0

# System of equations
def van_der_pol(t, z):
    x, y = z
    dxdt = y
    dydt = mu * (1 - x**2) * y - x
    return [dxdt, dydt]

# Time span
t_span = (0, 40)
t_eval = np.linspace(*t_span, 10000)

# Initial condition
z0 = [2, 0]

# Solve the system
sol = solve_ivp(van_der_pol, t_span, z0, t_eval=t_eval, method='RK45')

# Extract x and y
x = sol.y[0]
y = sol.y[1]

# Phase plot
plt.figure(figsize=(8, 6))
plt.plot(x, y)
plt.title('Phase Plot of Van der Pol Oscillator (μ = {})'.format(mu))
plt.xlabel('x')
plt.ylabel('dx/dt')
plt.grid(True)
plt.axis('equal')
plt.show()

