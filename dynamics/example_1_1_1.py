import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the Van der Pol oscillator system
def vanderpol_ode(t, z, mu):
    x, y = z
    dxdt = y
    dydt = mu * (1 - x**2) * y - x
    return [dxdt, dydt]

# Parameters
mu = 1.0  # Damping parameter
initial_conditions = [0.1, 0.0]  # [x0, y0]
t_span = (0, 30)  # Time interval for simulation
t_eval = np.linspace(t_span[0], t_span[1], 500) # Points at which to store the solution

# Solve the ODE
sol = solve_ivp(vanderpol_ode, t_span, initial_conditions, 
                args=(mu,), t_eval=t_eval, method='RK45')

# Plotting the results
plt.figure(figsize=(10, 5))

# Position vs. Time
plt.subplot(1, 2, 1)
plt.plot(sol.t, sol.y[0])
plt.xlabel("Time (t)")
plt.ylabel("Position (x)")
plt.title("Van der Pol Oscillator: Position vs. Time")
plt.grid(True)

# Phase Portrait (Velocity vs. Position)
plt.subplot(1, 2, 2)
plt.plot(sol.y[0], sol.y[1])
plt.xlabel("Position (x)")
plt.ylabel("Velocity (y)")
plt.title("Van der Pol Oscillator: Phase Portrait")
plt.grid(True)

plt.tight_layout()
plt.show()
