import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp

# Van der Pol parameters
mu = 2.0

# Van der Pol system
def van_der_pol(t, z):
    x, y = z
    dxdt = y
    dydt = mu * (1 - x**2) * y - x
    return [dxdt, dydt]

# Solve ODE
t_span = (0, 40)
t_eval = np.linspace(*t_span, 2000)
z0 = [2, 0]
sol = solve_ivp(van_der_pol, t_span, z0, t_eval=t_eval)

x = sol.y[0]
y = sol.y[1]

# Create plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-3, 3)
ax.set_ylim(-4, 4)
ax.set_xlabel('x')
ax.set_ylabel('dx/dt')
ax.set_title('Van der Pol Phase Animation (μ = {})'.format(mu))
ax.grid(True)
ax.axis('equal')

# Plot background vector field
X, Y = np.meshgrid(np.linspace(-3, 3, 20), np.linspace(-4, 4, 20))
U = Y
V = mu * (1 - X**2) * Y - X
ax.streamplot(X, Y, U, V, color='gray', density=0.6, linewidth=0.5)

# Initialize animated line
line, = ax.plot([], [], 'r-', lw=2)
point, = ax.plot([], [], 'ro')  # moving point

# Animation function
def update(frame):
    line.set_data(x[:frame], y[:frame])
    point.set_data(x[frame], y[frame])
    return line, point

# Animate
ani = FuncAnimation(fig, update, frames=len(x), interval=20, blit=True)

plt.show()
