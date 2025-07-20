import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

def van_der_pol(t, z, mu):
    x, y = z
    dxdt = y
    dydt = mu * (1 - x**2) * y - x
    return [dxdt, dydt]

# Parameter range for μ
mu_values = np.linspace(0.1, 10.0, 200)
maxima_by_mu = []

# Loop over μ
for mu in mu_values:
    # Solve for each μ
    t_span = (0, 100)
    t_eval = np.linspace(*t_span, 5000)
    z0 = [2.0, 0.0]
    sol = solve_ivp(lambda t, z: van_der_pol(t, z, mu),
                    t_span, z0, t_eval=t_eval, method='RK45')
    
    x = sol.y[0]
    
    # Discard transients (e.g., first half)
    x_steady = x[len(x)//2:]
    
    # Find peaks (local maxima)
    peaks, _ = find_peaks(x_steady, height=0)
    peak_values = x_steady[peaks]
    
    # Store each maximum for this μ
    for val in peak_values:
        maxima_by_mu.append((mu, val))

# Convert to arrays for plotting
mus, max_vals = zip(*maxima_by_mu)

# Plot bifurcation diagram
plt.figure(figsize=(10, 6))
plt.plot(mus, max_vals, 'k.', markersize=1)
plt.xlabel('μ')
plt.ylabel('x maxima (after transient)')
plt.title('Bifurcation Diagram of Van der Pol Oscillator')
plt.grid(True)
plt.show()
