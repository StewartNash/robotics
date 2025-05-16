import numpy as np
import matplotlib.pyplot as plt
import control

m1 = 40
m2 = 20
k1 = 400
k2 = 200
c1 = 20
c2 = 10

A = np.array([
	[0, 1, 0, 0],
	[-(k1 + k2)/m1, -(c1 + c2)/m1, k2/m1, c2/m1],
	[0, 0, 0, 1],
	[k2/m2, c2/m2, -k2/m2, -c2/m2]
])

B = np.array([
	[0, 0],
	[1/m1, 0],
	[0, 0],
	[0, 1/m2]
])

C = np.array([
	[1, 0, 0, 0],
	[0, 0, 1, 0]
])

D = np.array([
	[0, 0],
	[0, 0]
])

sys = control.StateSpace(A, B, C, D)
T = np.linspace(0, 10, 1000)
U = np.zeros((2, len(T)))
X0 = [0.2, 0, 0.2, 0]
#T, y, x = control.input_output_response(sys, T=T, U=U, X0=X0)
T, y = control.forced_response(sys, T=T, U=U, X0=X0)

plt.figure(figsize=(10, 5))
plt.plot(T, y[0], label="m1 position")
plt.plot(T, y[1], label="m2 position")
plt.title("Zero-input response of example 10-1")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("../manual/images/example_10_1.png", dpi=300, bbox_inches='tight')
plt.show()
