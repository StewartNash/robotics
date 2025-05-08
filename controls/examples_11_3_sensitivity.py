import numpy
import matplotlib.pyplot as plt
import control

K1 = 1
p1 = 0.5
p2 = 1

numerator = K1
denominator = numpy.poly([0, -p1, -p2])
GH = control.TransferFunction(numerator, denominator)

a = 1 # numerator = [1, a]
b = 2 # denominator = [1, b]
P = control.TransferFunction([1, a], [1, b]) # P_lead

compensated = P * GH

omega = numpy.logspace(-2, 2, 500)
GH_response = control.frequency_response(GH, omega)
compensated_response = control.frequency_response(compensated, omega)

H = GH_response.fresp[0, 0, :]
Hc =  compensated_response.fresp[0, 0, :]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(H.real.flatten(), H.imag.flatten(), label='Uncompensated', color='blue')
ax1.plot(Hc.real.flatten(), Hc.imag.flatten(), label='Compensated', color='orange')
theta = numpy.linspace(0, 2 * numpy.pi, 500)
ax1.plot(numpy.cos(theta), numpy.sin(theta), '--', color='gray', label='Unit Circle')
ax1.axhline(0, color='gray', linestyle='--')
ax1.axvline(0, color='gray', linestyle='--')
ax1.set_xlim(-10, 1)
ax1.set_ylim(-10, 1)
ax1.set_title("Full Nyquist Plot")
ax1.set_xlabel("Re")
ax1.set_ylabel("Im")
ax1.legend()
ax1.grid(True)
#ax1.axis("equal")

ax2.plot(H.real.flatten(), H.imag.flatten(), label="Uncompensated", color='blue')
ax2.plot(Hc.real.flatten(), Hc.imag.flatten(), label='Compensated', color='orange')
ax2.plot(numpy.cos(theta), numpy.sin(theta), '--', color='gray', label='Unit Circle')
ax2.axhline(0, color='gray', linestyle='--')
ax2.axvline(0, color='gray', linestyle='--')
ax2.set_xlim(-2, 2)
ax2.set_ylim(-2, 2)
ax2.set_title('Nyquist (detail)')
ax2.set_xlabel('Re')
ax2.set_ylabel('Im')
ax2.legend()
ax2.grid(True)
#ax2.axis('equal')

plt.tight_layout()
plt.show()
