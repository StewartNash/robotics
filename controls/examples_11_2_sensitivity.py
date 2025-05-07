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

plt.figure()
control.nyquist_plot(GH, label='Uncompensated', color='blue')
control.nyquist_plot(compensated, label='Compensated', color='orange')
plt.axhline(0, color='gray', linestyle='--')
plt.axvline(0, color='gray', linestyle='--')
plt.title('Nyquist Plot: Uncompensated vs. Compensated')
plt.xlabel('Re')
plt.ylabel('Im')
plt.grid(True)
plt.axis('equal')
plt.legend()
plt.show()
