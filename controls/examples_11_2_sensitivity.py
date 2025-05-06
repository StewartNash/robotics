import numpy
import matplotlib.pyplot as plt
import control

K1 = 1
p1 = 0.5
p2 = 1

numerator = K1
denominator = numpy.poly([0, -p1, -p2])
GH = control.TransferFunction(numerator, denominator)

#plt.figure()
#control.nyquist(GH, omega_limits=(0.01, 100), omega_num=1000)
#plt.show()

omega = numpy.logspace(-2, 2, 100)
magnitude, phase, w = control.freqresp(GH, omega)

plt.figure()
ax = plt.subplot(111, projection='polar')
ax.plot(phase[0], magnitude[0])
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.grid(True)
plt.show()
