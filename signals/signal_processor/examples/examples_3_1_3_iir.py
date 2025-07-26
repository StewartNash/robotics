from examples_3_1_1_iir import compute_filter_order, lowpass_computations, FilterType, FilterFamily
from examples_3_1_1_iir import butterworth_analog_poles, butterworth_digital_poles, frequency_scaling_parameter, iir_filter
from utility import polynomial_coefficients
import numpy as np

# Example 4.A.2
# Chebyshev lowpass (order = 6)
# Lowpass filter format

F = 2500	# samples/second
fp1 = 0.00	# Hz
fp2 = 500.00	# Hz
fs1 = 0.00	# Hz
fs2 = 750	# Hz
ripple_passband = 0.1737	# dB
ripple_stopband = 40.00		# dB

maximum_passband_attenuation = ripple_passband
minimum_stopband_attenuation = ripple_stopband
passband_edge_frequency = fp2
stopband_edge_frequency = fs2
sampling_frequency = F

K, A = lowpass_computations(maximum_passband_attenuation,
	minimum_stopband_attenuation,
	passband_edge_frequency,
	stopband_edge_frequency,
	sampling_frequency)

parameter_K = K
parameter_A = A
filter_family = FilterFamily.BUTTERWORTH

N = compute_filter_order(parameter_K, parameter_A, filter_family)
analog_poles = butterworth_analog_poles(10)
Ap = maximum_passband_attenuation
alpha = frequency_scaling_parameter(FilterFamily.BUTTERWORTH, fp2, F, N=N, Ap=Ap)
digital_poles = butterworth_digital_poles(analog_poles, alpha)
print(digital_poles)
denominator_coefficients = polynomial_coefficients(digital_poles)
numerator_coefficients = [(2, 1)] * len(denominator_coefficients)
print(denominator_coefficients)

b = numerator_coefficients
a = denominator_coefficients
fs = 1000  # Hz
t = np.arange(0, 0.1, 1 / F)
x = np.sin(2 * np.pi * 125 * t) + 0.5 * np.random.randn(len(t))  # 125 Hz sine + noise

t = t.tolist()
x = x.tolist()

# Filter the signal
y = iir_filter(x, b, a)

# Plot result
import matplotlib.pyplot as plt
plt.figure()
plt.plot(t, x, label="Input")
plt.plot(t, y, label="Filtered")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
