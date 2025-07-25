from examples_3_1_1_iir import compute_filter_order, lowpass_computations, FilterType, FilterFamily
from examples_3_1_1_iir import butterworth_analog_poles

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

print(butterworth_analog_poles(10))
