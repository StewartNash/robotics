import math
import numpy as np
import matplotlib.pyplot as plt
from examples_3_2_2_kaiser import kaiser_filter_order, kaiser_coefficients, kaiser_bandpass
from examples_3_2_2_kaiser import FilterType

def compute_ripple_parameters(minimum_stopband_attenuation, specified_passband_ripple):
	delta_stopband = 10 ** (-0.05 * minimum_stopband_attenuation)
	delta_passband = (10 ** (0.05 * specified_passband_ripple) - 1) / (10 ** (0.05 * specified_passband_ripple) + 1)
	delta = min(delta_passband, delta_stopband)
		
	return (delta, delta_stopband, delta_passband)
	
def compute_actual_passband_ripple(delta_passband):
	actual_passband_ripple = 20 * math.log10((1 + delta_passband) / (1 - delta_passband))
	
	return actual_passband_ripple
	
def transfer_function(z, impulse_response):
	return []
	
def magnitude_response(omega, impulse_response, sampling_frequency, filter_order):
	nk = filter_order
	frequency = omega / (2 * math.pi)
	period = 1 / sampling_frequency
	response_magnitude = impulse_response[0]
	for i in range(1, int((nk - 1) / 2) + 1):
		response_magnitude += impulse_response[i] * math.cos(2 * math.pi * frequency * i * period)
		
	return response_magnitude

#FREQUENCY_RESOLUTION = 100
FREQUENCY_RESOLUTION = 500

# Example 5.6
# -----------
# Desing a bandpass filter satisfying the following specifications
# actual_passband_ripple	(leq)	0.5	[dB]
# minimum_stopband_attenuation	(geq)	35	[dB]
# passband_frequency_low	(eq)	120	[Hz]
# passband_frequency_high	(eq)	180	[Hz]
# stopband_frequency_low	(eq)	60	[Hz]
# stopband_frequency_high	(eq)	240	[Hz]
# sampling_frequency		(eq)	600	[Hz]

actual_passband_ripple = 0.5
minimum_stopband_attenuation = 35
passband_frequency_low = 120
passband_frequency_high = 180
stopband_frequency_low = 60
stopband_frequency_high = 240
sampling_frequency = 600

transition_bandwidth_low = abs(passband_frequency_low - stopband_frequency_low)
transition_bandwidth_high = abs(passband_frequency_high - stopband_frequency_high)
transition_bandwidth = min(transition_bandwidth_low, transition_bandwidth_high)

(filter_order, delta, minimum_stopband_attenuation, parameter_d) = kaiser_filter_order(filter_type=FilterType.BANDPASS,
	passband_frequency_low=passband_frequency_low,
	passband_frequency_high=passband_frequency_high,
	stopband_frequency_low=stopband_frequency_low,
	stopband_frequency_high=stopband_frequency_high,
	sampling_frequency=sampling_frequency,
	specified_passband_ripple=actual_passband_ripple,
	minimum_stopband_attenuation=minimum_stopband_attenuation)

print("Filter order:", filter_order)

(kaiser_coeffs, alpha) = kaiser_coefficients(filter_order=filter_order,
	actual_passband_ripple=actual_passband_ripple)

impulse_response = kaiser_bandpass(passband_frequency_low=passband_frequency_low,
	passband_frequency_high = passband_frequency_high,
	sampling_frequency=sampling_frequency,
	transition_bandwidth=transition_bandwidth,
	filter_order=filter_order,
	kaiser_coeffs=kaiser_coeffs)	

maximum_frequency = sampling_frequency / 2
frequencies = np.linspace(0, maximum_frequency, FREQUENCY_RESOLUTION)
angular_frequencies = frequencies * 2 * math.pi
response_magnitude = [magnitude_response(w, impulse_response, sampling_frequency, filter_order) for w in angular_frequencies]
response_magnitude = [20 * math.log10(max(abs(m), 1e-10)) for m in response_magnitude]

#print(frequencies)
#print(response_magnitude)

plt.plot(frequencies, response_magnitude)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Response (dB)")
plt.title("Bandpass Filter Magnitude Response")
plt.show()

