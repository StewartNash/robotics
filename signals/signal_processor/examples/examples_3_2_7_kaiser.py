import math
import numpy as np
import enum
import matplotlib.pyplot as plt


class FilterType(enum.Enum):
    LOWPASS = 1
    HIGHPASS = 2
    BANDPASS = 3
    BANDSTOP = 4

def kaiser_coefficients(filter_order, minimum_stopband_attenuation):
    n = filter_order
    if minimum_stopband_attenuation <= 21.0:
        alpha = 0
    elif minimum_stopband_attenuation <= 50.0:
        alpha = 0.5842 * (minimum_stopband_attenuation - 21) ** 0.4 + 0.07886 * (minimum_stopband_attenuation - 21)
    else:
        alpha = 0.1102 * (minimum_stopband_attenuation - 8.7)
    
    kaiser_coeffs = [0.0] * n
    for i in range(n):
    	beta = alpha * math.sqrt(1 - (2 * i / (n - 1) - 1) ** 2)
    	kaiser_coeffs[i] = zo_mod_bessel_fk(beta, 31) / zo_mod_bessel_fk(alpha, 31)
    print("Kaiser coefficients")
    print("-------------------")
    for i in range(len(kaiser_coeffs)):
    	print(kaiser_coeffs[i])

    return (kaiser_coeffs, alpha)

def zo_mod_bessel_fk(x, N=25):
	output = 0
	for i in range(N + 1):
		output += ((x / 2) ** i / math.factorial(i)) ** 2
		
	return output

def kaiser_filter_order(filter_type,
    passband_frequency_low,
    passband_frequency_high,
    stopband_frequency_low,
    stopband_frequency_high,
    sampling_frequency,
    specified_passband_ripple,
    minimum_stopband_attenuation):

    transition_bandwidth_low  = abs(passband_frequency_high - stopband_frequency_high)
    transition_bandwidth_high = abs(passband_frequency_low - stopband_frequency_low)

    if filter_type == FilterType.LOWPASS:
        transition_bandwidth = transition_bandwidth_low
    elif filter_type == FilterType.HIGHPASS:
        transition_bandwidth = transition_bandwidth_high
    elif filter_type in (FilterType.BANDPASS, FilterType.BANDSTOP):
        transition_bandwidth = min(transition_bandwidth_low, transition_bandwidth_high)
    else:
        raise ValueError("Unsupported filter type")

    delta_stopband = 10 ** (-0.05 * minimum_stopband_attenuation)
    delta_passband = (10 ** (0.05 * specified_passband_ripple) - 1) / (10 ** (0.05 * specified_passband_ripple) + 1)
    delta = min(delta_passband, delta_stopband)
    print("delta", delta)
    print("actual stopband attenuation", -20 * math.log10(delta_stopband))
    minimum_stopband_attenuation = -20 * math.log10(delta)

    if minimum_stopband_attenuation <= 21:
        parameter_d = 0.9222
    else:
        parameter_d = (minimum_stopband_attenuation - 7.95) / 14.36
        
    print("D: ", parameter_d)
    filter_order = int(2 + parameter_d * sampling_frequency / transition_bandwidth)
    if filter_order % 2 == 0:
        filter_order += 1  # Make it odd

    return (filter_order, delta, minimum_stopband_attenuation, parameter_d)

def kaiser_lowpass(passband_frequency_high,
    stopband_frequency_low,
    sampling_frequency,
    filter_order,
    kaiser_coeffs):

    cutoff_frequency = 0.5 * (passband_frequency_high + stopband_frequency_low)
    impulse_response = []
    n = filter_order
    #print("hd(n)")
    #print("-----")
    #for i in range((filter_order - 1) // 2 + 1):
    for i in range(filter_order):
        if i == 0:
            sinc = 2 * cutoff_frequency / sampling_frequency
        else:
            arg = 2 * cutoff_frequency * (i - (n - 1) // 2) / sampling_frequency
            #sinc = math.sin(arg) / (i - (n - 1) // 2)
            sinc = 2 * cutoff_frequency * np.sinc(arg)
        impulse_response.append(sinc * kaiser_coeffs[i])
        #print(sinc)
    impulse_response /= np.sum(impulse_response)
    print("Impulse response")
    print("----------------")
    for i in range(len(impulse_response)):
    	print(impulse_response[i])
    
    return impulse_response
    
def magnitude_response(omega, impulse_response, sampling_frequency, filter_order):
    frequency = omega / (2 * math.pi)
    period = 1 / sampling_frequency
    response_magnitude = 0.0
    #real_component = 0.0
    #imaginary_component = 0.0
    for i in range(len(impulse_response)):
        response_magnitude += impulse_response[i] * np.exp(-1j * 2 * math.pi * frequency * i * period)
        #arg = 2 * math.pi * frequency * i * period
        #real_component += impulse_response[i] * math.cos(arg)
        #imaginary_component -= impulse_response[i] * math.sin(arg)
    #response_magnitude = math.sqrt(real_component ** 2 + imaginary_component ** 2)
    response_magnitude = abs(response_magnitude)
    
    return response_magnitude 

# -------------------------------
# Example 5.4: Lowpass FIR Filter
# -------------------------------

FREQUENCY_RESOLUTION = 500

# Filter specifications
actual_passband_ripple = 0.1  # dB
minimum_stopband_attenuation = 44  # dB
passband_frequency = 500  # Hz
stopband_frequency = 750  # Hz
sampling_frequency = 2500  # Hz

filter_order, delta, minimum_stopband_attenuation, parameter_d = kaiser_filter_order(
    filter_type=FilterType.LOWPASS,
    passband_frequency_low=passband_frequency,
    passband_frequency_high=passband_frequency,
    stopband_frequency_low=stopband_frequency,
    stopband_frequency_high=stopband_frequency,
    sampling_frequency=sampling_frequency,
    specified_passband_ripple=actual_passband_ripple,
    minimum_stopband_attenuation=minimum_stopband_attenuation)                                        
                    
kaiser_coeffs, alpha = kaiser_coefficients(filter_order, minimum_stopband_attenuation)
impulse_response = kaiser_lowpass(passband_frequency,
    stopband_frequency,
    sampling_frequency,
    filter_order,
    kaiser_coeffs)

maximum_frequency = sampling_frequency / 2
frequencies = np.linspace(0, maximum_frequency, FREQUENCY_RESOLUTION)
angular_frequencies = frequencies * 2 * math.pi

response_magnitude = [magnitude_response(w, impulse_response, sampling_frequency, filter_order)
                      for w in angular_frequencies]
response_magnitude_db = [20 * math.log10(max(abs(m), 1e-10)) for m in response_magnitude]

# Plot
plt.figure(figsize=(8,4))
plt.plot(frequencies, response_magnitude_db)
plt.xlim(0, sampling_frequency / 2)
plt.ylim(-100, 5)
plt.axvline(passband_frequency, color='green', linestyle='--', label='Passband Edge')
plt.axvline(stopband_frequency, color='red', linestyle='--', label='Stopband Start')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.title(f"Kaiser Lowpass FIR (N={filter_order}, A={minimum_stopband_attenuation} dB)")
plt.grid(True)
plt.legend()
plt.show()

