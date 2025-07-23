import math
import numpy as np
import matplotlib.pyplot as plt
import enum

class FilterType(enum.Enum):
    LOWPASS = 1
    HIGHPASS = 2
    BANDPASS = 3
    BANDSTOP = 4

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
    minimum_stopband_attenuation = -20 * math.log10(delta)

    if minimum_stopband_attenuation <= 21:
        parameter_d = 0.9222
    else:
        parameter_d = (minimum_stopband_attenuation - 7.95) / 14.36

    filter_order = int(2 + parameter_d * sampling_frequency / transition_bandwidth)
    if filter_order % 2 == 0:
        filter_order += 1  # Make it odd

    return (filter_order, delta, minimum_stopband_attenuation, parameter_d)

def kaiser_coefficients(filter_order, actual_passband_ripple):
    n = filter_order
    if actual_passband_ripple <= 21.0:
        alpha = 0
    elif actual_passband_ripple <= 50.0:
        alpha = 0.5842 * (actual_passband_ripple - 21) ** 0.4 + 0.07886 * (actual_passband_ripple - 21)
    else:
        alpha = 0.1102 * (actual_passband_ripple - 8.7)

    kaiser_coeffs = []
    for i in range((n - 1) // 2 + 1):
        beta = alpha * math.sqrt(1 - (2 * i / (n - 1)) ** 2)
        mod_bessel_fk_beta = 1
        mod_bessel_fk_alpha = 1
        for k in range(1, 31):
            mod_bessel_fk_beta += (((beta / 2) ** k) / math.factorial(k)) ** 2
            mod_bessel_fk_alpha += (((alpha / 2) ** k) / math.factorial(k)) ** 2
        kaiser_coeffs.append(mod_bessel_fk_beta / mod_bessel_fk_alpha)

    return (kaiser_coeffs, alpha)

def kaiser_lowpass(passband_frequency_high,
    stopband_frequency_high,
    sampling_frequency,
    filter_order,
    kaiser_coeffs):

    cutoff_frequency = passband_frequency_high  # Use passband edge as cutoff
    impulse_response = []
    for i in range((filter_order - 1) // 2 + 1):
        if i == 0:
            sinc = 2 * cutoff_frequency / sampling_frequency
        else:
            arg = 2 * math.pi * cutoff_frequency * i / sampling_frequency
            sinc = math.sin(arg) / (math.pi * i)
        impulse_response.append(sinc * kaiser_coeffs[i])

    # Make impulse response symmetric
    impulse_response = impulse_response[::-1][1:] + impulse_response
    return impulse_response

def magnitude_response(omega, impulse_response, sampling_frequency, filter_order):
    frequency = omega / (2 * math.pi)
    period = 1 / sampling_frequency
    response_magnitude = 0.0
    for i in range(len(impulse_response)):
        response_magnitude += impulse_response[i] * math.cos(2 * math.pi * frequency * i * period)
    return response_magnitude

# -------------------------------
# Example 5.4: Lowpass FIR Filter
# -------------------------------

FREQUENCY_RESOLUTION = 500

# Filter specs
actual_passband_ripple = 0.1  # dB
minimum_stopband_attenuation = 44  # dB
passband_frequency = 500  # Hz
stopband_frequency = 750  # Hz
sampling_frequency = 2500  # Hz

(filter_order, delta, minimum_stopband_attenuation, parameter_d) = kaiser_filter_order(
    filter_type=FilterType.LOWPASS,
    passband_frequency_low=passband_frequency,
    passband_frequency_high=passband_frequency,
    stopband_frequency_low=stopband_frequency,
    stopband_frequency_high=stopband_frequency,
    sampling_frequency=sampling_frequency,
    specified_passband_ripple=actual_passband_ripple,
    minimum_stopband_attenuation=minimum_stopband_attenuation)

print("Filter order:", filter_order)

(kaiser_coeffs, alpha) = kaiser_coefficients(
    filter_order=filter_order,
    actual_passband_ripple=actual_passband_ripple)

impulse_response = kaiser_lowpass(
    passband_frequency_high=passband_frequency,
    stopband_frequency_high=stopband_frequency,
    sampling_frequency=sampling_frequency,
    filter_order=filter_order,
    kaiser_coeffs=kaiser_coeffs)

# Compute and plot response
maximum_frequency = sampling_frequency / 2
frequencies = np.linspace(0, maximum_frequency, FREQUENCY_RESOLUTION)
angular_frequencies = frequencies * 2 * math.pi
response_magnitude = [magnitude_response(w, impulse_response, sampling_frequency, filter_order)
                      for w in angular_frequencies]
response_magnitude_db = [20 * math.log10(max(abs(m), 1e-10)) for m in response_magnitude]

plt.plot(frequencies, response_magnitude_db)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.title("Lowpass Filter Magnitude Response")
plt.grid(True)
plt.ylim([-100, 5])
plt.axvline(x=passband_frequency, color='green', linestyle='--', label='Passband Edge')
plt.axvline(x=stopband_frequency, color='red', linestyle='--', label='Stopband Start')
plt.legend()
plt.show()

