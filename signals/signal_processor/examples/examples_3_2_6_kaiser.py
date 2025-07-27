import math
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Bessel & Kaiser Window
# -------------------------
def I0_bessel(x, terms=50):
    """Approximate I0(x) using a power series."""
    result = 1.0
    y = (x/2)**2
    term = 1.0
    for k in range(1, terms):
        term *= y / (k*k)
        result += term
    return result

def kaiser_window(N, alpha):
    """Full Kaiser window (no SciPy)."""
    denom = I0_bessel(alpha)
    window = np.zeros(N)
    for n in range(N):
        ratio = (2*n)/(N-1) - 1.0
        window[n] = I0_bessel(alpha * math.sqrt(1 - ratio*ratio)) / denom
    return window

# -------------------------
# Filter Order & Alpha (Your Method)
# -------------------------
def kaiser_filter_order_textbook(passband_frequency,
                                 stopband_frequency,
                                 sampling_frequency,
                                 passband_ripple_db,
                                 stopband_atten_db):
    """Return (N, alpha) using textbook method you used originally."""
    delta_f = abs(stopband_frequency - passband_frequency)
    A = stopband_atten_db

    # D parameter
    if A <= 21:
        D = 0.9222
    else:
        D = (A - 7.95) / 14.36

    # Filter order (rounded to odd)
    N = int(2 + D * sampling_frequency / delta_f)
    if N % 2 == 0:
        N += 1

    # Kaiser alpha
    if A <= 21.0:
        alpha = 0
    elif A <= 50.0:
        alpha = 0.5842 * (A - 21) ** 0.4 + 0.07886 * (A - 21)
    else:
        alpha = 0.1102 * (A - 8.7)

    return N, alpha

# -------------------------
# Kaiser Lowpass Filter
# -------------------------
def kaiser_lowpass(passband_frequency_high,
                   stopband_frequency_low,
                   sampling_frequency,
                   filter_order,
                   alpha):
    # Normalized cutoff (0 to 0.5)
    fc = 0.5 * (passband_frequency_high + stopband_frequency_low) / sampling_frequency

    # Ideal sinc response centered at N/2
    M = filter_order - 1
    n = np.arange(filter_order)
    h_ideal = 2 * fc * np.sinc(2 * fc * (n - M/2))

    # Full Kaiser window
    window = kaiser_window(filter_order, alpha)

    # Windowed impulse response
    h = h_ideal * window

    # Normalize for unity DC gain
    h /= np.sum(h)
    return h

# -------------------------
# Manual Magnitude Response
# -------------------------
def magnitude_response_manual(impulse_response, frequencies, sampling_frequency):
    """Compute magnitude response manually (no FFT)."""
    N = len(impulse_response)
    response_mag = []
    for f in frequencies:
        omega = 2 * math.pi * f / sampling_frequency
        re = 0.0
        im = 0.0
        for n in range(N):
            re += impulse_response[n] * math.cos(omega * n)
            im -= impulse_response[n] * math.sin(omega * n)
        mag = math.sqrt(re*re + im*im)
        response_mag.append(mag)
    return np.array(response_mag)

# -------------------------
# Example Usage
# -------------------------
# Filter specs
fs = 2500
passband_frequency = 500
stopband_frequency = 750
passband_ripple_db = 0.1
stopband_atten_db = 44

# Determine order & alpha (your method)
N, alpha = kaiser_filter_order_textbook(passband_frequency,
                                        stopband_frequency,
                                        fs,
                                        passband_ripple_db,
                                        stopband_atten_db)

# Design filter
h = kaiser_lowpass(passband_frequency, stopband_frequency, fs, N, alpha)

# Frequency response
frequencies = np.linspace(0, fs/2, 500)
response_mag = magnitude_response_manual(h, frequencies, fs)
response_db = 20*np.log10(np.maximum(response_mag, 1e-12))

# Plot
plt.figure(figsize=(8,4))
plt.plot(frequencies, response_db)
plt.xlim(0, fs/2)
plt.ylim(-100, 5)
plt.axvline(passband_frequency, color='green', linestyle='--', label='Passband Edge')
plt.axvline(stopband_frequency, color='red', linestyle='--', label='Stopband Start')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.title(f"Kaiser Lowpass FIR (N={N}, A={stopband_atten_db} dB)")
plt.grid(True)
plt.legend()
plt.show()

