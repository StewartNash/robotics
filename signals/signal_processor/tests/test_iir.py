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

# s-plane poles and zeros (lowpass)
# [Re(zero), Im(zero), Re(pole), Im(pole)]
s_pz = [
	[0.0000000, 0.0000000, -0.1017847, 1.0379357],
	[0.0000000, 0.0000000, -0.2780809, 0.7598216],
	[0.0000000, 0.0000000, -0.3798656, 0.2781140]
]

# z-plane poles and zeros (lowpass)
# [Re(zero), Im(zero), Re(pole), Im(pole)]
z_pz = [
	[-1.0000000, 0.0000000, 0.2472978, 0.8758249],
	[-1.0000000, 0.0000000, 0.3740355, 0.6310338],
	[-1.0000000, 0.0000000, 0.5290679, 0.2421385]
]

z_pz_conjugate = [
	[-1.0000000, 0.0000000, 0.2472978, -0.8758249],
	[-1.0000000, 0.0000000, 0.3740355, -0.6310338],
	[-1.0000000, 0.0000000, 0.5290679, -0.2421385]
]

# second-order section coefficients
# [stage, numerator coefficient A1, numerator coefficient A2,
# denominator coefficient B1, denominator coefficient B2]
sec_coeff = [
	[2.0000000, 1.0000000, -0.4945955, 0.8282254],
	[2.0000000, 1.0000000, -0.7480710, 0.5381062],
	[2.0000000, 1.0000000, -1.0581359, 0.3385440]
]

# frequency response output
# [frequency, magnitude]
freq_resp = [
	[0.000, -0.173],
	[20.00, -0.166],
	[40.000, -0.145],
	[60.000, -0.115],
	[80.000, -0.079],
	[100.00, -0.044]
]
