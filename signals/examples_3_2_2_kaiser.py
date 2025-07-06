import math
import enum

# filter_order			nk, NK
# passband_frequency		fp
# passband_frequency_low	fp1, FP1
# passband_frequency_high	fp2, FP2
# stopband_frequency		fs, FS, f, F
# sampling_frequency		fs, FS, f, F (CORRECTION)
# stopband_frequency_low	fs1, fa1, FA1
# stopband_frequency_high	fs2, fa2, FA2
# transition_bandwidth_low	bt1
# transition_bandwidth_high	bt2
# transition_bandwidth		bt
# delta				_del
# delta_passband (ripple)	d1
# delta_stopband (ripple)	d2
# actual_passband_ripple	aap
# specified_passband_ripple	ap
# actual_stopband_attenuation	actual_minimum_stopband_attenuation
# minimum_stopband_attenuation	specified_minimum_stopband_attenuation
# minimum_stopband_attenuation	aa
# parameter_d			pard
# alpha				alp, ALP

# actual_passband_ripple	aa (CORRECTION)
# minimum_stopband_attenuation	aap (CORRECTION)
# kaiser_coeffs			wk, WK
# mod_bessel_fk_alpha		modified_bessel_first_kind_alpha
# mod_bessel_fk_beta		modified_bessel_first_kind_beta
# mod_bessel_fk_alpha		IOALP
# mod_bessel_fk_beta		IOBE
# beta				BE
#				KFAC (factorial)
# cutoff_frequency		wc, WC
# cutoff_frequency_low		wc1, WC1
# cutoff_frequency_high		wc2, WC2
# initial_impulse_response	H1
# impulse_response		H
# sinc_function			fnsx, FNSX
# argument			ARG

class FilterType(enum.Enum):
	LOW_PASS = 1
	HIGH_PASS = 2
	BAND_PASS = 3
	BAND_STOP = 4
	LOWPASS = 5
	HIGHPASS = 6
	BANDPASS = 7
	BANDSTOP = 8

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
	if filter_type == FilterType.LOW_PASS or filter_type == FilterType.LOWPASS:
		transition_bandwidth = transition_bandwidth_low
	if filter_type == FilterType.HIGH_PASS or filter_type == FilterType.HIGHPASS:
		transition_bandwidth = transition_bandwidth_high
	if filter_type == FilterType.BAND_PASS or filter_type == FilterType.BAND_STOP:
		transition_bandwidth = min(transition_bandwidth_low, transition_bandwidth_high)
	if filter_type == FilterType.BANDPASS or filter_type == FilterType.BANDSTOP:
		transition_bandwidth = min(transition_bandwidth_low, transition_bandwidth_high)		
	delta_stopband = 10 ** (-0.05 * minimum_stopband_attenuation)
	delta_passband = (10 ** (0.05 * specified_passband_ripple) - 1) / (10 ** (0.05 * specified_passband_ripple) + 1)
	delta = min(delta_passband, delta_stopband)
	#actual_passband_ripple = -20 * math.log10(delta)
	minimum_stopband_attenuation = -20 * math.log10(delta) # (CORRECTION)
	#if actual_passband_ripple <= 21:
	#	parameter_d = 0.9222
	#else:
	#	parameter_d = (actual_passband_ripple - 7.95) / 14.36
	if minimum_stopband_attenuation <= 21: # (CORRECTION)
		parameter_d = 0.9222
	else:
		parameter_d = (minimum_stopband_attenuation - 7.95) / 14.36		
	filter_order = int(2 + parameter_d * sampling_frequency / transition_bandwidth)
	if (filter_order / 2) == int(filter_order / 2):
		filter_order = filter_order + 1

	return (filter_order, delta, minimum_stopband_attenuation, parameter_d)

def kaiser_coefficients(filter_order, actual_passband_ripple):
	kaiser_coeffs = []
	
	nk = filter_order
	n = filter_order
	if actual_passband_ripple <= 21.0:
		alpha = 0
	else:
		alpha = 0.1102 * (actual_passband_ripple - 8.7)
	if actual_passband_ripple > 21.0 and actual_passband_ripple <= 50.0:
		alpha = (0.5842 * (actual_passband_ripple - 21) ** 0.4) + 0.07886 *  (actual_passband_ripple - 21)
	#KFAC = [1]
	#for k in range(2 - 1, 30 + 1 - 1):
	#	KFAC.append(KFAC[k - 1] * k)
	print ("Computing Kaiser coefficients")
	for i in range(0, int((nk - 1) / 2) + 1):
		print(str(i) + " out of " + str((nk - 1 ) / 2))
		beta = alpha * math.sqrt(1 - (2 * i / (nk - 1)) ** 2)
		mod_bessel_fk_beta = 1
		mod_bessel_fk_alpha = 1
		for k in range(1, 30 + 1):
			#mod_bessel_fk_beta = mod_bessel_fk_beta + (((beta / 2) ** k) / KFAC[k]) ** 2
			#mod_bessel_fk_alpha = mod_bessel_fk_alpha + (((alpha / 2) ** k) / KFAC[k]) ** 2
			mod_bessel_fk_beta = mod_bessel_fk_beta + (((beta / 2) ** k) / math.factorial(k)) ** 2
			mod_bessel_fk_alpha = mod_bessel_fk_alpha + (((alpha / 2) ** k) / math.factorial(k)) ** 2
		kaiser_coeffs.append(mod_bessel_fk_beta / mod_bessel_fk_alpha)
		
	#return (actual_passband_ripple, alpha, IOBE, IOALP)
	return (kaiser_coeffs, alpha)

def kaiser_lowpass(passband_frequency_high,
	stopband_frequency_high,
	sampling_frequency,
	filter_order,
	kaiser_coeffs):
	nk = filter_order
	cutoff_frequency = 0.5 * (passband_frequency_high + stopband_frequency_high)
	sinc_function = lambda x, y : math.sin(x * y * 2 * math.pi / sampling_frequency) / (x * y * 2 * math.pi / sampling_frequency)
	initial_impulse_response = 2 * cutoff_frequency / sampling_frequency
	impulse_response = [initial_impulse_response * kaiser_coeffs[0]]
	for i in range(1, int((nk - 1) / 2) + 1):
		impulse_response.append(initial_impulse_response * sinc_function(cutoff_frequency, i) * kaiser_coeffs[i])
	
	return impulse_response

def kaiser_highpass(passband_frequency_low,
	stopband_frequency_low,
	sampling_frequency,
	filter_order,
	kaiser_coeffs):
	nk = filter_order
	cutoff_frequency = 0.5 * (passband_frequency_low + stopband_frequency_low)
	sinc_function = lambda x, y : math.sin(x * y * 2 * math.pi / sampling_frequency) / (x * y * 2 * math.pi / sampling_frequency)
	initial_impulse_response = -2 * cutoff_frequency / sampling_frequency
	impulse_response = [(1 + initial_impulse_response) * kaiser_coeffs[0]]
	for i in range(1, int((nk - 1) / 2) + 1):
		impulse_response.append(initial_impulse_response * sinc_function(cutoff_frequency, i) * kaiser_coeffs[i])
	
	return impulse_response
	
def kaiser_bandpass(passband_frequency_low,
	passband_frequency_high,
	sampling_frequency,
	transition_bandwidth,
	filter_order,
	kaiser_coeffs):
	nk = filter_order
	cutoff_frequency_low = passband_frequency_low - transition_bandwidth / 2
	cutoff_frequency_high = passband_frequency_high + transition_bandwidth / 2
	impulse_response = [(2 / sampling_frequency) * (cutoff_frequency_high - cutoff_frequency_low) * kaiser_coeffs[0]]
	for i in range(1, int((nk - 1) / 2) + 1):
		argument = i * 2 * math.pi / sampling_frequency
		impulse_response.append(1 / (math.pi * i) * math.sin(cutoff_frequency_high * argument) - math.sin(cutoff_frequency_low * argument) * kaiser_coeffs[i])
		print("i = " + str(i) + " H(i) = " + str(impulse_response[i]))
		
	return impulse_response

def kaiser_bandstop(passband_frequency_low,
	passband_frequency_high,
	sampling_frequency,
	transition_bandwidth,
	filter_order,
	kaiser_coeffs):
	nk = filter_order
	cutoff_frequency_low = passband_frequency_low + transition_bandwidth / 2
	cutoff_frequency_high = passband_frequency_high - transition_bandwidth / 2
	impulse_response = [(2 * (cutoff_frequency_low - cutoff_frequency_high) / sampling_frequency + 1) * kaiser_coeffs[0]]
	for i in range(1, int((nk - 1) / 2) + 1):
		argument = i * 2 * math.pi / sampling_frequency
		impulse_response.append(1 / (math.pi * i) * (math.sin(cutoff_frequency_low * argument) - math.sin(cutoff_frequency_high * argument)) * kaiser_coeffs[i])
		
	return impulse_response

