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

# LOWPASS			LP
# HIGHPASS			HP
# BANDPASS			BP
# BANDSTOP			BS
# maximum_passband_attenuation	Ap	[positive dB]
# minimum_stopband_attenuation	As	[positive dB]
# passband_edge_frequency	fp
# stopband_edge_frequency	fs
# transition_bandwidth		fp <= f <= fs
# sampling_rate			F
# nyquist_frequency		F/2

# Filter Design Procedure
# 1. Enter filter specifications
# 	a. Filter type: LP, HP, BP or BS
# 	b. Filter parameters
# 		i. LP, HP: Ap, As, fp, fs, F
# 			A. LP: (fp < fs), (F > 2*fs)
# 			B. HP: (fp > fs), (F > 2*fp)
# 		ii. BP, BS: Ap, As, fp1, fp2, fs1, fs2, F/2
# 			A. BP: (fs1 < fp1 < fp2 < fs2), (F > 2*fs2)
# 			B. BS: (fp1 < fs1 < fs2 < fp2), (F > 2*fp2) 
# 2. Compute filter order, N (table 4.4)
# 3. Compute analog LP zeros
# 4. Compute analog LP poles
# 5. Compute digital poles and zeros
# 6. Compute second order section coefficients
# 7. Format coefficients as a function of section index, k
# 8. Compute coefficients B1,k for odd first-order section (N-odd only)
# 9. Determine second order section normalization coefficients


class FilterFamily(enum.Enum):
	BUTTERWORTH = 1
	CHEBYSHEV = 2
	ELLIPTIC = 3


class FilterType(enum.Enum):
	LOWPASS = 1
	HIGHPASS = 2
	BANDPASS = 3
	BANDSTOP = 4
	
def hplp_computations(maximum_passband_attenuation,
	minimum_stopband_attenuation,
	passband_edge_frequency,
	stopband_edge_frequency,
	sampling_frequency):
	
	Ap = maximum_passband_attenuation
	As = minimum_stopband_attenuation
	fp = passband_edge_frequency
	fs = stopband_edge_frequency
	omega_p = 2 * math.pi * fp
	omega_s = 2 * math.pi * fs
	parameter_epsilon = (10 ** (0.1 * Ap) - 1) ** 0.5
	parameter_lambda = (10 ** (0.1 * As) - 1) ** 0.5
	parameter_A = parameter_lambda / parameter_epsilon
	parameter_K0 = omega_p / omega_s
	
	return parameter_K0, parameter_A
	
def highpass_computations(maximum_passband_attenuation,
	minimum_stopband_attenuation,
	passband_edge_frequency,
	stopband_edge_frequency,
	sampling_frequency):
	
	parameter_K0, parameter_A = hplp_computations(maximum_passband_attenuation,
	minimum_stopband_attenuation,
	passband_edge_frequency,
	stopband_edge_frequency,
	sampling_frequency)
	parameter_K = 1 / parameter_K0
	
	A = parameter_A
	K0 = parameter_K0
	K = parameter_K
	
	return K, A
	
def lowpass_computations(maximum_passband_attenuation,
	minimum_stopband_attenuation,
	passband_edge_frequency,
	stopband_edge_frequency,
	sampling_frequency):
	
	parameter_K0, parameter_A = hplp_computations(maximum_passband_attenuation,
	minimum_stopband_attenuation,
	passband_edge_frequency,
	stopband_edge_frequency,
	sampling_frequency)
	parameter_K = parameter_K0
	
	A = parameter_A
	K0 = parameter_K0
	K = parameter_K
	
	return K, A
	
def bpbs_computations(fp1, fp2, fs1, fs2, F):
	KA = math.tan(math.pi * fp2 / F) - math.tan(math.pi * fp1 / F)
	KB = math.tan(math.pi * fp1 / F) * math.tan(math.pi * fp2 / F)
	KC = math.tan(math.pi * fs1 / F) * math.tan(math.pi * fs2 / F)
	K1 = (KA * math.tan(math.pi * fs1 / F)) / (KB - math.tan(math.pi * fs1 / F) ** 2)
	K2 = (KA * math.tan(math.pi * fs2 / F)) / (math.tan(math.pi * fs2 / F) ** 2 - KB)
	
	return K1, K2, KA, KB, KC
	
def bandpass_computations(fp1, fp2, fs1, fs2, F):
	K1, K2, KA, KB, KC = bpbs_computations(fp1, fp2, fs1, fs2, F)

	if KC >= KB:
		parameter_K = K1
	else: # KC < KB
		parameter_K = K2

	K = parameter_K	
	return K

def bandstop_computations(fp1, fp2, fs1, fs2, F):
	K1, K2, KA, KB, KC = bpbs_computations(fp1, fp2, fs1, fs2, F)

	if KC >= KB:
		parameter_K = 1 / K2
	else: # KC < KB
		parameter_K = 1 / K1
	
	K = parameter_K	
	return K
	
def compute_filter_order(parameter_K, parameter_A, filter_family):
	if filter_family == FilterFamily.BUTTERWORTH:
		filter_order = math.ceil(math.log(parameter_A) / math.log(1 / parameter_K))
	elif filter_family == FilterFamily.CHEBYSHEV:
		filter_order = math.ceil(math.acosh(parameter_A) / math.acosh(1 / parameter_K))
	elif filter_family == FilterFamily.ELLIPTIC:
		q, q0 = elliptic_computations(parameter_K)
		filter_order = math.ceil(math.log(16 * A) / math.log(1 / q))
	else:
		filter_order = 0
		
	N = filter_order
	return N

def elliptic_computations(parameter_K):
	K = parameter_K
	parameter_q0 = (1 - (1 - K ** 2)) ** 0.25 / (2 * (1 + (1 - K ** 2) ** 0.25))
	q0 = parameter_q0	
	parameter_q = q0 + 2 * qo ** 5 + 15 * q0 ** 9 + 150 * q0 ** 13
	q = parameter_q
	
	return q, q0
	
#def butterworth_analog_poles(filter_order, passband_angular_frequency=0, normalized=True):
def butterworth_analog_poles(filter_order, passband_angular_frequency=0):
	N = filter_order
	
	if N % 2: # Odd
		k = [i + 1 for i in range(int((N + 1) / 2))]
	else: # Even
		k = [i + 1 for i in range(int(N / 2))]
	
	poles = [(-math.sin(math.pi * (2 * k_ - 1) / (2 * N)), math.cos(math.pi * ((2 * k_) - 1) / (2 * N))) for k_ in k]
	#if passband_angular_frequency and not normalized:
	if passband_angular_frequency:
		omega = passband_angular_frequency
		poles = [(s[0] * omega ** (- 1 / N), s[1] * omega ** (- 1 / N)) for s in poles]
	
	return poles

