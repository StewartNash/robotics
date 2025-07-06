import math

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

# Kaiser order (NK) calculation
# Type: 1 - lowpass, 2 - highpass, 3 - bandpass, 4 - bandstop
# Passbands: FP1, FP2
# Stopbands: FA1, FA2
# Sampling rate: F
# AP in dB, AA in dB
# see design procedure section 5.3.2 steps 1 - 5
def kaiser_filter_order(filter_type, fp1, fp2, fa1, fa2, f, ap, aa):
	bt1  = abs(fp2 - fa2)
	bt2 = abs(fp1 - fa1)
	if filter_type == 1:
		bt = bt1
	if filter_type == 2:
		bt = bt2
	if filter_type == 3 or filter_type == 4:
		if bt1 < bt2:
			bt = bt1
		else:
			bt = bt2
	d2 = 10 ** (-0.05 * aa)
	d1 = (10 ** (0.05 * ap) - 1) / (10 ** (0.05 * ap) + 1)
	if d1 < d2:
		_del = d1
	else:
		_del = d2
	aap = -20 * math.log10(_del)) / math.log10(10)
	if aap <= 21:
		pard = 0.9222
	else:
		pard = (aap - 7.95) / 14.36
	nk = int(2 + pard * fs / bt)
	if (nk/2) == int(nk/2):
		nk = nk + 1
	return nk
	
# Compute Kaiser coefficients WK=AK eq. 5.52
# AAP from 15090, ALP using eq. 5.49
# IOBE and IOALP using eq. 5.41
# WK using eq. 5.39
# See design procedure section 5.3.2 steps 3, 6
# To complete eq. 5.52 branch to subroutine based on "type" eqs. given in section 5.3.3
def kaiser_coefficients(nk, AAP):
	AAP = 0
	#ALP = 0
	WK = []    
	
	n = nk
	if AAP <= 21.0:
		ALP = 0
	else:
		ALP = 0.1102 * (AAP - 8.7)
	if AAP > 21.0 and AAP <= 50.0:
		ALP = (0.5842 * (AAP - 21) ** 0.4) + 0.07886 *	(AAP - 21)
	KFAC = [1]
	for k in range(2 - 1, 30 + 1 - 1):
		KFAC.append(KFAC[k - 1] * k)
	print ("Computing Kaiser coefficients")
	for i in range(0, (NK - 1) / 2 + 1):
		print(str(i) + " out of " + str((nk - 1 ) / 2))
		BE = ALP * math.sqrt(1 - (2 * i / (NK - 1)) ** 2)
		IOBE = 1
		IOALP = 1
		for k in range(1, 30 + 1):
			IOBE = IOBE + (((BE / 2) ** k) / KFAC[k]) ** 2
			IOALP = IOALP + (((ALP / 2) ** k) / KFAC[k]) ** 2
		WK.append(IOBE / IOALP)		
		
	return (AAP, ALP, IOBE, IOALP)
	
# Kaiser Lowpass Subroutine Eqs. 5.52, 5.56, & 5.57
def kaiser_lowpass(FP2, FA2, FS, WK):
	WC = 0.5 * (FP2 + FA2)
	fnsx = lambda x, y : math.sin(x * y * 2 * math.pi / FS) / (x * y * 2 * math.pi / FS)
	H1 = 2 * WC / FS
	H = [H1 * WK[0]]
	for i in range(1, (NK - 1) / 2 + 1):
		H.append(H1 * fnsx(WC, i) * WK[i + 1 - 1]
	
	return H

# Kaiser Highpass Subroutine Eqs. 5.52, 5.58, & 5.59
def kaiser_highpass(FP1, FA1, FS, WK):
	WC = 0.5 * (FP1 + FA1)
	fnsx = lambda x, y : math.sin(x * y * 2 * math.pi / FS) / (x * y * 2 * math.pi / FS)
	H1 = -2 * WC / FS
	H = [(1 + H1) * WK[0]]
	for i in range(1, (NK - 1) / 2 + 1):
		H.append(H1 * fnsx(WC, i) * WK[i])
	
	return H
	
# Kaiser Bandpass Subroutine Eqs. 5.52, 5.60, & 5.61
def kaiser_bandpass(FP1, FP2, BT, WK):
	WC1 = FP1 - BT / 2
	WC2 = FP2 + BT / 2
	H = [(2 / FS) * (WC2 - WC1) * WK[0]]
	for i in range(1, (NK - 1) / 2 + 1):
		ARG = i * 2 * math.pi / FS
		H.append(1 / (math.pi * i) * math.sin(WC2 * ARG) - math.sin(WC1 * ARG) * WK[i])
		print("i = " + str(i) + " H(i) = " + str(H[i])
		
	return H
	
# Kaiser bandstop subroutine Eqs. 5.52, 5.62, & 5.63
def kaiser_bandstop(FP1, FP2, BT):
	WC1 = FP1 + BT / 2
	WC2 = FP2 - BT / 2
	H = [(2 * (WC1 - WC2) / FS + 1) * WK[0]]
	for i in range(1, (NK - 1) / 2 + 1):
		ARG = i * 2 * math.pi / FS
		H.append(1 / (math.pi * i) * (math.sin(WC1 * ARG) - math.sin(WC2 * ARG)) * WK[i]
		
	return H
	
