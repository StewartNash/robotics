import math

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
def kaiser_coefficients():
    AAP = 0
    ALP = 0
    IOBE = 0
    IOALP = 0
    return (AAP, ALP, IOBE, IOALP)
