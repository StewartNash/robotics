
def polynomial_coefficients(roots):
	"""
	Compute polynomial coefficients for z^-1 form given roots.
	Roots can be complex. Returns coefficients:
	c[0] + c[1] z^-1 + ... + c[N] z^-N
	"""
	coefficients = [(1.0, 0.0)]
	for r in roots:
		new_coefficients = [(0.0, 0.0)] * (len(coefficients) + 1)
		for i in range(len(coefficients)):
			new_coefficients[i] = complex_add(new_coefficients[i], coefficients[i])
			new_coefficients[i + 1] = complex_subtract(new_coefficients[i + 1], complex_multiply(r, coefficients[i]))
		coefficients = new_coefficients
	
	return coefficients
	
def complex_multiply(first, second):
	a = first[0] # Real
	b = first[1] # Imaginary
	c = second[0] # Real
	d = second[1] # Imaginary
	
	return (a * c - b * d,  a * d + b * c)
	
def complex_add(first, second):
	a = first[0] # Real
	b = first[1] # Imaginary
	c = second[0] # Real
	d = second[1] # Imaginary
	
	return (a + c,  b + d)

def complex_subtract(first, second):
	a = first[0] # Real
	b = first[1] # Imaginary
	c = second[0] # Real
	d = second[1] # Imaginary
	
	return (a - c,  b - d)
	
def real_to_complex(number_list):
	output = [(r, 0.0)  for r in number_list]
	
	return output
