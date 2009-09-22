# Miscellaneous code
# Robert Lasenby 2009

def unzip (l):
	"""Transposes the first two levels of a list of list (resp tuples etc.)"""
	c = zip(*l)
	return [list(a) for a in c]
