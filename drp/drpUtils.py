# Miscellaneous code
# Robert Lasenby 2009

from itertools import izip_longest

def unzip (l):
	"""Transposes the first two levels of a list of list (resp tuples etc.)"""
	c = zip(*l)
	return [list(a) for a in c]

def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)
