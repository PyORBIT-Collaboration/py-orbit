def multiDimArray(*dims):
	"""
	Method. Creates multi-dimensional arrays, such as a[i][k][j].
	Some examples of the use of this function:
		a = multiDimArray(5,10,2)
		a = multiDimArray(*[5,10,2])
		a[1][2][1] = 0
		By default all elements are initialized to 0.
	"""
	res = []
	if len(dims) == 1:
		for j in xrange(dims[0]):
			res.append(0.)
	else:
		dims_rest = dims[1:len(dims)]
		for j in xrange(dims[0]):
			res.append(multiDimArray(*dims_rest))
	return res

