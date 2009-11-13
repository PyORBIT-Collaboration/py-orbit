def multiDimDoubleArray(*dims):
	"""
	Method. Creates multi-dimensional arrays with doubles, such as a[i][k][j].
	Some examples of the use of this function:
		a = multiDimArray(5,10,2)
		a = multiDimArray(*[5,10,2])
		a[1][2][1] = 0.
		By default all elements are initialized to 0.
	"""
	res = []
	if len(dims) == 1:
		for j in xrange(dims[0]):
			res.append(0.)
	else:
		dims_rest = dims[1:len(dims)]
		for j in xrange(dims[0]):
			res.append(multiDimDoubleArray(*dims_rest))
	return res
		
def multiDimIntArray(*dims):
	"""
	Method. Creates multi-dimensional arrays with integers, such as a[i][k][j].
	Some examples of the use of this function:
		a = multiDimArray(5,10,2)
		a = multiDimArray(*[5,10,2])
		a[1][2][1] = 0
		By default all elements are initialized to 0.
	"""
	res = []
	if len(dims) == 1:
		for j in xrange(dims[0]):
			res.append(0)
	else:
		dims_rest = dims[1:len(dims)]
		for j in xrange(dims[0]):
			res.append(multiDimIntArray(*dims_rest))
	return res
