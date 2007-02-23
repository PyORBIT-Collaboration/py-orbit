"""
The module includes python utilities that can be useful everywhere.
"""

"""
This function creates multi-dimensional arrays like a[i][k][j].
The examples of using this function:
a = multiDimArray(5,10,2)
a = multiDimArray(*[5,10,2])
a[1][2][1] = 0
By default all elements are initialized by 0.
"""
def multiDimArray(*dims):
	res = []
	if len(dims) == 1:
		for j in xrange(dims[0]):
			res.append(0.)
	else:
		dims_rest = dims[1:len(dims)]
		for j in xrange(dims[0]):
			res.append(multiDimArray(*dims_rest))
	return res

"""
This function finializes the execution of the ORBIT's script.
"""
def orbitFinalize(message = None):
	import orbit_mpi
	import sys
	import traceback
	if(orbit_mpi.world_rank() == 0):
		print "ORBIT message:",message
		traceback.print_stack()
	sys.exit(1)
