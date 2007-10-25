def orbitFinalize(message = None):
	"""
	Method. Finalizes the execution of the ORBIT script.
	"""
	import orbit_mpi
	import sys
	import traceback
	if(orbit_mpi.world_rank() == 0):
		print "ORBIT message: ", message
		traceback.print_stack()
	print "STOP"
	sys.exit(1)

