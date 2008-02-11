def orbitFinalize(message = None):
	"""
	Method. Finalizes the execution of the ORBIT script.
	"""
	import orbit_mpi
	import sys
	import traceback
	if(orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD) == 0):
		print "ORBIT message: ", message
		traceback.print_stack()
	print "STOP"
	sys.exit(1)

