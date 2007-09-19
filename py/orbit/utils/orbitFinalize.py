"""
This function finalizes the execution of the ORBIT script.
"""

def orbitFinalize(message = None):
	import orbit_mpi
	import sys
	import traceback
	if(orbit_mpi.world_rank() == 0):
		print "ORBIT message: ", message
		traceback.print_stack()
	sys.exit(1)
