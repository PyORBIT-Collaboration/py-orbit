import sys
import orbit_mpi

cpu  = orbit_mpi.cpu()
rank = orbit_mpi.world_rank()
size = orbit_mpi.world_size()
tm = orbit_mpi.time()
print "rank=",rank," size=",size, " name=", cpu," time=",tm

orbit_mpi.finalize("Error Message!")
