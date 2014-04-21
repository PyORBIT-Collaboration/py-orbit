import math

#pyORBIT MPI module import
import orbit_mpi
from orbit_mpi import mpi_datatype

from bunch import Bunch

def bunch_orbit_to_pyorbit_nHarm(ringLength, nHarm, kineticEnergy, \
	name_of_orbit_mpi_bunch_file, pyOrbitBunch = None):
	"""
	Translates ORBIT_MPI bunch to pyORBIT bunch, incorporating
	RF harmonic number, and returns it.
	PyORBIT bunch needs the ring length (m) and energy, mass, and
	charge of the synchronous particle, but ORBIT_MPI does not 
	have it. So, this information is specified in pyOrbitBunch or
	it will be proton by default.
	Lines in bunch files:
	ORBIT_MPI: x[mm] xp[mrad] y[mm] yp[mrad] phi[rad] dE[GeV].
	pyORBIT:   x[m]  xp[rad]  y[m]  yp[rad]  z[m]     dE[GeV]
	"""
	zfac = ringLength / (2 * math.pi * nHarm)
	if(pyOrbitBunch == None):  pyOrbitBunch = Bunch()
	# Take the MPI Communicator from bunch: it could be different
	# from MPI_COMM_WORLD
	comm = pyOrbitBunch.getMPIComm()
	rank = orbit_mpi.MPI_Comm_rank(comm)
	size = orbit_mpi.MPI_Comm_size(comm)	
	main_rank = 0
	# We will operate on file only at the CPU with rank = 0
	file_in = None
	if(rank == main_rank):
		file_in = open(name_of_orbit_mpi_bunch_file,"r")

	pyOrbitBunch.getSyncParticle().kinEnergy(kineticEnergy)

	# Here we assign ln_nonempty = 1 for each CPU if the line
	# read by CPU with the rank = 0 was non-empty
	# Otherwise, ln_nonempty will be zero everywhere
	ln = None
	ln_nonempty = 0

	if(rank == main_rank):
		ln_nonempty = 0
		ln = file_in.readline().strip()
		if(len(ln) > 0):
			ln_nonempty = 1
	ln_nonempty = orbit_mpi.MPI_Bcast(ln_nonempty, \
		mpi_datatype.MPI_INT,main_rank,comm)

	var_arr = (0., 0., 0., 0., 0., 0.)
	
	n_count = 1

	while(ln_nonempty == 1):
		# The rank of CPU that will get the next particle.
		# Loops through all CPUs. 
		recv_rank = n_count % size	
		if(rank == main_rank):
			res_arr = ln.strip().split()
			x  =  float(res_arr[0]) / 1000.
			xp =  float(res_arr[1]) / 1000.
			y  =  float(res_arr[2]) / 1000.
			yp =  float(res_arr[3]) / 1000.
			z  = -float(res_arr[4]) * zfac
			dE =  float(res_arr[5])
			val_arr = (x, xp, y, yp, z, dE)
			# Send the information if rank = 0 is not
			# going to keep this particle
			if(recv_rank != main_rank):
				orbit_mpi.MPI_Send(val_arr, \
					mpi_datatype.MPI_DOUBLE, \
					recv_rank, 111, comm)
			else:
				pyOrbitBunch.addParticle(val_arr[0], \
					val_arr[1], val_arr[2], \
					val_arr[3], val_arr[4], \
					val_arr[5])
		if(rank == recv_rank and rank != main_rank):
			val_arr = orbit_mpi.MPI_Recv(\
				mpi_datatype.MPI_DOUBLE, \
				main_rank, 111, comm)
			pyOrbitBunch.addParticle(val_arr[0], \
				val_arr[1], val_arr[2], \
				val_arr[3], val_arr[4], \
				val_arr[5])
		# Let's again find out if we want to proceed 
		if(rank == main_rank):
			ln_nonempty = 0
			ln = file_in.readline().strip()
			if(len(ln) > 0):
				ln_nonempty = 1
		ln_nonempty = orbit_mpi.MPI_Bcast(ln_nonempty, \
			mpi_datatype.MPI_INT, main_rank, comm)				
		n_count += 1		

	if(rank == main_rank): file_in.close()
	return pyOrbitBunch

