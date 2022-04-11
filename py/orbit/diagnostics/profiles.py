import math

from bunch import Bunch

#pyORBIT MPI module import
import orbit_mpi
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

def profiles(Bunch, coord, histogram, steps = 100, Min = 1.0, Max = -1.0):
	"""
        Returns a profile for one of the following Bunch coordinates:
	x[m] xp[rad] y[m] yp[rad] z[m] dE[GeV]
	"""

	b = Bunch

	# Take the MPI Communicator from bunch: It could be
        # different from MPI_COMM_WORLD

	comm = Bunch.getMPIComm()
	rank = orbit_mpi.MPI_Comm_rank(comm)
	size = orbit_mpi.MPI_Comm_size(comm)
	main_rank = 0

	# n_parts_arr - array of size of the number of CPUs,
	# contains the number of macroparticles on each CPU

	n_parts_arr = [0]*size
	n_parts_arr[rank] = b.getSize()
	n_parts_arr = orbit_mpi.MPI_Allreduce(n_parts_arr, \
                mpi_datatype.MPI_INT,mpi_op.MPI_SUM,comm)

        partdat = []

	if(rank == main_rank):
		for i in range(n_parts_arr[rank]):
                        if coord == "x":
                                partdat.append(b.x(i))
                        if coord == "px":
                                partdat.append(b.px(i))
                        if coord == "y":
                                partdat.append(b.y(i))
                        if coord == "py":
                                partdat.append(b.py(i))
                        if coord == "z":
                                partdat.append(b.z(i))
                        if coord == "dE":
                                partdat.append(b.dE(i))

	# That is just for case.
        # Actually, MPI_Barrier command is not necessary.

	orbit_mpi.MPI_Barrier(comm)

	val_arr = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

	for i_cpu in range(1,size):

		# Again, that is just for case.
                # Actually, MPI_Barrier command is not necessary.
		orbit_mpi.MPI_Barrier(comm)

		for i in range(n_parts_arr[i_cpu]):
			if(rank == main_rank):

				#get the coordinate array
				(x, px, y, py, z, dE) = \
                                orbit_mpi.MPI_Recv(mpi_datatype.MPI_DOUBLE, \
                                i_cpu, 222, comm)
                                if coord == "x":
                                        partdat.append(x)
                                if coord == "px":
                                        partdat.append(px)
                                if coord == "y":
                                        partdat.append(y)
                                if coord == "py":
                                        partdat.append(py)
                                if coord == "z":
                                        partdat.append(z)
                                if coord == "dE":
                                        partdat.append(dE)

			elif(rank == i_cpu):
				#send the coordinate array
				x  = b.x(i)
				px = b.px(i)
				y  = b.y(i)
				py = b.py(i)
				z  = b.z(i)
				dE = b.dE(i)
				val_arr = (x, px, y, py, z, dE)
				orbit_mpi.MPI_Send(val_arr, \
                                mpi_datatype.MPI_DOUBLE, main_rank, 222, comm)

        l = len(partdat)
        m = min(partdat)
        M = max(partdat)

        c = (M + m) / 2.0
        d = (M - m) * 1.1 / 2.0
        M = c + d
        m = c - d

        if Max > M:
                M = Max
        if Min < m:
                m = Min

        dx = (M - m) / steps

        grid = [m]
        prof = [0]
        for i in range(1, steps + 1):
                x = m + i * dx
                grid.append(x)
                prof.append(0)
        grid.append(M)
        prof.append(0)

        for n in range(l):
                i = (partdat[n] - m) / dx
                i = int(i)
                if i < 0:
                        pass
                elif i > range(steps):
                        pass
                else:
                        frac = (partdat[n] - m) / dx % 1
                        prof[i] = prof[i] + (1.0 - frac)
                        prof[i + 1] = prof[i + 1] + frac

        sum = 0.0
        for i in range(steps + 1):
                sum = sum + prof[i]

	file_out = histogram
	if(rank == main_rank):
		file_out = open(histogram, "w")
                
                file_out.write("Min = " + str(m) + "  Max = " + \
                               str(M) + " steps = " + str(steps) + "\n")
                file_out.write("nParts = " + str(l) + " HistSum = " + \
                               str(sum) + "\n\n")

                for i in range(steps + 1):
                        file_out.write(str(grid[i]) + "   " + \
                                       str(prof[i]) + "\n")

	if(rank == main_rank):
		file_out.close()
