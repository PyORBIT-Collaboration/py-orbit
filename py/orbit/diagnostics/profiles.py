import os
import math

from spacecharge import Grid1D
from orbit_utils import BunchExtremaCalculator

#pyORBIT MPI module import
import orbit_mpi
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

def profiles(bunch, coord, histogram, steps = 100, Min = 1.0, Max = -1.0):
	"""
	Returns a profile as Grid1D object for one of the following bunch coordinates:
	x[m] xp[rad] y[m] yp[rad] z[m] dE[GeV]
	"""
	
	grid1D = Grid1D(steps)
	
	# Take the MPI Communicator from bunch: It could be different from MPI_COMM_WORLD
	comm = bunch.getMPIComm()
	rank = orbit_mpi.MPI_Comm_rank(comm)
	size = orbit_mpi.MPI_Comm_size(comm)
	main_rank = 0
	
	nParts = bunch.getSizeGlobal()
	if(nParts < 2):
		if(rank == main_rank):
			file_out = open(histogram, "w")
			file_out.write("bunch has only number of particles =" +str(nParts)+os.linesep)
			file_out.close()
			return grid1D

	coord_index = -1
	if coord == "x" : coord_index = 0
	if coord == "px": coord_index = 1
	if coord == "y" : coord_index = 2
	if coord == "py": coord_index = 3
	if coord == "z" : coord_index = 4
	if coord == "dE": coord_index = 5
	if(coord_index < 0):
		st  = "Error: profiles(bunch,coord,...)" + os.linesep
		st += "coord is not x,xp,y,yp,z, or dE." + os.linesep
		st += "Stop." + os.linesep
		orbit_mpi.finalize(st)
		return None
	
	#---- calculations Min Max values
	bunch_extrema_cal = BunchExtremaCalculator()

	if(Min >= Max):
		#---- We are not going to use 
		bunch_extrema_cal = BunchExtremaCalculator()
		(xMin,xMax,yMin,yMax,zMin,zMax) = bunch_extrema_cal.extremaXYZ(bunch)
		(xpMin,xpMax,ypMin,ypMax,dE_Min,dE_Max) = bunch_extrema_cal.extremaXpYpdE(bunch)
		valMinMax_arr = []
		valMinMax_arr.append([xMin,xMax])
		valMinMax_arr.append([xpMin,xpMax])
		valMinMax_arr.append([yMin,yMax])
		valMinMax_arr.append([ypMin,ypMax])
		valMinMax_arr.append([zMin,zMax])
		valMinMax_arr.append([dE_Min,dE_Max])
		[Min,Max] = valMinMax_arr[coord_index]
		
	if(Min == Max):
		if(rank == main_rank):
			file_out = open(histogram, "w")
			file_out.write("No histogram min=" + str(Min) + " max=" + str(Max) + os.linesep)
			file_out.close()
			return grid1D
	
	#---- perform bunch binning
	grid1D.setGridZ(Min,Max)
	grid1D.binBunch(bunch,coord_index)
	grid1D.synchronizeMPI(comm)
	sum_value = grid1D.getSum()
	
	#--- dump the histogram to the file
	if(rank == main_rank):
		file_out = open(histogram, "w")
		st  = "% "+" Min, Max = %15.8g , %15.8g "%(Min,Max)
		st += " nSteps = %5d "%steps
		st += " hist_sum = %15.8g"%sum_value
		file_out.write(st + os.linesep)
		for ind in range(steps):
			val = grid1D.getGridZ(ind)
			rho = grid1D.getValueOnGrid(ind)
			st = " %15.8g  %15.8g "%(val,rho)
			file_out.write(st + os.linesep)
		file_out.close()

	return grid1D