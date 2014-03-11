#!/usr/bin/env python

"""
This is not a parallel version! 
"""
import math
import random
import sys
from bunch import Bunch
#from mpi import orbit_mpi
import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

class InjectParts:
	""" 
	This routine injects particles into a bunch with user specified distribution
	functions.
	"""
	
	def __init__(self, nparts, bunch, lostbunch, injectregion, xDistFunc, yDistFunc, lDistFunc, nmaxmacroparticles = -1, injectturninterval=1):
		self.nparts = nparts
		self.bunch = bunch
		self.lostbunch = lostbunch
		self.injectregion = injectregion
		self.xDistFunc = xDistFunc
		self.yDistFunc = yDistFunc
		self.lDistFunc = lDistFunc
		self.nmaxmacroparticles	= nmaxmacroparticles
		self.injectturninterval = injectturninterval
    
	
	def addParticles(self):
		(xmin,xmax,ymin,ymax) = self.injectregion
	
		rank = 0
		numprocs = 1
		
		mpi_init = orbit_mpi.MPI_Initialized()
		comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
		
		if(mpi_init):
			rank = orbit_mpi.MPI_Comm_rank(comm)
			numprocs = orbit_mpi.MPI_Comm_size(comm)
		
		nPartsGlobal = self.bunch.getSizeGlobal()
		
		if(self.nmaxmacroparticles > 0):
			if(nPartsGlobal > self.nmaxmacroparticles):
				return
				
				#if((nTurnsDone % injectTurnInterval) != 0):
	#return
		
		x_rank0 = []
		xp_rank0 = []
		y_rank0 = []
		yp_rank0 = []
		z_rank0 = []
		dE_rank0 = []
		ninjected_rank0 = 0
		x_local = []
		xp_local = []
		y_local = []
		yp_local = []
		z_local = []
		dE_local = []
		ninjectedlocal  = 0

		if(rank == 0):
			for i in xrange(int(self.nparts)):
				(x,px) = self.xDistFunc.getCoordinates()
				(y,py) = self.yDistFunc.getCoordinates()
				(z,dE) = self.lDistFunc.getCoordinates()

				if((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)):
					ninjectedlocal = ninjectedlocal + 1
					x_rank0.append(x)
					xp_rank0.append(px)
					y_rank0.append(y)
					yp_rank0.append(py)
					z_rank0.append(z)
					dE_rank0.append(dE)
					#self.bunch.addParticle(x,px,y,py,z,dE)
				else:
					self.lostbunch.addParticle(x,px,y,py,z,dE)
					#self.bunch.compress()

		if(mpi_init):
			ninjected = orbit_mpi.MPI_Bcast(ninjectedlocal, mpi_datatype.MPI_INT, 0, comm)
			x_local = orbit_mpi.MPI_Bcast(x_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
			xp_local = orbit_mpi.MPI_Bcast(xp_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
			y_local = orbit_mpi.MPI_Bcast(y_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
			yp_local = orbit_mpi.MPI_Bcast(yp_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
			z_local = orbit_mpi.MPI_Bcast(z_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
			dE_local = orbit_mpi.MPI_Bcast(dE_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)

		n_remainder = ninjected % numprocs;
		n_inj_local = ninjected/numprocs;
	
		#inject the equal number of particles on each CPU
	
		i_start = rank * n_inj_local
		i_stop = (rank+1) * n_inj_local 
		for i in xrange(i_start, i_stop):
			self.bunch.addParticle(x_local[i],xp_local[i],y_local[i],yp_local[i],z_local[i],dE_local[i])
				
		n_max_index = numprocs * n_inj_local

		for i in xrange (n_remainder - 1):
			i_cpu = int( numprocs * random.random())
			orbit_mpi.MPI_Bcast(i_cpu, mpi_datatype.MPI_INT, 0,comm)
			if(rank == i_cpu):
				self.bunch.addParticle(x_local[i + 1 + n_max_index], xp_local[i + 1 + n_max_index], y_local[i + 1 + n_max_index],
									   yp_local[i + 1 + n_max_index],z_local[i + 1 + n_max_index], dE_local[i + 1 + n_max_index])
				n_inj_local = n_inj_local + 1
					
		#nInjectedMacros += n_inj_local;
		
		self.bunch.compress()
		self.lostbunch.compress()
			
	def addParticlesOld(self):
		(xmin,xmax,ymin,ymax) = self.injectregion
		
		nPartsGlobal = self.bunch.getSizeGlobal()
		
		if(self.nmaxmacroparticles > 0):
			if(nPartsGlobal > self.nmaxmacroparticles):
				return
		
		#if((nTurnsDone % self.injectturninterval) != 0):
				#	return

		for i in xrange(int(self.nparts)):
			(x,px) = self.xDistFunc.getCoordinates()
			(y,py) = self.yDistFunc.getCoordinates()
			(z,dE) = self.lDistFunc.getCoordinates()
			if((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)):
				self.bunch.addParticle(x,px,y,py,z,dE)
			else:
				self.lostbunch.addParticle(x,px,y,py,z,dE)
		
		self.bunch.compress()
		self.lostbunch.compress()
