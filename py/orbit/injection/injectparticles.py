#!/usr/bin/env python

"""
InjectParts parts class was modified on 2022.02.11 by A. Shishlo.
It was already parallel from the last modification.
The ParticleId and Initial coordinates were added as particle's attributes
to the new and lost particles. To assign these particles' attributes user 
should add them to the initial (empty) bunch before the injection starts.
"""

import math
import random
import sys

from bunch import Bunch

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

class InjectParts:
	""" 
	This routine injects particles into a bunch with user specified distribution
	functions.
	"""
	
	def __init__(self, nparts, bunch, lostbunch, injectregion, xDistFunc, yDistFunc, lDistFunc, nmaxmacroparticles = -1):
		self.nparts = nparts
		self.npartsfloat = float(nparts)
		self.bunch = bunch
		self.lostbunch = lostbunch
		self.injectregion = injectregion
		self.xDistFunc = xDistFunc
		self.yDistFunc = yDistFunc
		self.lDistFunc = lDistFunc
		self.nmaxmacroparticles	= nmaxmacroparticles
    
	
	def addParticles(self):
		"""
		Performs injections.
		"""
		#---- User can skip injection if self.nparts is zero
		if(self.nparts == 0): return
		random.seed(100)
		
		(xmin,xmax,ymin,ymax) = self.injectregion
	
		#adjusts number of particles injected according to varying pattern width
		if self.lDistFunc.name == "JohoLongitudinalPaint":
			self.lDistFunc.getCoordinates()
			self.npartsfloat = self.lDistFunc.frac_change*self.npartsfloat
			self.nparts = int(round(self.npartsfloat))
		
		rank = 0
		numprocs = 1
		
		mpi_init = orbit_mpi.MPI_Initialized()
		comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
		
		if(mpi_init):
			rank = orbit_mpi.MPI_Comm_rank(comm)
			numprocs = orbit_mpi.MPI_Comm_size(comm)
		
		nPartsGlobal = self.bunch.getSizeGlobal()
		
		if(self.nmaxmacroparticles > 0):
			if(nPartsGlobal >= self.nmaxmacroparticles):
				return
		
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
			for i in range(int(self.nparts)):
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
				else:
					self.addLostParticle(self.bunch,self.lostbunch,x,px,y,py,z,dE)

		nPartsLostGlobal = self.lostbunch.getSizeGlobal()
		nPartsTotalGlobal = nPartsGlobal + nPartsLostGlobal
		
		ninjected = orbit_mpi.MPI_Bcast(ninjectedlocal, mpi_datatype.MPI_INT, 0, comm)
		x_local = orbit_mpi.MPI_Bcast(x_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
		xp_local = orbit_mpi.MPI_Bcast(xp_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
		y_local = orbit_mpi.MPI_Bcast(y_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
		yp_local = orbit_mpi.MPI_Bcast(yp_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
		z_local = orbit_mpi.MPI_Bcast(z_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)
		dE_local = orbit_mpi.MPI_Bcast(dE_rank0, mpi_datatype.MPI_DOUBLE, 0, comm)

		n_remainder = ninjected % numprocs;
		n_inj_local = ninjected/numprocs;

		#---- inject the equal number of particles on each CPU
		i_start = rank * n_inj_local
		i_stop = (rank+1) * n_inj_local
		for i in range(i_start, i_stop):
			particleId = nPartsTotalGlobal + i
			self.addInjectedParticle(self.bunch,x_local[i],xp_local[i],y_local[i],yp_local[i],z_local[i],dE_local[i],particleId)
			
		#---- inject the reminder of the particles
		n_max_index = numprocs * n_inj_local

		nPartsGlobal = self.bunch.getSizeGlobal()
		nPartsLostGlobal = self.lostbunch.getSizeGlobal()
		nPartsTotalGlobal = nPartsGlobal + nPartsLostGlobal

		for i in range(n_remainder):
			i_cpu = random.randint(0,numprocs-1)
			i_cpu = orbit_mpi.MPI_Bcast(i_cpu, mpi_datatype.MPI_INT, 0,comm)
			if(rank == i_cpu):
				particleId = nPartsTotalGlobal + i
				self.addInjectedParticle(self.bunch,x_local[i + n_max_index], xp_local[i + n_max_index], 
					                                  y_local[i + n_max_index], yp_local[i + n_max_index], 
					                                  z_local[i + n_max_index], dE_local[i + n_max_index],particleId)
				#---- here n_inj_local is just for information for debugging
				n_inj_local = n_inj_local + 1
		
		self.bunch.compress()
		self.lostbunch.compress()

	def setnparts(self, nparts):
		"""
		Sets the number of partcicles that will be injected.
		"""
		self.nparts = nparts
		self.npartsfloat = float(nparts)

	def addInjectedParticle(self,bunch,x,px,y,py,z,dE,particleId):
		"""
		This method adds the injected particle to the bunch.
		It also will check if bunch has ParticleIdNumber attribute,
		in this case, an added particle will be assigned this Id.
		If bunch has InitailCoords attributes they will be assigned
		as well.
		"""
		bunch.addParticle(x,px,y,py,z,dE)
		if(bunch.hasPartAttr("ParticleIdNumber") != 0):
			particle_index = bunch.getSize() - 1
			bunch.partAttrValue("ParticleIdNumber", particle_index, 0, particleId)
		if(bunch.hasPartAttr("ParticleInitialCoordinates") != 0):
			particle_index = bunch.getSize() - 1
			bunch.partAttrValue("ParticleInitialCoordinates",particle_index, 0,x)
			bunch.partAttrValue("ParticleInitialCoordinates",particle_index, 1,px)
			bunch.partAttrValue("ParticleInitialCoordinates",particle_index, 2,y)
			bunch.partAttrValue("ParticleInitialCoordinates",particle_index, 3,py)
			bunch.partAttrValue("ParticleInitialCoordinates",particle_index, 4,z)
			bunch.partAttrValue("ParticleInitialCoordinates",particle_index, 5,dE)
	
	def addLostParticle(self,bunch,lostbunch,x,px,y,py,z,dE):
		"""
		This method adds the lost particle to the lost particles bunch.
		It also will check if lost bunch has ParticleIdNumber attribute,
		in this case, an added particle will be assigned Id = -1 .
		If bunch has InitailCoords attributes they will be assigned
		as well.
		"""
		#---- check bunch for particle attributes
		if(bunch.hasPartAttr("ParticleIdNumber") != 0 and 
			lostbunch.hasPartAttr("ParticleIdNumber") == 0):
			lostbunch.addPartAttr("ParticleIdNumber")
		if(bunch.hasPartAttr("ParticleInitialCoordinates") != 0 and 
			lostbunch.hasPartAttr("ParticleInitialCoordinates") == 0):
			lostbunch.addPartAttr("ParticleInitialCoordinates")	
		#----------------------------------------------------
		
		lostbunch.addParticle(x,px,y,py,z,dE)
		
		if(lostbunch.hasPartAttr("ParticleIdNumber") != 0):
			particle_index = lostbunch.getSize() - 1
			particleId = -1
			lostbunch.partAttrValue("ParticleIdNumber", particle_index, 0, particleId)
		
		if(lostbunch.hasPartAttr("ParticleInitialCoordinates") != 0):
			particle_index = lostbunch.getSize() - 1
			lostbunch.partAttrValue("ParticleInitialCoordinates",particle_index, 0,x)
			lostbunch.partAttrValue("ParticleInitialCoordinates",particle_index, 1,px)
			lostbunch.partAttrValue("ParticleInitialCoordinates",particle_index, 2,y)
			lostbunch.partAttrValue("ParticleInitialCoordinates",particle_index, 3,py)
			lostbunch.partAttrValue("ParticleInitialCoordinates",particle_index, 4,z)
			lostbunch.partAttrValue("ParticleInitialCoordinates",particle_index, 5,dE)	

