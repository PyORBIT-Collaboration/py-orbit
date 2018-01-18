"""
This package is a collection of diagnostics nodes fo linacs.
"""

import os
import math
import sys

#---- MPI module function and classes
import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

# import from orbit Python utilities
from orbit.utils import orbitFinalize
from orbit.utils import phaseNearTargetPhase, phaseNearTargetPhaseDeg
from orbit.utils import speed_of_light

# import from orbit c++ utilities
from orbit_utils import Polynomial
from orbit_utils import Function

# from LinacAccLattice import Sequence
from LinacAccLatticeLib import Sequence
from LinacAccNodes import BaseLinacNode, LinacNode

class LinacBPM(BaseLinacNode):
	"""
	The linac BPM representation. It calculates the average x,y, and phases
	of all particles in the bunch, and the amplitude signal which is an amplitude
	of the Fourier harmonics. To calculate this it needs the frequency [Hz] of BPM's
	electronics and the normalization factor.
	At this moment, all functionality is implemented on the Python level, but in
	the future the operations with the particles' coordinates should be moved to 
	C++ level.
	"""
	def __init__(self, frequency = 805.0e+6, name = "BPM"):
		BaseLinacNode.__init__(self,name)
		self.setType("BPM")
		self.frequency = 805.0e+6
		self.phase_hist_arr = []
		self.phase_step = 1.0 # in deg
		np = int(360./self.phase_step)
		self.phase_step = 360./np
		for ind in range(np):
			self.phase_hist_arr.append(0.)
		self.x_avg = 0.
		self.y_avg = 0.
		self.phase_avg = 0. # in deg
		self.phase_max = 0. # in deg
		self.synch_pahse = 0. # in deg
		self.rms_phase = 0. # in deg
		self.norm_coeff = 1.0
		self.fourier_amp = 0.
		self.fourier_phase = 0. # in deg
		self.amp = 0.
		
	def setNormalizationCoefficient(self,norm_coeff):
		self.norm_coeff = norm_coeff
	
	def getNormalizationCoefficient(self):
		return self.norm_coeff
		
	def setFrequency(self, frequency):
		self.frequency = frequency
		
	def getFrequency(self):
		return self.frequency
		
	def _cleanPhaseHist(self):
		for ind in range(len(self.phase_hist_arr)):
			self.phase_hist_arr[ind] = 0.
			
	def getAvgX(self):
		return self.x_avg
			
	def getAvgY(self):
		return self.y_avg
		
	def getPhaseRMS(self):
		return self.rms_phase
		
	def getAvgPhase(self):
		return self.phase_avg
		
	def getSynchPhase(self):
		return self.synch_pahse
		
	def getAmplitude(self):
		return self.amp
		
	def getFourierAmplitude(self):
		return self.fourier_amp
		
	def getFourierPhase(self):
		return self.fourier_phase
		
	def getPeakPhase(self):
		return self.phase_max
		
	def dumpPhaseHistorgam(self, file_out_name = None):
		nHistP = len(self.phase_hist_arr)
		fl_out = None
		if(file_out_name != None):
			fl_out = open(file_out_name,"w")
		for ind in range(nHistP):
			val = self.phase_hist_arr[ind]
			phase = (ind + 0.5)*self.phase_step - 180.
			st = " %3d "%ind
			st += " %7.2f  %12.5g "%(phase,val)
			if(fl_out != None):
				fl_out.write(st+"\n")
			else:
				print st
		if(fl_out != None): fl_out.close()
		
	def track(self, paramsDict):
		"""
		The LinacBPM class implementation of the AccNode class track(paramsDict) method.
		"""
		bunch = paramsDict["bunch"]
		self.analyzeBunch(bunch)

	def analyzeBunch(self,bunch):
		"""
		Here we assume that the macrosize is the same for each 
		particle
		"""
		comm = bunch.getMPIComm()	
		self._cleanPhaseHist()
		self.x_avg = 0.
		self.y_avg = 0.
		self.phase_avg = 0. # in deg
		self.phase_max = 0. # in deg
		self.synch_pahse = 0. # in d
		self.fourier_amp = 0.
		self.fourier_phase = 0. # in deg
		self.amp = 0.
		nPartsGlobal = bunch.getSizeGlobal()
		if( nPartsGlobal == 0): return		
		x_avg = 0.
		y_avg = 0.
		phase_avg = 0.
		beta = bunch.getSyncParticle().beta()
		z_to_phase = - 360.*self.frequency/(speed_of_light*beta)
		synch_phase = 360.*bunch.getSyncParticle().time()*self.frequency
		self.synch_pahse = phaseNearTargetPhaseDeg(synch_phase,0.)
		nParts = bunch.getSize()
		for ind in range(nParts):
			x_avg += bunch.x(ind)
			y_avg += bunch.y(ind)
			phase_avg += z_to_phase*bunch.z(ind)
		#---- for parallel case
		(x_avg,y_avg,phase_avg) = orbit_mpi.MPI_Allreduce((x_avg,y_avg,phase_avg),mpi_datatype.MPI_DOUBLE,mpi_op.MPI_SUM,comm)
		x_avg /= nPartsGlobal
		y_avg /= nPartsGlobal
		phase_avg /= nPartsGlobal
		phase2_avg = 0.
		for ind in range(nParts):
			phase2_avg += (z_to_phase*bunch.z(ind)-phase_avg)**2
		phase2_avg /= nPartsGlobal
		self.rms_phase = math.sqrt(phase2_avg)
		self.x_avg = x_avg
		self.y_avg = y_avg
		phase_avg += synch_phase		
		self.phase_avg = phaseNearTargetPhaseDeg(phase_avg,0.)
		#------- fill out histogram
		nHistP = len(self.phase_hist_arr)
		for ind in range(nParts):
			phase_ind = int((180. + phaseNearTargetPhaseDeg(z_to_phase*bunch.z(ind),0.))/self.phase_step)
			if(phase_ind < 0): phase_ind = 0
			if(phase_ind >= nHistP): phase_ind = nHistP - 1
			self.phase_hist_arr[phase_ind] += 1.0
		phase_hist_arr = orbit_mpi.MPI_Allreduce(self.phase_hist_arr,mpi_datatype.MPI_DOUBLE,mpi_op.MPI_SUM,comm)
		phase_hist_arr = list(phase_hist_arr)
		#---- find the position of the max value
		total_sum = 0.
		ind_max = -1
		max_val = 0
		for ind in range(nHistP):
			val = phase_hist_arr[ind]
			if(val > max_val):
				ind_max = ind
				max_val = val
			self.phase_hist_arr[ind] = val
			total_sum += val
		self.phase_max = (ind_max + 0.5)*self.phase_step
		self.phase_max -= 180. 
		self.phase_max += synch_phase
		self.phase_max = phaseNearTargetPhaseDeg(self.phase_max,0.)
		#---- calculate Fourier amplitude
		sin_sum = 0.
		cos_sum = 0.
		grad_to_rad_coeff = math.pi/180.
		phase_step = self.phase_step*grad_to_rad_coeff
		#---- normalization
		n_local_coeff = 1./(total_sum*phase_step)
		for ind in range(nHistP):
			phase_hist_arr[ind] *= 	n_local_coeff
		for ind in range(nHistP):
			self.phase_hist_arr[ind] = phase_hist_arr[ind]
		#--- Fourier amplitude and phase
		for ind in range(nHistP):
			val = phase_hist_arr[ind]
			phase = (ind + 0.5)*self.phase_step*grad_to_rad_coeff
			sin_sum += val*math.sin(phase)
			cos_sum += val*math.cos(phase)
		sin_sum *= self.phase_step*grad_to_rad_coeff
		cos_sum *= self.phase_step*grad_to_rad_coeff
		self.fourier_amp = math.sqrt(sin_sum**2 + cos_sum**2)/math.pi
		self.amp = self.norm_coeff*self.fourier_amp
		self.fourier_phase = math.atan2(cos_sum,sin_sum)*180./math.pi + 90.
		self.fourier_phase += synch_phase
		self.fourier_phase = phaseNearTargetPhaseDeg(self.fourier_phase,0.)
		
			
			
		