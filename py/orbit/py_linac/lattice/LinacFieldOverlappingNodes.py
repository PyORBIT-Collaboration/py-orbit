#!/usr/bin/env python

#--------------------------------------------------------
# This is a collection of classes for describing the nodes
# with the overlapping quad or rf+quad fields.
#--------------------------------------------------------

import math
import sys
import os

from orbit.utils import orbitFinalize, phaseNearTargetPhase, phaseNearTargetPhaseDeg

#---- base linac nodes
from LinacAccNodes import AbstractRF_Gap
from LinacAccNodes import BaseLinacNode, Drift, Quad

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# from linac import the RF gap classes
from linac import RfGapThreePointTTF
from linac import RfGapThreePointTTF_slow

class AxisField_and_Quad_RF_Gap(AbstractRF_Gap):
	"""
	The class represents the part of the RF gap. It uses the part (longitudinally) 
	of the RF axis field that can be overlapped by quadrupole magnet field.
	It is almost copy of the AxisFieldRF_Gap class, but it describes the part of 
	the whole RF gap because the axis field could cover the dipole correctors or
	BPMs, and it has to be cut in peaces. It also could overlap the neighboring 
	quadrupoles, therefore it should account for the tracking inside quad fields. 
	
	The method setAsFirstRFGap(True/False) will be called by RF cavity.
	
	User have to provide the instance of the AxisFieldRF_Gap at the constructor.
	The instance of this class has the axis field for the whole gap.
	"""

	def __init__(self, axis_field_rf_gap):
		"""
		Constructor for the axis field RF gap. 
		The axis_field_rf_gap is the instance of the AbstractRF_Gap class.
		E0L parameter is in GeV. Phases are in radians.
		"""
		AbstractRF_Gap.__init__(self,axis_field_rf_gap.getName())
		self.setType("GAP&Q")
		self.axis_field_rf_gap = axis_field_rf_gap	
		self.addParam("E0TL",self.axis_field_rf_gap.getParam("E0TL"))
		self.addParam("mode",self.axis_field_rf_gap.getParam("mode"))
		self.addParam("gap_phase",self.axis_field_rf_gap.getParam("gap_phase"))
		self.addParam("rfCavity",self.axis_field_rf_gap.getParam("rfCavity"))
		self.addParam("E0L",self.axis_field_rf_gap.getParam("E0L"))
		self.addParam("EzFile",self.axis_field_rf_gap.getParam("EzFile"))
		self.setPosition(self.axis_field_rf_gap.getPosition())
		#---- aperture parameters
		if(axis_field_rf_gap.hasParam("aperture") and axis_field_rf_gap.hasParam("aprt_type")):
			self.addParam("aperture",axis_field_rf_gap.getParam("aperture"))
			self.addParam("aprt_type",axis_field_rf_gap.getParam("aprt_type"))
		#---- axis field related parameters
		self.z_step = 0.01
		self.z_min = 0.
		self.z_max = 0.
		#---- gap_phase_vs_z_arr keeps [pos,phase] pairs after the tracking
		self.gap_phase_vs_z_arr = []
		#---- The position of the particle during the run. 
		#---- It is used for the path length accounting.
		self.part_pos = 0.
		#---- The RF gap model - three points model
		self.cppGapModel = RfGapThreePointTTF()
		#---- If we going to use the longitudinal magnetic field component of quad
		self.useLongField = False
		#---- quadrupole field sources
		#----quads_fields_arr is an array of [quad, fieldFunc, z_center_of_field]
		self.quads_fields_arr = []
		#---- If it is true then the this tracking will be in the reversed lattice
		self.reversed_lattice = False		
		
	def setLinacTracker(self, switch = True):
		"""
		This method will switch RF gap model to slower one where transformations 
		coefficients are calculated for each particle in the bunch.
		"""
		BaseLinacNode.setLinacTracker(self,switch)
		AbstractRF_Gap.setLinacTracker(self,switch)
		if(switch):
			self.cppGapModel = RfGapThreePointTTF_slow()			
		else:
			self.cppGapModel = RfGapThreePointTTF()

	def setUseLongitudinalFieldOfQuad(self, use):
		"""
		If we going to use the longitudinal magnetic field component of quad
		"""
		self.useLongField = use
		
	def getUseLongitudinalFieldOfQuad(self):
		"""
		If we going to use the longitudinal magnetic field component of quad
		"""
		return self.useLongField

	def getAxisFieldRF_Gap(self):
		"""
		It returns the  AxisFieldRF_Gap instance for this gap.
		"""
		return self.axis_field_rf_gap
		
	def getBaseRF_Gap(self):
		"""
		It returns the BaseRF_Gap instance for this gap.
		"""		
		return self.axis_field_rf_gap.baserf_gap
		
	def reverseOrderNodeSpecific(self):
		"""
		This method is used for a lattice reversal and a bunch backtracking
		This is a node type specific method. The implementation of the abstract
		method of AccNode class from the top level lattice package. 
		"""
		self.reversed_lattice = not self.reversed_lattice
		#---- Here the order of quads does not matter
		for quad_ind in range(len(self.quads_fields_arr)):
			[quad, fieldFunc, z_center_of_field] = self.quads_fields_arr[quad_ind]
			self.quads_fields_arr[quad_ind][2] = - z_center_of_field
		(self.z_min,self.z_max) = (-self.z_max,-self.z_min)  
			
	def isNodeInReversedLattice(self):
		"""
		Returns True if this node in the reversed lattice or False if it otherwise.
		"""
		return self.reversed_lattice		
		
	def addQuad(self, quad, fieldFunc, z_center_of_field):
		"""
		Adds the quad with the field function and the position.
		The position of the quad is relative to the center of 
		the parent rf gap node.
		"""
		self.quads_fields_arr.append([quad, fieldFunc, z_center_of_field])		
		
	def getQuads(self):
		"""
		Returns the list of quads in this node.
		"""
		quads = []
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			quads.append(quad)
		return quads		
		
	def getPosAndQuad_Arr(self):
		"""
		Return the array with pairs: the quad and the position of its center.
		"""
		quad_arr = []
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			quad_arr.append([quad,z_center_of_field])
		return quad_arr

	def getTotalField(self,z_from_center):
		"""
		Returns the combined field of all overlapping quads.
		z_from_center - is a distance from the center of the parent RF gap node.
		"""
		z = z_from_center
		G = 0.
		if(z < self.z_min or z > self.z_max): return G
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			if(fieldFunc != None):
				gl = quad.getParam("dB/dr")*quad.getLength()
				G += gl*fieldFunc.getFuncValue(z - z_center_of_field)
			else:
				G += quad.getTotalField(z - z_center_of_field)
		return G		

	def getTotalFieldDerivative(self,z_from_center):
		"""
		Returns the combined derivative of the field of all overlapping quads.
		z_from_center - is a distance from the center of the parent RF gap node.
		"""
		z = z_from_center
		GP = 0.
		if(z < self.z_min or z > self.z_max): return GP
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			if(fieldFunc != None):
				gl = quad.getParam("dB/dr")*quad.getLength()
				GP += gl*fieldFunc.getFuncDerivative(z - z_center_of_field)
			else:
				GP += 0.
		return GP

	def getZ_Step(self):
		"""
		Returns the longitudinal step during the tracking.
		"""
		return self.z_step
		
	def setZ_Step(self,z_step):
		if(self.axis_field_rf_gap.axis_field_func == None):
			msg  = "Class AxisFieldRF_Gap: You have to get the axis field from a file first!"
			msg += os.linesep
			msg += "Call readAxisFieldFile(dir_location,file_name) method first!"
			msg += "Stop."
			orbitFinalize(msg)				
		length = self.getLength()
		nParts = int(length*1.0000001/z_step)
		if(nParts < 1): nParts = 1
		self.z_step = length/nParts
		#---- this will set the even distribution of the lengths between parts
		self.setnParts(nParts)

	def getZ_Min_Max(self):
		"""
		Returns the tuple (z_min,z_max) with the limits of the axis field.
		These parameters define the length of the node. The center of the node
		is at 0.
		"""
		return (self.z_min,self.z_max)

	def setZ_Min_Max(self,z_min,z_max):
		"""
		Sets the actual longitudinal sizes of the node. It is used for small correction
		of the length to avoid fields overlapping from neighbouring gaps.
		"""
		self.z_min = z_min
		self.z_max = z_max
		length = self.z_max - self.z_min
		self.setLength(length)
		self.setZ_Step(self.z_step)

	def getEzFiled(self,z):
		"""
		Returns the Ez field on the axis of the RF gap in V/m. 
		"""
		rfCavity = self.getRF_Cavity()
		E0L = 1.0e+9*self.getParam("E0L")
		rf_ampl = rfCavity.getAmp()
		Ez = self.getEzFiledInternal(z,rfCavity,E0L,rf_ampl)
		return Ez

	def getEzFiledInternal(self,z,rfCavity,E0L,rf_ampl):
		"""
		Returns the Ez field on the axis of the RF gap in V/m. 
		"""
		if(self.reversed_lattice): z *= -1
		Ez = E0L*rf_ampl*self.axis_field_rf_gap.axis_field_func.getY(z)	
		return Ez

	def getRF_Cavity(self):
		"""
		Returns the parent RF Cavity.
		"""
		return self.getParam("rfCavity")
				
	def track(self, paramsDict):
		"""
		The AxisFieldRF_Gap class implementation of
		the AccNode class track(probe) method.
		User have to track the design bunch first to setup all gaps arrival time. 
		"""
		rfCavity = self.getRF_Cavity()
		if(not rfCavity.isDesignSetUp()):
			sequence = self.getSequence()
			accLattice = sequence.getLinacAccLattice()
			msg  = "The AxisFieldRF_Gap class. "
			msg += "You have to run trackDesign on the LinacAccLattice"
			msg += "first to initialize all RF Cavities' phases!"
			msg += os.linesep
			if(accLattice != None):				
				msg = msg + "Lattice =" + accLattice.getName()				
				msg = msg + os.linesep
			if(sequence != None):				
				msg = msg + "Sequence =" + sequence.getName()				
				msg = msg + os.linesep
			msg = msg + "RF Cavity =" + rfCavity.getName()				
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			orbitFinalize(msg)
		#-----------------------------------------
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		part_length = self.getLength(index)		
		bunch = paramsDict["bunch"]
		charge = bunch.charge()
		syncPart = bunch.getSyncParticle()	
		eKin_in = syncPart.kinEnergy()
		momentum = syncPart.momentum()
		E0L = 1.0e+9*self.getParam("E0L")
		modePhase = self.axis_field_rf_gap.baserf_gap.getParam("mode")*math.pi
		frequency = rfCavity.getFrequency()	
		rf_ampl = rfCavity.getAmp()
		arrival_time = syncPart.time()
		designArrivalTime = rfCavity.getDesignArrivalTime()
		phase_shift = rfCavity.getPhase() - rfCavity.getDesignPhase()
		phase = rfCavity.getFirstGapEtnrancePhase() + phase_shift
		#----------------------------------------
		phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi + phase,2.0*math.pi)		
		if(index == 0):
			self.part_pos = self.z_min 
			self.gap_phase_vs_z_arr = [[self.part_pos,phase],]
		zm = self.part_pos
		z0 = zm + part_length/2
		zp = z0 + part_length/2
		Em = self.getEzFiledInternal(zm,rfCavity,E0L,rf_ampl)
		E0 = self.getEzFiledInternal(z0,rfCavity,E0L,rf_ampl)
		Ep = self.getEzFiledInternal(zp,rfCavity,E0L,rf_ampl)			
		#------- track through a quad
		G = self.getTotalField((zm+z0)/2)
		dB_dz = 0.
		if(self.useLongField == True): dB_dz = self.getTotalFieldDerivative((zm+z0)/2)
		if(abs(G) != 0.):
			kq = G/bunch.B_Rho()
			#------- track through a quad
			step = part_length/2
			self.tracking_module.quad1(bunch,step/4.0, kq)
			self.tracking_module.quad2(bunch,step/2.0)
			self.tracking_module.quad1(bunch,step/2.0, kq)
			self.tracking_module.quad2(bunch,step/2.0)
			self.tracking_module.quad1(bunch,step/4.0, kq)
			if(abs(dB_dz) != 0.):
				self.tracking_module.quad3(bunch,step,dB_dz)
		else:
			self.tracking_module.drift(bunch,part_length/2)
		self.part_pos += part_length/2
		#call rf gap model to track the bunch
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)		
		#print s
		#---- this part is the debugging ---STOP---
		self.cppGapModel.trackBunch(bunch,part_length/2,Em,E0,Ep,frequency,phase+delta_phase+modePhase)
		#------- track through a quad		
		G = self.getTotalField((z0+zp)/2)
		dB_dz = 0.
		if(self.useLongField == True): dB_dz = self.getTotalFieldDerivative((z0+zp)/2)		
		if(abs(G) != 0.):
			kq = G/bunch.B_Rho()
			step = part_length/2
			self.tracking_module.quad1(bunch,step/4.0, kq)
			self.tracking_module.quad2(bunch,step/2.0)
			self.tracking_module.quad1(bunch,step/2.0, kq)
			self.tracking_module.quad2(bunch,step/2.0)
			self.tracking_module.quad1(bunch,step/4.0, kq)
			if(abs(dB_dz) != 0.):
				self.tracking_module.quad3(bunch,step,dB_dz)
		else:
			self.tracking_module.drift(bunch,part_length/2)
		#---- advance the particle position
		self.part_pos += part_length/2
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)
		#print s
		#---- this part is the debugging ---STOP---
		#---- Calculate the phase at the center
		if(index == (nParts - 1)):
			pos_old = self.gap_phase_vs_z_arr[0][0]			
			phase_gap = self.gap_phase_vs_z_arr[0][1]
			ind_min = -1
			for ind in range(1,len(self.gap_phase_vs_z_arr)):
				[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
				if(math.fabs(pos) >= math.fabs(pos_old)):
					ind_min = ind -1
					phase_gap = self.gap_phase_vs_z_arr[ind_min][1]
					phase_gap = phaseNearTargetPhase(phase_gap,0.)
					self.gap_phase_vs_z_arr[ind_min][1] = phase_gap
					break
				pos_old = pos
			self.setGapPhase(phase_gap)
			#---- wrap all gap part's phases around the central one
			if(ind_min > 0):
				for ind in range(ind_min-1,-1,-1):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind+1]
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)
				for ind in range(ind_min+1,len(self.gap_phase_vs_z_arr)):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind-1]				
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)

	def trackDesign(self, paramsDict):
		"""
		The method is tracking the design synchronous particle through the RF Gap.
		If the gap is a first gap in the cavity we put the arrival time as 
		a cavity parameter. The pair of the cavity design phase and this arrival time 
		at the first gap are used during the real bunch tracking.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		part_length = self.getLength(index)
		bunch = paramsDict["bunch"]
		syncPart = bunch.getSyncParticle()
		eKin_in = syncPart.kinEnergy()
		#---- parameter E0L is in GeV, but cppGapModel = RfGapThreePointTTF() uses fields in V/m
		E0L = 1.0e+9*self.getParam("E0L")
		modePhase = self.axis_field_rf_gap.baserf_gap.getParam("mode")*math.pi
		rfCavity = self.getRF_Cavity()
		rf_ampl = rfCavity.getDesignAmp()
		arrival_time = syncPart.time()
		frequency = rfCavity.getFrequency()
		phase = rfCavity.getFirstGapEtnrancePhase()
		#---- calculate the entance phase
		if(self.isFirstRFGap() and index == 0):
			rfCavity.setDesignArrivalTime(arrival_time)
			phase = self.axis_field_rf_gap.calculate_first_part_phase(bunch)
			rfCavity.setFirstGapEtnrancePhase(phase)
			rfCavity.setFirstGapEtnranceDesignPhase(phase)
			rfCavity.setDesignSetUp(True)
			rfCavity._setDesignPhase(rfCavity.getPhase())
			rfCavity._setDesignAmp(rfCavity.getAmp())
			#print "debug firs gap first part phase=",phase*180./math.pi," arr time=",arrival_time
		else:
			first_gap_arr_time = rfCavity.getDesignArrivalTime()
			#print "debug name=",self.getName()," delta_phase=",frequency*(arrival_time - first_gap_arr_time)*360.0," phase=",phase*180/math.pi
			phase = math.fmod(frequency*(arrival_time - first_gap_arr_time)*2.0*math.pi+phase,2.0*math.pi)		
		#print "debug design name=",self.getName()," arr_time=",arrival_time," phase=",phase*180./math.pi," E0TL=",E0TL*1.0e+3," freq=",frequency
		if(index == 0):
			self.part_pos = self.z_min 
			self.gap_phase_vs_z_arr = [[self.part_pos,phase],]
		zm = self.part_pos
		z0 = zm + part_length/2
		zp = z0 + part_length/2
		Em = self.getEzFiledInternal(zm,rfCavity,E0L,rf_ampl)
		E0 = self.getEzFiledInternal(z0,rfCavity,E0L,rf_ampl)
		Ep = self.getEzFiledInternal(zp,rfCavity,E0L,rf_ampl)
		#---- advance the particle position
		self.tracking_module.drift(bunch,part_length/2)
		self.part_pos += part_length/2	
		#call rf gap model to track the bunch
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)		
		#print s
		#---- this part is the debugging ---STOP---
		self.cppGapModel.trackBunch(bunch,part_length/2,Em,E0,Ep,frequency,phase+delta_phase+modePhase)
		self.tracking_module.drift(bunch,part_length/2)
		#---- advance the particle position
		self.part_pos += part_length/2
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)
		#print s
		#---- this part is the debugging ---STOP---
		#---- Calculate the phase at the center
		if(index == (nParts - 1)):
			pos_old = self.gap_phase_vs_z_arr[0][0]			
			phase_gap = self.gap_phase_vs_z_arr[0][1]
			ind_min = -1
			for ind in range(1,len(self.gap_phase_vs_z_arr)):
				[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
				if(math.fabs(pos) >= math.fabs(pos_old)):
					ind_min = ind -1
					phase_gap = self.gap_phase_vs_z_arr[ind_min][1]
					phase_gap = phaseNearTargetPhase(phase_gap,0.)
					self.gap_phase_vs_z_arr[ind_min][1] = phase_gap
					break
				pos_old = pos
			self.setGapPhase(phase_gap)
			#---- wrap all gap part's phases around the central one
			if(ind_min > 0):
				for ind in range(ind_min-1,-1,-1):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind+1]
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)
				for ind in range(ind_min+1,len(self.gap_phase_vs_z_arr)):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind-1]				
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)
		

class OverlappingQuadsNode(BaseLinacNode):
	"""
	The class represent the set of quads with the overlapping fields.
	"""
	def __init__(self, name = "OverlappingQuads"):
		"""
		Constructor. Creates the OverlappingQuadsNode instance.
		"""
		BaseLinacNode.__init__(self,name = "OVRLPQ")
		self.setType("OVRLPQ")
		self.setnParts(1)
		#----quads_fields_arr is an array of [quad, fieldFunc, z_center_of_field]
		self.quads_fields_arr = []
		#---- current z position
		self.z_value = 0.
		#---- z-step - the step in longitudinal direction during the tracking in [m]
		self.z_step = 0.01 
		#---- the initial length is 0.
		self.setLength(0.)
		#---- If we going to use the longitudinal magnetic field component of quad
		self.useLongField = False
		#---- If it is true then the this tracking will be in the reversed lattice
		self.reversed_lattice = False
		
	def setUseLongitudinalFieldOfQuad(self, use):
		"""
		If we going to use the longitudinal magnetic field component of quad
		"""
		self.useLongField = use
		
	def getUseLongitudinalFieldOfQuad(self):
		"""
		If we going to use the longitudinal magnetic field component of quad
		"""
		return self.useLongField		
		
	def reverseOrderNodeSpecific(self):
		"""
		This method is used for a lattice reversal and a bunch backtracking
		This is a node type specific method. The implementation of the abstract
		method of AccNode class from the top level lattice package.  
		"""
		self.reversed_lattice = not self.reversed_lattice
		#---- Here the order of quads does not matter
		for quad_ind in range(len(self.quads_fields_arr)):
			[quad, fieldFunc, z_center_of_field] = self.quads_fields_arr[quad_ind]
			self.quads_fields_arr[quad_ind][2] = - z_center_of_field + self.getLength()
			
	def isNodeInReversedLattice(self):
		"""
		Returns True if this node in the reversed lattice or False if it otherwise.
		"""
		return self.reversed_lattice
		
	def addQuad(self, quad, fieldFunc, z_center_of_field):
		"""
		Adds the quad with the field function and the position.
		The position of the quad is relative to the beginning of this OverlappingQuadsNode.
		"""
		self.quads_fields_arr.append([quad, fieldFunc, z_center_of_field])
				
	def setZ_Step(self,z_step):
		"""
		Sets the longitudinal step for the tracking along the node.
		"""
		self.z_step = z_step
		
	def getZ_Step(self):
		"""
		Returns the longitudinal step for the tracking along the node.
		"""		
		return self.z_step
		
	def getZ_Min_Max(self):
		"""
		Returns the tuple (z_min,z_max) with the limits of z coordinate from the center.
		These parameters define the length of the node. The center of the node
		is at 0.
		"""
		L2 = self.getLength()/2
		return (-L2,L2)	
		
	def getQuads(self):
		"""
		Returns the list of quads in this node.
		"""
		quads = []
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			quads.append(quad)
		return quads
	
	def getCentersOfField(self):
		"""
		Returns the array of centers of the quads in this node.
		"""
		centers_arr = []
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			centers_arr.append(z_center_of_field)
		return centers_arr

	def initialize(self):
		"""
		The OverlappingQuadsNode class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		length = self.getLength()
		lengthStep = length/nParts
		for i in xrange(nParts):
			self.setLength(lengthStep,i)

	def track(self, paramsDict):
		"""
		The  OverlappingQuadsNode class implementation of the AccNode class track(paramDict) method.
		"""
		index = self.getActivePartIndex()	
		length = self.getLength(index)
		if(index == 0): self.z_value = - self.getLength()/2
		bunch = paramsDict["bunch"]
		charge = bunch.charge()
		momentum = bunch.getSyncParticle().momentum()		
		n_steps = int(length/self.z_step)+1
		z_step = length/n_steps
		for z_ind in range(n_steps):
			z = self.z_value + z_step*(z_ind+0.5)
			G = self.getTotalField(z)
			dB_dz = 0.
			if(self.useLongField == True): dB_dz = self.getTotalFieldDerivative(z)			
			kq = G/bunch.B_Rho()
			if(abs(kq) == 0.):
				self.tracking_module.drift(bunch,z_step)
				continue
			#------- track through a combined quad
			self.tracking_module.quad1(bunch,z_step/4.0, kq)
			self.tracking_module.quad2(bunch,z_step/2.0)
			self.tracking_module.quad1(bunch,z_step/2.0, kq)
			self.tracking_module.quad2(bunch,z_step/2.0)
			self.tracking_module.quad1(bunch,z_step/4.0, kq)
			if(abs(dB_dz) != 0.):
				self.tracking_module.quad3(bunch,z_step, dB_dz)			
		self.z_value += length
		
	def getTotalField(self,z_from_center):
		"""
		Returns the combined field of all overlapping quads.
		z_from_center - is a distance from the center of the node.
		z - is a distance from the beginning of the node.
		"""
		z = z_from_center + self.getLength()/2
		G = 0.
		if(z < 0. or z > self.getLength()): return G
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			if(fieldFunc != None):
				gl = quad.getParam("dB/dr")*quad.getLength()
				G += gl*fieldFunc.getFuncValue(z - z_center_of_field)
			else:
				G += quad.getTotalField(z - z_center_of_field)
		return G
		
	def getTotalFieldDerivative(self,z_from_center):
		"""
		Returns the combined derivative of the field of all overlapping quads.
		z_from_center - is a distance from the center of the node.
		z - is a distance from the beginning of the node.
		"""
		z = z_from_center + self.getLength()/2
		GP = 0.
		if(z < 0. or z > self.getLength()): return GP
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			if(fieldFunc != None):
				gl = quad.getParam("dB/dr")*quad.getLength()
				GP += gl*fieldFunc.getFuncDerivative(z - z_center_of_field)
			else:
				GP += 0.
		return GP	

