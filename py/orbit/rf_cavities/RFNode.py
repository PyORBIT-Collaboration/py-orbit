"""
Module. Includes classes for RF accelerator nodes.
"""

import sys
import os
import math

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

#import RF cavity classes
from rfcavities import Frequency_Cav
from rfcavities import Harmonic_Cav

class Base_RFNode(DriftTEAPOT):

	def __init__(self, length, name = "base_rfnode"):
		"""
			Constructor. Creates Base RF Cavity TEAPOT element.
			It will never be called.
		"""
		DriftTEAPOT.__init__(self, name)
		self.setType("base rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		#self..trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		#self..trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class Frequency_RFNode(Base_RFNode):

	def __init__(self, RFFreq, RFE0TL, RFPhase,\
		length, name = "frequency_rfnode"):
		"""
			Constructor. Creates Frequency
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.frequencynode = Frequency_Cav(RFFreq, RFE0TL, RFPhase)
		self.setType("frequency rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.frequencynode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.frequencynode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase,\
		length, name = "harmonic_rfnode"):
		"""
			Constructor. Creates Harmonic
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			RFHNum, RFVoltage, RFPhase)
		self.setType("harmonic rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class BRhoDep_Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, accelDict, bunch,\
		length, name = "brho_timedep_harmonic_rfnode"):
		"""
			Constructor. Creates BRho Time Dependent
			Harmonic RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.Z2Phi = ZtoPhi
		self.localDict = accelDict
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		BRho_tuple = self.localDict["BRho"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		mass = bunch.mass()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		eold = keold + mass
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		BRho = interp(time, n_tuple,\
			time_tuple, BRho_tuple)
		pcnew = 0.299792458 * charge * BRho
		enew = math.sqrt(pcnew * pcnew + mass * mass)
		kenew = enew - mass
		bunch.getSyncParticle().kinEnergy(kenew)
		dESync = enew - eold
		Zsync = syncZ(ZtoPhi, gammaTrans, gamma, charge,\
			dESync, RFHNum, RFVoltage, RFPhase)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			RFHNum, RFVoltage, RFPhase)
		self.setType("harmonic rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		BRho_tuple = self.localDict["BRho"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		mass = bunch.mass()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		eold = keold + mass
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		BRho = interp(time, n_tuple,\
			time_tuple, BRho_tuple)
		pcnew = 0.299792458 * charge * BRho
		enew = math.sqrt(pcnew * pcnew + mass * mass)
		kenew = enew - mass
		bunch.getSyncParticle().kinEnergy(kenew)
		dESync = enew - eold
		ZtoPhi = self.Z2Phi
		Zsync = syncZ(ZtoPhi, gammaTrans, gamma, charge,\
			dESync, RFHNum, RFVoltage, RFPhase)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		BRho_tuple = self.localDict["BRho"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		mass = bunch.mass()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		eold = keold + mass
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		BRho = interp(time, n_tuple,\
			time_tuple, BRho_tuple)
		pcnew = 0.299792458 * charge * BRho
		enew = math.sqrt(pcnew * pcnew + mass * mass)
		kenew = enew - mass
		bunch.getSyncParticle().kinEnergy(kenew)
		dESync = enew - eold
		ZtoPhi = self.Z2Phi
		Zsync = syncZ(ZtoPhi, gammaTrans, gamma, charge,\
			dESync, RFHNum, RFVoltage, RFPhase)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class SyncPhaseDep_Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, accelDict, bunch,\
		length, name = "syncphase_timedep_harmonic_rfnode"):
		"""
			Constructor. Creates SyncPhase Time Dependent
			Harmonic RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.Z2Phi = ZtoPhi
		self.localDict = accelDict
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			RFHNum, RFVoltage, RFPhase)
		self.setType("harmonic rf node")
		self.setLength(0.0)
		
		print "Params1 = ", time, dESync, keold, kenew
		print "Params2 = ", RFVoltage, RFPhase, SyncPhase, Zsync

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		ZtoPhi = self.Z2Phi
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length
		
		print "Params1 = ", time, dESync, keold, kenew
		print "Params2 = ", RFVoltage, RFPhase, SyncPhase, Zsync

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		ZtoPhi = self.Z2Phi
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

		print "Params1 = ", time, dESync, keold, kenew
		print "Params2 = ", RFVoltage, RFPhase, SyncPhase, Zsync

class Barrier_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, RFVoltage, RFPhasep,\
		RFPhasem, dRFPhasep, dRFPhasem,\
		length, name = "barrier_rfnode"):
		"""
			Constructor. Creates Barrier
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.barriernode = Barrier_Cav(ZtoPhi, RFVoltage,\
			RFPhasep, RFPhasem, dRFPhasep, dRFPhasem)
		self.setType("barrier rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.barriernode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.barriernode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class SyncPhaseDep_Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, accelDict, bunch,\
		length, name = "syncphase_timedep_harmonic_rfnode"):
		"""
			Constructor. Creates SyncPhase Time Dependent
			Harmonic RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.Z2Phi = ZtoPhi
		self.localDict = accelDict
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			RFHNum, RFVoltage, RFPhase)
		self.setType("harmonic rf node")
		self.setLength(0.0)
		
		print "Params1 = ", time, dESync, keold, kenew
		print "Params2 = ", RFVoltage, RFPhase, SyncPhase, Zsync

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		ZtoPhi = self.Z2Phi
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length
		
		print "Params1 = ", time, dESync, keold, kenew
		print "Params2 = ", RFVoltage, RFPhase, SyncPhase, Zsync

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		ZtoPhi = self.Z2Phi
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

		print "Params1 = ", time, dESync, keold, kenew
		print "Params2 = ", RFVoltage, RFPhase, SyncPhase, Zsync

def interp(x, n_tuple, x_tuple, y_tuple):
	"""
	Linear interpolation: Given n-tuple + 1 points,
	x_tuple and y_tuple, routine finds y = y_tuple
	at x in x_tuple. Assumes x_tuple is increasing array.
	"""
	if x <= x_tuple[0]:
		y = y_tuple[0]
		return y
	if x >= x_tuple[n_tuple]:
		y = y_tuple[n_tuple]
		return y
	dxp = x - x_tuple[0]
	for n in range(n_tuple):
		dxm = dxp
		dxp = x - x_tuple[n + 1]
		dxmp = dxm * dxp
		if dxmp <= 0:
			break
	y = (-dxp * y_tuple[n] + dxm * y_tuple[n + 1]) /\
		(dxm - dxp)
	return y

def syncZ(ZtoPhi, gammaTrans, gamma, charge,\
	dESync, RFHNum, RFVoltage, RFPhase):
	"""
	Calculates position of synchronous particle.
	"""
	Zsync = 0
	if abs(dESync) > abs(charge * RFVoltage):
		return Zsync
	PhaseTot = math.asin(dESync / (charge * RFVoltage))
	if gamma > gammaTrans and gammaTrans > 0:
		PhaseTot = math.pi - PhaseTot
	Zsync = -(PhaseTot - math.pi * RFPhase / 180.0)\
		/ (RFHNum * ZtoPhi)
	return Zsync

