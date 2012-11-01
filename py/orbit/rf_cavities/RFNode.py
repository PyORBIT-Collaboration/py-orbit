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
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync, RFHNum,\
			RFVoltage, RFPhase)
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

