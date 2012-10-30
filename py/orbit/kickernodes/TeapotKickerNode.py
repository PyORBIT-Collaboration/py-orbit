"""
This module is a foil node class for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

# import injection class
from orbit.kickernodes import waveforms
from orbit.kickernodes import XKicker, YKicker


class TeapotXKickerNode(DriftTEAPOT):
	""" 
	The kicker node class for TEAPOT lattice
	"""
	def __init__(self, bunch, strength, waveform, name = "kicker"):
		"""
		Constructor. Creates the Kicker TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.kicker = XKicker(bunch, strength, waveform)
		self.setType("XKicker")
		self.setLength(0.0)

	def track(self, paramsDict):
		"""
		The kicker-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.kicker.kick()		
			
class TeapotYKickerNode(DriftTEAPOT):
	""" 
	The kicker node class for TEAPOT lattice
	"""
	def __init__(self, bunch, strength, waveform, name = "kicker"):
		"""
		Constructor. Creates the Kicker TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.kicker = YKicker(bunch,strength, waveform)
		self.setType("YKicker")
		self.setLength(0.0)
	
	def track(self, paramsDict):
		"""
		The kicker-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.kicker.kick()		
		
