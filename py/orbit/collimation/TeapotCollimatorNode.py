"""
This module is a collimator node class for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

# import Collimator class
from collimator import Collimator

class TeapotCollimatorNode(DriftTEAPOT):
	""" 
	The collimator node class for TEAPOT lattice
	"""
	def __init__(self, length, ma, density_fac, shape, a, b, c, d, angle, name = "collimator no name"):
		"""
		Constructor. Creates the Collimator TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.collimator = Collimator(length,ma,density_fac,shape,a,b,c,d,angle)
		self.setType("collimator teapot")
		self.setLength(length)

	def track(self, paramsDict):
		"""
		The collimator-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		self.collimator.collimateBunch(bunch, lostbunch)
		
			

