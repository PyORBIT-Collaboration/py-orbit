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

class TeapotCollimatorNode(DriftTEAPOT):
	""" 
	The collimator node class for TEAPOT lattice
	"""
	def __init__(self, name = "collimator no name"):
		"""
		Constructor. Creates the Collimator TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.setType("collimator teapot")

	def track(self, paramsDict):
		"""
		The collimator-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here
		print "debug tracking the bunch through the collimator name=",self.getName()," part ind=",self.getActivePartIndex()," length=",length
	
	

