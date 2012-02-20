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


class TeapotInjectionlNode(DriftTEAPOT):
	""" 
	The foil node class for TEAPOT lattice
	"""
	def __init__(self, xmin, xmax, ymin, ymax, thick, name = "foil no name"):
		"""
		Constructor. Creates the Foil TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.foil = Foil(xmin, xmax, ymin, ymax, thick)
		self.setType("foil teapot")
		self.setLength(0.0)
		# The user choice of scattering routine. Defualt (0) is full scatter
		self.scatterChoice = 0

	def track(self, paramsDict):
		"""
		The foil-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		if(self.scatterChoice == 0):
			self.foil.traverseFoilFullScatter(bunch, lostbunch)
		else:
			self.foil.traverseFoilSimpleScatter(bunch)
		#put the track method here
		print "debug tracking the bunch through the foil name=",self.getName()," part ind=",self.getActivePartIndex()," length=",length
	
	
	def setScatterChoice(self, choice):
		self.scatterChoice = choice
		
