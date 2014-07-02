"""
This module is a Aperture node class for TEAPOT lattice
"""

import os
import math

from bunch import Bunch
import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT
                                     
# import injection class
from aperture import Aperture

#The aperture class simply removes particles from on bunch and places them in the lostbunch.

#This create an aperture node. The shape variable should be a number, 1 is circle, 2 is elipse, and 3 is rectangle. a is the first dimension, either the radius for a circle, or the half length in the x diminsion. b is the y half length of the aperture and does nothing for a circle.  c is the x offset and d is the y offset of the aperture.
class TeapotApertureNode(DriftTEAPOT):
	def __init__(self, shape, a, b, pos = 0, c = 0, d = 0, name = "aperture"):
		DriftTEAPOT.__init__(self,name)
		self.shape = shape
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.pos = pos
		self.Aperture = Aperture(self.shape, self.a, self.b, self.c, self.d, self.pos)
	
	def track(self, paramsDict):
		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		self.Aperture.checkBunch(bunch, lostbunch)

	def setPosition(self, pos):
		self.pos = pos
		self.Aperture.setPosition(self.pos)

#This create a circular aperture node. a the radius for a circle. c is the x offset and d is the y offset of the aperture.
class CircleApertureNode(DriftTEAPOT):
	def __init__(self, a, pos = 0, c = 0, d = 0, name = "aperture"):
		DriftTEAPOT.__init__(self,name)
		self.shape = 1
		self.a = a
		self.b = 1
		self.c = c
		self.d = d
		self.pos = pos
		self.Aperture = Aperture(self.shape, self.a, self.b, self.c, self.d, self.pos)
	
	def track(self, paramsDict):

		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		self.Aperture.checkBunch(bunch, lostbunch)
	
	def setPosition(self, pos):
		self.pos = pos
		self.Aperture.setPosition(self.pos)


#This create an elpitical aperture node. a is the the half length in the x diminsion. b is the y half length of the aperture.  c is the x offset and d is the y offset of the aperture.
class EllipseApertureNode(DriftTEAPOT):
	def __init__(self, a, b, pos = 0, c = 0, d = 0,  name = "aperture"):
		DriftTEAPOT.__init__(self,name)
		self.shape = 2
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.pos = pos
		self.Aperture = Aperture(self.shape, self.a, self.b, self.c, self.d, self.pos)
	
	def track(self, paramsDict):
		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		self.Aperture.checkBunch(bunch, lostbunch)

	def setPosition(self, pos):
		self.pos = pos
		self.Aperture.setPosition(self.pos)

#This create an rectangular aperture node. a is the the half length in the x diminsion. b is the y half length of the aperture.  c is the x offset and d is the y offset of the aperture.
class RectangleApertureNode(DriftTEAPOT):
	def __init__(self, a, b, pos = 0, c = 0, d = 0, name = "aperture"):
		DriftTEAPOT.__init__(self,name)
		self.shape = 3
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.pos = pos
		self.Aperture = Aperture(self.shape, self.a, self.b, self.c, self.d, self.pos)
	
	def track(self, paramsDict):
		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		self.Aperture.checkBunch(bunch, lostbunch)

	def setPosition(self,pos):
		self.pos = pos
		self.Aperture.setPosition(self.pos)
