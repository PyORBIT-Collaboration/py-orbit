#!/usr/bin/env python

#--------------------------------------------------------
# This is a collection of aperture classes for linac lattices
#--------------------------------------------------------

import math
import sys
import os

from orbit.py_linac.lattice import BaseLinacNode, Quad

# import injection C++ class
from aperture import Aperture

class LinacApertureNode(BaseLinacNode):
	"""
	The aperture classes removes particles from bunch and places them in the lostbunch
	if their coordinates are not inside the aperture:
	The shape variable could be:
	1 is circle (a is a radius)
	2 is elipse (a and b are a half-axises)
	3 is rectangle (a and b are a half-horizontal and vertical sizes)
	c and d parameters are x and y offsets of the center
	"""
	def __init__(self, shape, a, b, pos = 0., c = 0., d = 0., name = "aperture"):
		BaseLinacNode.__init__(self,name)
		self.shape = shape
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.pos = pos
		self.aperture = Aperture(self.shape, self.a, self.b, self.c, self.d, self.pos)
	
	def track(self, paramsDict):
		bunch = paramsDict["bunch"]
		if(paramsDict.has_key("lostbunch")):
			lostbunch = paramsDict["lostbunch"]
			self.aperture.checkBunch(bunch, lostbunch)
		else:
			self.aperture.checkBunch(bunch)

	def trackDesign(self, paramsDict):
		"""
		This method does nothing for the aperture case.
		"""
		pass

	def setPosition(self, pos):
		self.pos = pos
		self.Aperture.setPosition(self.pos)

	def getPosition(self):
		return self.pos

class CircleLinacApertureNode(LinacApertureNode):
	"""
	The curcular aperture shape = 1
	"""
	def __init__(self, a, pos = 0., c = 0., d = 0., name = "aperture"):
		LinacApertureNode.__init__(self,1,a,a,pos,c,d,name)


class EllipseLinacApertureNode(LinacApertureNode):
	"""
	The ellipse aperture shape = 2
	"""
	def __init__(self, a, b, pos = 0., c = 0., d = 0., name = "aperture"):
		LinacApertureNode.__init__(self,2,a,b,pos,c,d,name)
		

class RectangleLinacApertureNode(LinacApertureNode):
	"""
	The rectangle aperture shape = 3
	"""
	def __init__(self, a, b, pos = 0., c = 0., d = 0., name = "aperture"):
		LinacApertureNode.__init__(self,3,a,b,pos,c,d,name)





