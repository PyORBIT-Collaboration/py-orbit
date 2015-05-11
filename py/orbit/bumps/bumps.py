#!/usr/bin/env python

"""
This is not a parallel version! 
"""
import math
import random
import sys
from bunch import Bunch

class simpleBump:
	""" 
	This routine adds a transverse particle coordinate bump.
	"""
	
	def __init__(self, bunch, xbump, xpbump, ybump, ypbump):
		
		self.bunch = bunch
		self.xbump = xbump
		self.ybump = ybump
		self.xpbump = xpbump
		self.ypbump = ypbump
            
	def bump(self):
		
		nparts = self.bunch.getSize();
		
		for i in range(nparts):
			newx = self.bunch.x(i) + self.xbump
			self.bunch.x(i, newx)
			newxp = self.bunch.xp(i) + self.xpbump
			self.bunch.xp(i, newxp)
			newy = self.bunch.y(i) + self.ybump
			self.bunch.y(i, newy)
			newyp = self.bunch.yp(i) + self.ypbump
			self.bunch.yp(i, newyp)

	def getLength(self):
		return 0
