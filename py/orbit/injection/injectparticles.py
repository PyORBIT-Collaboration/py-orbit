#!/usr/bin/env python

"""
This is not a parallel version! 
"""
import math
import random
import sys
from bunch import Bunch

class InjectParts:
	""" 
	This routine injects particles into a bunch with user specified distribution
	functions.
	"""
	
	def __init__(self, nparts, bunch, lostbunch, injectregion, xDistFunc, yDistFunc, lDistFunc):
		self.nparts = nparts
		self.bunch = bunch
		self.lostbunch = lostbunch
		self.injectregion = injectregion
		self.xDistFunc = xDistFunc
		self.yDistFunc = yDistFunc
		self.lDistFunc = lDistFunc
            
	def addParticles(self):
		(xmin,xmax,ymin,ymax) = self.injectregion
		
		for i in xrange(int(self.nparts)):
			(x,px) = self.xDistFunc.getCoordinates()
			(y,py) = self.yDistFunc.getCoordinates()
			(z,dE) = self.lDistFunc.getCoordinates()
			if((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)):	
				self.bunch.addParticle(x,px,y,py,z,dE)
			else:
				self.lostbunch.addParticle(x,px,y,py,z,dE)
				
		self.bunch.compress()
		self.lostbunch.compress()