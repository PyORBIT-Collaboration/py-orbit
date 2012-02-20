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
	
	def __init__(self, nparts, bunch, lostfoilbunch, foilparams, xDistFunc, yDistFunc, lDistFunc):
		self.nparts = nparts
		self.bunch = bunch
		self.lostfoilbunch = lostfoilbunch
		self.foilparams = foilparams
		self.xDistFunc = xDistFunc
		self.yDistFunc = yDistFunc
		self.lDistFunc = lDistFunc
		self.(xmin,xmax,ymin,ymax) = foilparams

	def addParticles(self):
		(x,px) = self.xDistFunc
		(y,py) = self.yDistFunc
		(z,dE) = self.lDistFunc
		
		for i in xrange(nParts):
			if((x > self.xmin) && (x < self.xmax) && (y > self.ymin) && (y < self.ymax)): 
				self.bunch.addParticle(x,px,y,py,z,dE)
			else:
				print 'paticle outside of foil aperture '
			
		self.bunch.compress()
