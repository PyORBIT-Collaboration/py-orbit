#!/usr/bin/env python

"""
This is not a parallel version! 
"""
import math
import random
import sys
from bunch import Bunch

class XKicker:
	""" 
	This routine injects particles into a bunch with user specified distribution
	functions.
	"""
	
	def __init__(self, bunch, strength, xWaveform):
		
		self.bunch = bunch
		self.waveform = xWaveform
		self.strength = strength
            
	def kick(self):
		
		nparts = self.bunch.getSize();
		kickfactor = self.waveform.getKickFactor()
		xkick = self.strength*kickfactor
		for i in range(nparts):
			newxp = self.bunch.xp(i) + xkick
			self.bunch.xp(i, newxp)

class YKicker:
	""" 
	This routine injects particles into a bunch with user specified distribution
	functions.
	"""
	
	def __init__(self, bunch, strength, yWaveform):
		
		self.bunch = bunch
		self.waveform = yWaveform
		self.strength = strength
	
	def kick(self):
		
		nparts = self.bunch.getSize();
		kickfactor = self.waveform.getKickFactor();
		ykick = self.strength*kickfactor
		for i in range(nparts):
			newyp = self.bunch.yp(i) + ykick
			self.bunch.yp(i, newyp)