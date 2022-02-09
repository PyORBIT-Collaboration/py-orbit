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
	This routine injects particles into a bunch with user-
        specified distribution functions.
	"""

	def __init__(self, bunch, kx, xWaveform):
		self.bunch = bunch
		self.waveform = xWaveform
		self.kx = kx

	def kick(self):
		nparts = self.bunch.getSize();
		strength = self.waveform.getStrength()
		xkick = self.kx * strength
		for i in range(nparts):
			newxp = self.bunch.xp(i) + xkick
			self.bunch.xp(i, newxp)

class YKicker:
	"""
	This routine injects particles into a bunch with user-
        specified distribution functions.
	"""

	def __init__(self, bunch, ky, yWaveform):
		self.bunch = bunch
		self.waveform = yWaveform
		self.ky = ky

	def kick(self):
		nparts = self.bunch.getSize();
		strength = self.waveform.getStrength();
		ykick = self.ky * strength
		for i in range(nparts):
			newyp = self.bunch.yp(i) + ykick
			self.bunch.yp(i, newyp)
