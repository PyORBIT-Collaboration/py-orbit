#!/usr/bin/env python

"""
This is not a parallel version! 
"""

import math
import random
import sys

class rootTWaveform:
	""" 
	This class has horizontal sqrt(time) waveform.
	"""
	def __init__(self, syncpart, lattlength, duration, startamp, endamp):
		self.name = "RootTWaveform"
		self.syncpart = syncpart
		self.lattlength = lattlength
		self.duration = duration
		self.startamp = startamp
		self.endamp = endamp
	
	def getKickFactor(self):
		time = self.syncpart.time()
		amp = self.startamp - self.endamp;
		factor = amp * ( 1 - math.sqrt(time/self.duration)) + self.endamp
		return factor
	


class flatTopWaveform:
	""" 
	This class has the flat top (const) waveform.
	"""
	def __init__(self, amp):
		self.name = "FlatTopWaveform"
		self.amp = amp
	
	def getKickFactor(self):
		return self.amp

	