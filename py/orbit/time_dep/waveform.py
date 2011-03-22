import sys
import os
import math

from orbit.utils import NamedObject, TypedObject

class Waveform(NamedObject, TypedObject):
	def __init__(self, name = "no name"):
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, "waveform")		
		
class KickerWaveform(Waveform):
	def __init__(self):
		self.setType("kicker waveform")
		self.kx = None
		self.ky = None
		
	def getKx(self):
		return self.kx
		
	def getKy(self):
		return self.ky
		
	def setConstant(self, kx,ky):
		self.kx = kx
		self.ky = ky
		
	def setSquareRoot(self,t,t1,t2,kx,ky,dx,dy):
		dt = sqrt((t-t1)/(t2-t1))
		self.kx = kx*(1-dt)+dx
		self.ky = ky*(1-dt)+dy
		
class MagnetWaveform(Waveform):
	def __init__(self):
		self.setType("magnet waveform")
		self.strength = None
	
	def getStrength(self):
		return self.strength
		
	def setConstant(self, strength):
		self.strength = strength
				
	def setLinear(self, t,t1,t2,s1,s2):
		if t<t1: self.strength = t1
		elif t>t1: self.strength = t2
		else:
		  dt = sqrt((t-t1)/(t2-t1))
		  self.strength = s1+dt*(s2-s1)

