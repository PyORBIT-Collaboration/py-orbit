import sys
import os
import math

from orbit.utils import NamedObject, TypedObject

class Waveform(NamedObject, TypedObject):
	def __init__(self, name = "no name"):
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, "base waveform")
		
class KickerWaveform(Waveform):
	def __init__(self, name = "no name"):
		Waveform.__init__(self, name)
		self.setType("kicker waveform")
		self.kx = None
		self.ky = None
		
	def setKx(self,kx):
		self.kx = kx
	
	def getKx(self):
		return self.kx
		
	def setKy(self,ky):
		self.ky = ky
		
	def getKy(self):
		return self.ky

class ConstantKickerWaveform(KickerWaveform):
	def __init__(self, name = "no name"):
		KickerWaveform.__init__(self, name)
		
	def initialize(self,kx,ky):
		self.setKx(kx)
		self.setKy(ky)
		
	def calc(self, time):		
		pass
		
class SquareRootKickerWaveform(KickerWaveform):
	def __init__(self, name = "no name"):
		KickerWaveform.__init__(self, name)
		self.__initial = []
		
	def initialize(self,t1,t2,kx,ky,dx,dy):				
		self.__initial = (t1,t2,dx,dy)
		self.setKx(kx)
		self.setKy(ky)
		
	def calc(self, time):
		t = time
		(t1,t2,dx,dy) = self.__initial
		dt = sqrt((t-t1)/(t2-t1))
		self.setKx(self.getKx()*(1-dt)+dx)
		self.setKy(self.getKy()*(1-dt)+dy)		
		
class MagnetWaveform(Waveform):
	def __init__(self, name = "no name"):
		Waveform.__init__(self, name)
		self.setType("magnet waveform")
		self.strength = None
	
	def setStrength(self,strength):
		self.strength = strength

	def getStrength(self):
		return self.strength

class ConstantMagnetWaveform(MagnetWaveform):
	def __init__(self, name = "no name"):
		MagnetWaveform.__init__(self, name)
		
	def initialize(self,strength):
		self.setStrength(strength)
	
	def calc(self,time):
		pass

class LinearMagnetWaveform(MagnetWaveform):
	def __init__(self, name = "no name"):
		MagnetWaveform.__init__(self, name)
		self.__initial = []
		
	def initialize(self,t1,t2,s1,s2):
		self.__initial = (t1,t2,s1,s2)
		
	def calc(self, t):		
		(t1,t2,s1,s2) = self.__initial
		if t<t1: self.strength = s1
		elif t>t2: self.strength = s2
		else:
		  dt = math.sqrt((t-t1)/(t2-t1))
		  self.setStrength(s1+dt*(s2-s1))
