#!/usr/bin/env python

"""
This is not a parallel version! 
"""


import math
import random
import sys

class JohoTransverse:
	""" 
	This class has the Joho distribution function generators in each plane.
	"""
	def __init__(self, order, alpha, beta, emitlim, centerpos=0, centermom=0, tailfraction=0, tailfactor=1):
		self.name = "JohoTransverse"
		self.order = order
		self.alpha = alpha
		self.beta = beta
		self.emitlim = emitlim
		self.centerpos = centerpos
		self.centermom = centermom
		self.tailfraction = tailfraction
		self.tailfactor = tailfactor
		self.__initialize()
		
	def __initialize(self):
		self.emit = self.emitlim * 2./(1. + self.order)
		self.pos = math.sqrt(self.emit)
		self.gamma = (1. + self.alpha * self.alpha)/self.beta
		self.mom = math.sqrt(self.emit * self.gamma)
		self.coschi = math.sqrt(1./(1. + (self.alpha * self.alpha)))
		self.sinchi = -self.alpha * self.coschi
		self.poslength = math.sqrt(self.emitlim * self.beta)
		self.momlength = math.sqrt(self.emitlim * self.gamma)
		self.orderinv = 1./self.order
		self.emitrms = 0.5 * self.emitlim/(1. + self.order)

	def getCoordinates(self):
	
		s1 = random.random()
		s2 = random.random()
		a = math.sqrt(1 - pow(s1, self.orderinv))
		al = 2. * math.pi * s2
		u = a * math.cos(al)
		v = a * math.sin(al)
		dpos = self.poslength * u
		dmom = self.momlength * (u * self.sinchi + v * self.coschi)
		if(self.tailfraction > 0.):
			if(random.random() < self.tailfraction):
				dpos *= self.tailfactor
				dmom *= self.tailfactor
		
				
		pos = self.centerpos + dpos
		mom = self.centermom + dmom
		return (pos,mom)
	
class JohoLongitudinal:
	""" 
        This class has the Joho distribution function generators in each plane.
        """
	def __init__(self, order, zlim, dElim, nlongbunches=0, deltazbunch=0, deltaznotch=0, tailfraction=0, tailfactor=1):
		self.name = "JohoLongitudinal"
		self.order = order
		self.zlim = zlim
		self.dElim = dElim
		self.nlongbunches = nlongbunches
		self.deltazbunch = deltazbunch
		self.deltaznotch = deltaznotch
		self.tailfraction = tailfraction
		self.tailfactor = tailfactor

	def getCoordinates(self):
		
		
		orderinv = 1./self.order
		s1 = random.random()
		s2 = random.random()
		a = math.sqrt(1. - pow(s1, orderinv))
		al = 2. * math.pi * s2
		u = a * math.cos(al)
		v = a * math.sin(al)
		zinj = self.zlim * u
		dEinj = self.dElim * v;
		factor = 360/248.0;
		if(self.tailfraction > 0.):
			if(random.random() < self.tailfraction):
				zinj *= self.ltailfraction
				dEinj *= self.ltailfraction
		
		if(self.nlongbunches > 1):
			ibunch = int(1 + self.nlongbunches * random.random())
			if self.nlongbunches < ibunch:
				ibunch = self.nlongbunches
			
			offset = (2. * ibunch - self.nlongbunches - 1)/2.
			ztemp = offset * self.deltazbunch
		
			if(self.deltaznotch != 0.):  
				while (ztemp < self.deltaznotch/2.) & (ztemp > -self.deltaznotch/2.):
					ibunch = int(1 + self.nlongbunches * random.random())
					if self.nlongbunches < ibunch:
						ibunch = self.nlongbunches
					offset = (2. * ibunch - self.nlongbunches - 1)/2.
					ztemp = offset * self.deltazbunch
			zinj += ztemp

		return (zinj,dEinj)



