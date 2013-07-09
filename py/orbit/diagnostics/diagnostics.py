#!/usr/bin/env python

"""
This is not a parallel version! 
"""


import math
import random
import sys
from bunch import BunchTwissAnalysis
from orbit.utils.consts import speed_of_light

class StatLats:
	""" 
	This class gathers delivers the statistical twiss parameters
	"""
	def __init__(self, filename):
		self.file_out = open(filename,"a")
		self.bunchtwissanalysis = BunchTwissAnalysis()
	
	def writeStatLats(self, s, bunch, lattlength = 0):
		self.bunchtwissanalysis.analyzeBunch(bunch)
		emitx = self.bunchtwissanalysis.getEmittance(0)
		betax = self.bunchtwissanalysis.getBeta(0)
		alphax = self.bunchtwissanalysis.getAlpha(0)
		betay = self.bunchtwissanalysis.getBeta(1)
		alphay = self.bunchtwissanalysis.getAlpha(1)
		emity = self.bunchtwissanalysis.getEmittance(1)
		dispersionx = self.bunchtwissanalysis.getDispersion(0, bunch)
		ddispersionx = self.bunchtwissanalysis.getDispersionDerivative(0, bunch)
		dispersiony = self.bunchtwissanalysis.getDispersion(1, bunch)
		ddispersiony = self.bunchtwissanalysis.getDispersionDerivative(1, bunch)
		
		sp = bunch.getSyncParticle()
		time = sp.time()
		if lattlength > 0:
			time = sp.time()/(lattlength/(sp.beta() * speed_of_light))

		self.file_out.write(str(s) + "\t" +  str(time) + "\t" + str(emitx)+ "\t" + str(emity)+ "\t" + str(betax)+ "\t" + str(betay)+ "\t" + str(alphax)+ "\t" + str(alphay) + "\t" + str(dispersionx) + "\t" + str(ddispersionx) + "\t" + str(dispersiony) + "\t" + str(ddispersiony) + "\n" )
		
	def closeStatLats(self):
		self.file_out.close()


class StatLatsSetMember:
	"""
	This class delivers the statistical twiss parameters
	"""
	def __init__(self, file):
		self.file_out = file
		self.bunchtwissanalysis = BunchTwissAnalysis()
	
	def writeStatLats(self, s, bunch, lattlength = 0):
		
		self.bunchtwissanalysis.analyzeBunch(bunch)
		emitx = self.bunchtwissanalysis.getEmittance(0)
		betax = self.bunchtwissanalysis.getBeta(0)
		alphax = self.bunchtwissanalysis.getAlpha(0)
		betay = self.bunchtwissanalysis.getBeta(1)
		alphay = self.bunchtwissanalysis.getAlpha(1)
		emity = self.bunchtwissanalysis.getEmittance(1)
		dispersionx = self.bunchtwissanalysis.getDispersion(0, bunch)
		ddispersionx = self.bunchtwissanalysis.getDispersionDerivative(0, bunch)
		dispersiony = self.bunchtwissanalysis.getDispersion(1, bunch)
		ddispersiony = self.bunchtwissanalysis.getDispersionDerivative(1, bunch)
		
		sp = bunch.getSyncParticle()
		time = sp.time()

		if lattlength > 0:
			time = sp.time()/(lattlength/(sp.beta() * speed_of_light))
	

		self.file_out.write(str(s) + "\t" +  str(time) + "\t" + str(emitx)+ "\t" + str(emity)+ "\t" + str(betax)+ "\t" + str(betay)+ "\t" + str(alphax)+ "\t" + str(alphay) + "\t" + str(dispersionx) + "\t" + str(ddispersionx) + "\t" + str(dispersiony) + "\t" + str(ddispersiony) + "\n")
	
	def closeStatLats(self):
		self.file_out.close()

class Moments:
	"""
		This class delivers the beam moments
	"""
	def __init__(self, filename, order):
		self.file_out = open(filename,"a")
		self.bunchtwissanalysis = BunchTwissAnalysis()
		self.order = order

	def writeMoments(self, s, bunch, lattlength = 0):
		
		sp = bunch.getSyncParticle()
		time = sp.time()
		if lattlength > 0:
			time = sp.time()/(lattlength/(sp.beta() * speed_of_light))
								 
		self.bunchtwissanalysis.analyzeBunch(bunch)
		self.bunchtwissanalysis.computeBunchMoments(bunch, self.order)
				
		self.file_out.write(str(s) + "\t" +  str(time) + "\t")						 
		for i in range(0,self.order+1):
			for j in range(0,i+1):
				self.file_out.write(str(self.bunchtwissanalysis.getBunchMoment(i-j,j)) + "\t")
		self.file_out.write("\n")
	
	def closeMoments(self):
		self.file_out.close()


class MomentsSetMember:
	"""
		This class delivers the beam moments
	"""
	def __init__(self, file, order):
		self.file_out = file
		self.order = order
		self.bunchtwissanalysis = BunchTwissAnalysis()
		
	def writeMoments(self, s, bunch, lattlength = 0 ):
		
		sp = bunch.getSyncParticle()
		time = sp.time()
	
		if lattlength > 0:
			time = sp.time()/(lattlength/(sp.beta() * speed_of_light))
	
		self.bunchtwissanalysis.analyzeBunch(bunch)
		self.bunchtwissanalysis.computeBunchMoments(bunch, self.order)
				
		self.file_out.write(str(s) + "\t" +  str(time) + "\t")
		for i in range(0,self.order+1):
			for j in range(0,i+1):
				self.file_out.write(str(self.bunchtwissanalysis.getBunchMoment(i-j,j)) + "\t")
		self.file_out.write("\n")

			
	def closeMoments(self):
		self.file_out.close()
