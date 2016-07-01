#!/usr/bin/env python

"""
This is not a parallel version! 

"""
import math
import random
import sys
from bunch import Bunch


# ATTENTION !!! The python packet numpy and scipy are required
import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize, leastsq
from orbit.teapot import TEAPOT_MATRIX_Lattice

class orbit:
	def __init__(self, lattice, bunch):
		self.lattice = lattice
		self.bunch = bunch
		
	def offset(self):
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)

		OTM = np.zeros((6, 6))
		kickOTM = np.zeros((6))
		mt = matrix_lattice.getOneTurnMatrix()
		for i in range(6):
			for j in range(6):
				OTM[i][j] = mt.get(i,j)
			for i in range(6):
				kickOTM[i] = mt.get(i,6)

		#z0 fulfill: z0 = Mz0 with M as one turn matrix
		# Mz + k = z => z = (I-M)^-1 k
		z0, _, _, _ = np.linalg.lstsq(np.eye(6, dtype=int) - OTM,kickOTM)
		return z0
		
	def get_orbit(self):
	
		z0 =  self.offset()
		
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		OrbitX, OrbitY = matrix_lattice.getRingOrbit(z0)
		
		return OrbitX, OrbitY

class correction:
	""" 
	This routine corrected the closed orbit. Simple method. 
	"""
	
	def __init__(self, lattice, bunch):
		
		self.lattice = lattice
		self.bunch = bunch
		self.kicker_x = []
		self.kicker_y = []
		
		self.bpm_x = []
		self.bpm_y = []
		
		self.solution_x = 0.0
		self.solution_y = 0.0
		
	def get_solution_x(self):
		return self.solution_x
		
	def get_solution_y(self):
		return self.solution_y
		
	def orbit_corr(self):
	
		self.kicker_x, self.kicker_y = self.find_elements("monitor teapot", "kick teapot")
		self.bpm_x, self.bpm_y = self.find_elements("kick teapot", "monitor teapot")

		print "Numbers of monitor found", len(self.bpm_x)
		print "Numbers of corrector found", len(self.kicker_x)
		

		p0 = np.zeros(len(self.kicker_x))
		plsq_x = leastsq(self.optimzation_function, p0, args=(self.kicker_x,self.bpm_x))
		plsq_y = leastsq(self.optimzation_function, p0, args=(self.kicker_y,self.bpm_y))
		solution_x = plsq_x[0]
		solution_y = plsq_y[0]
		
		m = 0
		nodes = self.lattice.getNodes()
		for node in nodes:
			if node.getType() == "kick teapot":
				node.setParam("kx",solution_x[m])
				node.setParam("ky",solution_y[m])
				m = m + 1
		
	#======================green function===============
	def cos_function(self,phi,tune):
		if phi < 0:
			phi = phi + tune
		return  math.cos((tune/2-phi)* 2.0*math.pi)
	#======================green function===============
	
	#======================green function===============
	def sin_function(self,phi,tune):
		if phi < 0:
			phi = phi + tune
		return  math.sin((tune/2-phi)* 2.0*math.pi)
	#======================green function===============
	
		
	#======================get x information with green function===============
	def get_x(self,kick,j,bpm,kicker):
		mu_s = bpm[j][0]
		beta_s = bpm[j][1]
		alpha_s = bpm[j][2]
		tmp = 0
		tune = bpm[j][7]
		for i in range(len(kick)):
			mu_i = kicker[i][0]
			beta_i = kicker[i][1]
			alpha_i = kicker[i][2]
			tmp = tmp + math.sqrt(beta_i) * kick[i] * self.cos_function((mu_s-mu_i),tune)
		tmp = tmp*math.sqrt(beta_s)/2/math.sin(tune/2*2*math.pi)
		return tmp
	#======================get position information with green function===============


	def optimzation_function(self,kick,kicker,bpm):
		err = np.zeros(len(bpm))
		for j in range(len(bpm)):
			x = self.get_x(kick,j,bpm,kicker) 
			err[j] =  (x + bpm[j][6])**2
		return err

	def find_elements(self,count_el, find_el):
		
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
		(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
						
		OrbitX, OrbitY = orbit(self.lattice,self.bunch).get_orbit()
		bpm_x = []
		bpm_y = []
		eps_length = 1e-6
		m = 0
		pos_old = pos = 0.0
		nodes = self.lattice.getNodes()
		for node in nodes:
			print node.getType()
			if node.getType() == count_el:
				m = m + 1
			if node.getType() == find_el:
				for j in range(len(arrPosBetaX)):
					if (round(self.lattice.getNodePositionsDict()[node][1],4)==round(arrPosBetaX[j][0],4)):
						pos_old = pos
						pos = self.lattice.getNodePositionsDict()[node][1]
						if(abs(pos_old-pos) > eps_length):
							bpm_x.append([muX[j][1], arrPosBetaX[j][1], arrPosAlphaX[j][1], pos,m,node.getName(),OrbitX[j][1],muX[-1][1]])
							bpm_y.append([muY[j][1], arrPosBetaY[j][1], arrPosAlphaY[j][1], pos,m,node.getName(),OrbitY[j][1],muY[-1][1]])
		return bpm_x,bpm_y