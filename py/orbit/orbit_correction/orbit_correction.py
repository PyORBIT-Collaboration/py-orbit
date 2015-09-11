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
		
	def orbit_corr(self):
	
		self.kicker_x, self.kicker_y = self.find_elements("monitor", "kicker")
		self.bpm_x, self.bpm_y = self.find_elements("kicker", "monitor")

		print "Numbers of monitor found", len(self.bpm_x)
		print "Numbers of corrector found", len(self.kicker_x)
		
		#for i in range(len(self.kicker_x)):
		#	print self.kicker_x[i]
	
		#for i in range(len(self.bpm_x)):
		#	print self.bpm_x[i]
		

		p0 = np.zeros(len(self.kicker_x))
		plsq_x = leastsq(self.optimzation_function, p0, args=(self.kicker_x,self.bpm_x))
		plsq_y = leastsq(self.optimzation_function, p0, args=(self.kicker_y,self.bpm_y))
		solution_x = plsq_x[0]
		solution_y = plsq_y[0]
		

		m = 0
		nodes = self.lattice.getNodes()
		for node in nodes:
			if node.getType() == "kicker":
				node.setParam("kx",solution_x[m])
				node.setParam("ky",solution_y[m])
				m = m + 1
	
	def get_position(self,j,n,kick, bpm,kicker):
		mu_s = bpm[j][0]
		beta_s = bpm[j][1]
		alpha_s = bpm[j][2]
		tmp1 = tmp2 = 0
		for i in range(n):
			mu_i = kicker[i][0]
			beta_i = kicker[i][1]
			alpha_i = kicker[i][2]
			phi = mu_s-mu_i
			if phi < 0:
				phi = phi + tuneX
			tmp1 = tmp1 + math.sqrt(beta_i*beta_s) * kick[i] * math.sin(2.0*math.pi*(phi)) #
			tmp2 = tmp2 + kick[i]* math.sqrt(beta_i/beta_s) * (math.cos(2.0*math.pi*(phi))  - alpha_s * math.sin(2.0*math.pi*(phi)))
		return tmp1


	def optimzation_function(self,kick,kicker,bpm):
		err = np.zeros(len(bpm))
		# which coo. x or xp
		for j in range(len(bpm)):
			x = self.get_position(j,bpm[j][6],kick, bpm,kicker)
			err[j] = (bpm[j][3] + x)**2
		return err

	def find_elements(self,count_el, find_el):
		
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
		(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
		
		bpm_x = []
		bpm_y = []
		eps_length = 1e-6
		m = 0
		pos_old = pos = 0.0
		nodes = self.lattice.getNodes()
		for node in nodes:
			if node.getType() == count_el:
				m = m + 1
			if node.getType() == find_el:
				for j in range(len(arrPosBetaX)):
					if (round(self.lattice.getNodePositionsDict()[node][1],4)==round(arrPosBetaX[j][0],4)):
						pos_old = pos
						pos = self.lattice.getNodePositionsDict()[node][1]
						if(abs(pos_old-pos) > eps_length):
							mt = matrix_lattice.makeMatrix(pos) 
							bpm_x.append([muX[j][1], arrPosBetaX[j][1], arrPosAlphaX[j][1], mt.get(0,6), mt.get(1,6) ,pos,m,node.getName()])
							bpm_y.append([muY[j][1], arrPosBetaY[j][1], arrPosAlphaY[j][1], mt.get(2,6), mt.get(3,6) ,pos,m,node.getName()])
		return bpm_x,bpm_y