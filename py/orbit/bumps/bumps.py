#!/usr/bin/env python

"""
This is not a parallel version! 
"""
import math
import random
import sys
from bunch import Bunch

import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize, leastsq
from orbit.teapot import TEAPOT_MATRIX_Lattice


class simpleBump:
	""" 
	This routine adds a transverse particle coordinate bump.
	"""
	
	def __init__(self, bunch, xbump, xpbump, ybump, ypbump):
		
		self.bunch = bunch
		self.xbump = xbump
		self.ybump = ybump
		self.xpbump = xpbump
		self.ypbump = ypbump
            
	def bump(self):
		
		nparts = self.bunch.getSize();
		
		for i in range(nparts):
			newx = self.bunch.x(i) + self.xbump
			self.bunch.x(i, newx)
			newxp = self.bunch.xp(i) + self.xpbump
			self.bunch.xp(i, newxp)
			newy = self.bunch.y(i) + self.ybump
			self.bunch.y(i, newy)
			newyp = self.bunch.yp(i) + self.ypbump
			self.bunch.yp(i, newyp)

	def getLength(self):
		return 0



class close_orbit_bumps:
	""" 
	This routine try to find the kicks for a closed orbit bump.
	"""
	def __init__(self):
		self.tuneX = 0.0
		self.variables = []
		self.constraints = []
		
	def lattice_function(self,lattice,bunch,variables, constraints):
		
		self.variables = variables
		self.constraints = constraints
		
		nv = np.array(self.variables).shape[0]
		nc = np.array(self.constraints).shape[0]
		
		if not  nv == nc:
			print "Number of variables and constraints must be equal"
			print "Stop."
			sys.exit(0)
			
	
		matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)
		(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
		self.tuneX =  muX[-1][1]
		nodes = lattice.getNodes()
	
		mux_kick = []
		betax_kick = []
		alphax_kick = []
		for n in range(nv):
			for node in nodes:
				if node.getName() == self.variables[n][0][0]:
					for j in range(len(arrPosBetaX)):
						if (round(lattice.getNodePositionsDict()[node][1],4)==round(arrPosBetaX[j][0],4)):
							mux_kick.append(muX[j][1])
							betax_kick.append(arrPosBetaX[j][1])
							alphax_kick.append(arrPosAlphaX[j][1])
							self.variables[n].insert(3,[muX[j][1],arrPosBetaX[j][1],arrPosAlphaX[j][1]])
				
		mux_bpm = []
		betax_bpm = []
		alphax_bpm = []						
		for n in range(nc):
			for node in nodes:
				if node.getName() == self.constraints[n][0][0]:
					for j in range(len(arrPosBetaX)):
						if (round(lattice.getNodePositionsDict()[node][1],4)==round(arrPosBetaX[j][0],4)):
							mux_bpm.append(muX[j][1])
							betax_bpm.append(arrPosBetaX[j][1])
							alphax_bpm.append(arrPosAlphaX[j][1])
							self.constraints[n].insert(3,[muX[j][1],arrPosBetaX[j][1],arrPosAlphaX[j][1]])
							
		
				
	def get_lattice_tune(self):
		return self.tuneX
		
	def get_variables(self):
		return self.variables
		
	def get_constraints(self):
		return self.constraints
	
	#======================green function===============
	def cos_function(self,phi):
		if phi < 0:
			phi = phi + self.tuneX 
		return  math.cos((self.tuneX/2-phi)* 2.0*math.pi)
	#======================green function===============
	
	#======================green function===============
	def sin_function(self,phi):
		if phi < 0:
			phi = phi + self.tuneX 
		return  math.sin((self.tuneX/2-phi)* 2.0*math.pi)
	#======================green function===============
	
	
	#======================get xp with green function===============
	def get_xp(self,kick,k):	
		mu_s = self.constraints[k][2][0]
		beta_s = self.constraints[k][2][1]
		alpha_s = self.constraints[k][2][2]
	
		tmp =  0
		for i in range(0,4):
			mu_i = self.variables[i][2][0]
			beta_i = self.variables[i][2][1]
			alpha_i = self.variables[i][2][2]
			phi = mu_s-mu_i
			if phi < 0:
				phi = phi + self.tuneX
			tmp = tmp + kick[i]* math.sqrt(beta_i) * ( self.sin_function((mu_s-mu_i))  - alpha_s * self.cos_function((mu_s-mu_i)) )
		return tmp/math.sqrt(beta_s)/2/math.sin(self.tuneX/2*2*math.pi)
	#======================get posstion change due matrix===============
	
	#======================get x information with green function===============
	def get_x(self,kick,k):
		mu_s = self.constraints[k][2][0]
		beta_s = self.constraints[k][2][1]
		alpha_s = self.constraints[k][2][2]
		tmp = 0
		for i in range(len(kick)):
			mu_i = self.variables[i][2][0]
			beta_i = self.variables[i][2][1]
			alpha_i = self.variables[i][2][2]
			tmp = tmp + math.sqrt(beta_i) * kick[i] * self.cos_function((mu_s-mu_i))
		tmp = tmp*math.sqrt(beta_s)/2/math.sin(self.tuneX/2*2*math.pi)
		return tmp
	#======================get position information with green function===============
	
	#======================optimzation function===============	
	def optimzation_function(self,kick):
		nv = np.array(self.variables).shape[0]
		nc = np.array(self.constraints).shape[0]
		
		err = np.zeros(nc)
				
		# which coo. x or xp
		for i in range(nc):
			if self.constraints[i][1][0] == "x":
				x = self.get_x(kick,i)
				err[i] = self.constraints[i][1][1]-x
			elif self.constraints[i][1][0] == "xp":
				xp = self.get_xp(kick,i)
				err[i] = self.constraints[i][1][1]-xp
			else:
				print "Unknown coordinate"
				print "Stop."
				sys.exit(0)
		return LA.norm(err)
	#======================optimzation function===============
	
	def get_start_value(self):
		nv = np.array(self.variables).shape[0]
		p0 = np.zeros(nv)
		for i in range(nv):
			p0[i] = self.variables[i][1][1]
		return p0
		
		
	def find_kicks(self,method="Simplex"):
		p0 = self.get_start_value()
		if method == "Simplex":
			res = minimize(self.optimzation_function, p0, method='Nelder-Mead', options={'xtol': 1e-8, 'disp': None,'maxiter' : 20000})
			return res.x
		if method == "LeastSq":
			res = minimize(self.optimzation_function, p0, method='SLSQP')
			return res.x
		else:
			print "Unknown routine"
			print "Routines are Simplex or LeastSq"
			print "Stop."
			sys.exit(0)
	#======================optimzation function===============
