#!/usr/bin/env python

#--------------------------------------------------------
# This is a collection of classes for describing quad 
# overlapping fields.
#--------------------------------------------------------
# The classes for the quads with the overlapping fields uses
# the field distribution from the following paper:
#    M.Berz, B. Erdelyn, K.Makino
#    Fringe Field Effects in Small Rings of Large Acceptance
#    Phys. Rev STAB, V3, 124001(2000)
#
#--------------------------------------------------------

import math
import sys
import os

from orbit.py_linac.lattice import BaseLinacNode, Drift, Quad
from orbit.py_linac.lattice import AxisFieldRF_Gap
from orbit.py_linac.lattice import AxisField_and_Quad_RF_Gap

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

from orbit_utils import Function

class AbstractQuadFieldSourceFunction:
	"""
	It is an abstract class describing the quadrupole magnetic 
	field as a function of the longitudinal coordinate.
	"""
	def __init__(self):
		pass
	
	def getLimitsZ(self):
		pass
	
	def getFuncValue(self,z):
		""" 
		Returns the quad's normalized field distribution 
		at the distance z from the center.
		the distribution should be normalized to 1 over the
		integration along the longitudinal coordinate.
		"""
		pass
	
class SimpleQuadFieldFunc(AbstractQuadFieldSourceFunction):
	"""
	It is an implementation of the QuadFieldSourceFunction class for
	a simple quad with constant field between (-L/2;+L/2).
	"""
	def __init__(self,quad):
		self.quad = quad
		
	def getLimitsZ(self):
		"""
		Returns (+L/2,-L/2) as longitudinal limits of the quad field.
		"""
		L = self.quad.getLength()
		return (-L/2,L/2)
		
	def getFuncValue(self,z):
		""" 
		Returns the quad's normalized field distribution at the distance z from the center 
		"""
		L = self.quad.getLength()
		if(abs(z) <= L/2):
			return 1./L
		return 0.
		
class EngeFunction(AbstractQuadFieldSourceFunction):
	""" 
	The Enge function with parameters from Berz's paper 
	M.Berz, B. Erdelyn, K.Makino
  Fringe Field Effects in Small Rings of Large Acceptance
  Phys. Rev STAB, V3, 124001(2000)	
	"""
	def __init__(self, length_param, acceptance_diameter_param, cutoff_level = 0.01):
		self.length = length_param
		self.acceptance_diameter = acceptance_diameter_param
		self.a_arr = [0.296471,4.533219,-2.270982,1.068627,-0.036391,0.022261]	
		self.normalization= 1.0
		self.n_func_points = 500
		#-----find cut-off z value
		self.cutoff_z = self.acceptance_diameter		
		step = self.acceptance_diameter
		self.cutoff_level = cutoff_level
		self.cutoff_z = self._findCutOff(step, cutoff_level)
		#------------------------------------------------------
		self.func = Function()
		self._normalize()
		
	def setEngeCoefficients(self,a_arr):
		"""
		Sets new values for Enge function's coeffients.
		"""
		self.a_arr = a_arr
		step = self.length/2
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)	
		self._normalize()
		
	def setCutOffLevel(self,cutoff_level):
		""" Sets the cutoff level for quad's  field """
		step = self.length/2
		self.cutoff_level = cutoff_level
		self.cutoff_z = self._findCutOff(step, cutoff_level)
		self._normalize()
		
	def setCutOffZ(self,cutoff_z):
		""" Sets the cutoff distance from the center of quad's field"""
		self.cutoff_z = cutoff_z
		self._normalize()
		
	def setLength(self,length):
		""" Sets the length of quad's field"""
		self.length = length
		step = self.length/2.0
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)
		self._normalize()
		
	def setAcceptanceDiameter(self,acceptance_diameter):
		""" Sets the acceptance diameter of the quad """
		self.acceptance_diameter = acceptance_diameter
		step = self.length/2.0
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)
		self._normalize()
		
	def setNumberOfPoints(self,n_func_points):
		""" Sets the number of points in the field function """
		self.n_func_points = n_func_points
		step = self.length/2.0
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)
		self._normalize()
		
	def getCuttOffZ(self):
		""" Returns the cutoff distance from the center of quad's field"""
		return self.cutoff_z 
		
	def getNumberOfPoints(self):
		""" Returns the number of points in the field function """
		return self.n_func_points
		
	def _getTrueEngeFunc(self, x):
		""" Returns the quad's field at the distance x from the center """
		# x is the distance from the center of the magnet with the iron length l """
		x = (math.fabs(x) - self.length/2.0)/self.acceptance_diameter
		sum_exp = self.a_arr[0]
		x0 = x
		for i in range(1,len(self.a_arr)):
			sum_exp += self.a_arr[i]*x0
			x0 *= x
		if(abs(sum_exp) > 30.): sum_exp = 30.0*sum_exp/abs(sum_exp)
		return self.normalization/(1.0+math.exp(sum_exp))

	def _findCutOff(self,step, cutoff_level):
		""" Finds the distance from the center where the field is less than cutoff level """
		self.normalization = 1.0
		init_val = self._getTrueEngeFunc(0.)
		z = step
		val = self._getTrueEngeFunc(z)/init_val
		if(val <= cutoff_level):
			return z
		while(val > cutoff_level):
			z += step
			val = self._getTrueEngeFunc(z)/init_val
		z0 = z - step
		z1 = z
		step_z = step/self.n_func_points
		val0 =  self._getTrueEngeFunc(z0)/init_val
		val1 =  self._getTrueEngeFunc(z1)/init_val		
		while(abs(z0-z1) > step_z):
			z_new = (z0+z1)/2.0
			val_new = self._getTrueEngeFunc(z_new)/init_val
			if(val_new <= cutoff_level):
				z1 = z_new
				val1 = val_new
			else:
				z0 = z_new
				val0 = val_new			
		self.cutoff_z = (z0+z1)/2.0
		return self.cutoff_z
				
	def _normalize(self):
		""" Normalizes the quad field function to the integral of 1 """
		self.normalization = 1.0
		step = self.cutoff_z/(self.n_func_points - 1)
		self.func.clean()
		sum_int = 0.
		for ind in range(self.n_func_points):
			z = step*ind
			val = self._getTrueEngeFunc(z)
			self.func.add(z,val)
			sum_int += val
		sum_int -= (self._getTrueEngeFunc(0.) + self._getTrueEngeFunc(step*(self.n_func_points - 1)))/2.0
		sum_int *= 2.0*step
		self.normalization = 1.0/sum_int
		self.func.setConstStep(1)
		
	def getFuncValue(self,z):
		""" Returns the quad's field at the distance z from the center """
		if(abs(z) >= self.func.getMaxX()): return 0.
		return self.normalization*self.func.getY(abs(z))
		
	def getLimitsZ(self):
		""" Returns the tuple with min and max Z value for this field """
		z_max = self.func.getMaxX()
		return (-z_max,z_max)
				
#-----------------------------------------------------------------------		
#-----Test of the Enge Function ----------------	
#-----------------------------------------------------------------------
if __name__ == "__main__":	
	#---- MEBT quads ----
	length_param = 0.066
	acceptance_diameter_param = 0.0363
	#---- DTL Permanent Quad
	length_param = 0.035
	acceptance_diameter_param = 0.025
	#--------------------
	cutoff_level = 0.01	
	func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
	z_max = func.getCuttOffZ()
	np = func.getNumberOfPoints()
	np = 50
	step = z_max/(np-1)
	for ind in range(2*np-1):
		z = ind*step-z_max
		val = func.getFuncValue(z)
		print " %12.5e  %12.5e "%(z*1000.,val)
	print "z limits=",func.getLimitsZ()
	func.setCutOffZ(0.5*func.getLimitsZ()[1])
	print "new z limits=",func.getLimitsZ()



