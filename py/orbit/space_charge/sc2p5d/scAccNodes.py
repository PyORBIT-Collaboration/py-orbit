"""
Module. Includes classes for all 2.5D space charge accelerator nodes.
"""

import sys
import os
import math

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

class SC2p5D_AccNode(AccNodeBunchTracker):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalc2p5Drb wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2p5 SC accelerator node element.
		"""
		AccNodeBunchTracker.__init__(self,name)
		self.setType("SC2p5D")
		self.boundary = None
		self.sc_length = 0.
		self.switcher = True
		self.sc_calculator = sc_calculator
		
	def setBoundary(self, boundary):
		"""
		Assigns the boundary for this space-charge node.
		"""
		self.boundary = boundary
		
	def getBoundary(self):
		"""
		Returns the boundary of this space-charge node.
		"""		
		return self.boundary
		
	def setLengthOfSC(self, sc_length):
		"""
		Defines the path length that will be used in SC kick calculations.
		"""
		self.sc_length = sc_length
		
	def getLengthOfSC(self):
		"""
		Returns the path length that will be used in SC kick calculations.
		"""		
		return self.sc_length
		
	def setCalculationOn(self, switcher):
		"""
		Sets the boolean parameter that define if the calculations will be performed. True of False.
		"""
		self.switcher = True
		
	def getCalculationOn(self):
		"""
		Returns the boolean parameter that define if the calculations will be performed. True of False.
		"""
		return 	self.switcher	
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		bunch = paramsDict["bunch"]
		if(self.boundary != None):
			self.sc_calculator.trackBunch(bunch,self.sc_length,self.boundary)
		else:
			self.sc_calculator.trackBunch(bunch,self.sc_length)
		
