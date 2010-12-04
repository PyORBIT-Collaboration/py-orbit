"""
Module. Includes abstract classes for all types of space charge accelerator nodes.
"""

import sys
import os
import math

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

class SC_Base_AccNode(AccNodeBunchTracker):
	"""
	The subclass of the AccNodeBunchTracker class. It is a base class for SC nodes. It uses the c++ space charge calculator to calculate the .
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the Space Charge (SC) accelerator node element.
		"""
		AccNodeBunchTracker.__init__(self,name)
		self.setType("SC_Base")
		self.sc_length = 0.
		self.switcher = True
		self.sc_calculator = sc_calculator
		
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
		
	def getCalculator(self):
		"""
		Returns the space-charge calculator.
		"""
		return self.sc_calculator
		
	def trackDesign(self, paramsDict):
		"""
		This method is for Linac Nodes compatibility. It is empty and should not be used for Space Charge calculations.
		"""
		pass	
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator. Each subclass should implement this method.
		"""
		pass
	

