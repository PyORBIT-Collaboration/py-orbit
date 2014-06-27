"""
Module. Includes classes for the 2D slice-by-slice space charge accelerator nodes.
"""

import sys
import os
import math

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the base SC AccNode class
from orbit.space_charge.scAccNodes import SC_Base_AccNode

class SC2DSliceBySlice_AccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalcSliceBySlice2D wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2D slice-by-slice SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("SC2DSliceBySlice")
		self.boundary = None
		
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
				
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		if(self.switcher != True): return
		bunch = paramsDict["bunch"]
		if(self.boundary != None):
			self.sc_calculator.trackBunch(bunch,self.sc_length,self.boundary)
		else:
			self.sc_calculator.trackBunch(bunch,self.sc_length)
		

