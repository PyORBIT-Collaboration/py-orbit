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

#import the base DirectForce AccNode class
from orbit.space_charge.scAccNodes import SC_Base_AccNode

class DirectForce2p5D_AccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalc2p5D wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2p5 SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("DirectForce2p5D")
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		if(self.switcher != True): return
		bunch = paramsDict["bunch"]
		self.sc_calculator.trackBunch(bunch,self.sc_length)


