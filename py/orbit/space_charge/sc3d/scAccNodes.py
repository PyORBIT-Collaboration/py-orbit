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

#import the base SC AccNode class
from orbit.space_charge.scAccNodes import SC_Base_AccNode

class SC3D_AccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalc3D wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2p5 SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("SC3D")
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		if(self.switcher != True): return
		bunch = paramsDict["bunch"]
		self.sc_calculator.trackBunch(bunch,self.sc_length)
		
class SC_UniformEllipses_AccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalcUnifEllipse wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the Uniform Ellipses SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("UnifEllsSC")
		self.pipe_radius = 1.0
						
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Uniform Ellipses Space Charge calculator.
		"""
		#print "debug ????????????? SC node=",self.getName()," scL=",self.getLengthOfSC()
		if(self.switcher != True): return
		bunch = paramsDict["bunch"]
		self.sc_calculator.trackBunch(bunch,self.sc_length)		
		

