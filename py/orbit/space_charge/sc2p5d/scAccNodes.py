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

class SC2p5D_AccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalc2p5D wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2p5 SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("SC2p5D")
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
		
class SC2p5Drb_AccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeCalc2p5Drb wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2p5 SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("SC2p5Drb")
		self.pipe_radius = 1.0
		
	def setPipeRadius(self, pipe_radius):
		"""
		Assigns the pipe_radius for this space-charge node.
		"""
		self.pipe_radius = pipe_radius
		
	def getPipeRadius(self):
		"""
		Returns the pipe_radius of this space-charge node.
		"""		
		return self.pipe_radius
				
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		#print "debug ????????????? SC node=",self.getName()," scL=",self.getLengthOfSC()," pipe_r=",self.pipe_radius
		if(self.switcher != True): return
		bunch = paramsDict["bunch"]
		self.sc_calculator.trackBunch(bunch,self.sc_length,self.pipe_radius)
		


