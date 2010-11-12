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

class SC_Base_AccNode(AccNodeBunchTracker):
	"""
	The subclass of the AccNodeBunchTracker class. It is a base class for SC 2.5D nodes. It uses the c++ space charge calculator to calculate the .
	"""
	def __init__(self, sc_calculator, name = "no name"):			
		"""
		Constructor. Creates the 2p5 SC accelerator node element.
		"""
		AccNodeBunchTracker.__init__(self,name)
		self.setType("SC2p5D_Base")
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
		

