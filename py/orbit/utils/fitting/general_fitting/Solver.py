"""
This is a collection of classes for a general fitting problem.
At this moment there is only one fitting algorithm - simplex.
"""

import os
import math
import sys

from orbit.utils import NamedObject

# import the finalization function 
from orbit.utils import orbitFinalize

#====================================================================
#       class Solver
#====================================================================

class Solver:
	""" 
	The class is a main class of the general fitting package.
	It keeps references to all other components
	scoreboard
	search algorithm
	solve stopper
	scorer
	
	"""
	def __init__(self):
		"""
		Constructor of the main class of the general fitting package.
		"""
		self.scoreboard = Scoreboard(self)
		self.search_algorithm = SimplexAlgorithms()
		self.stopper = SolveStopperFactory.runForeverStopper()
		self.scorer = None
		self.parameter_arr = []
		self.is_running = False
		
	def getScoreboard(self):
		"""
		This method returns the Scoreboard.
		"""
		return self.scoreboard	
		
	def setScorer(self, scorer):
		"""
		This method sets the scorer.
		"""
		self.scorer = scorer
		
	def getScorer(self):
		"""
		This method returns the scorer.
		"""
		return self.scorer		
		
	def setAlgorithm(self, search_algorithm):
		"""
		This method sets the search algorithm.
		"""
		self.search_algorithm = search_algorithm
		
	def getAlgorithm(self):
		"""
		This method returns the search algorithm.
		"""
		return self.search_algorithm		
		
	def setStopper(self, stopper):
		"""
		This method sets the stopper.
		"""
		self.stopper = stopper
		
	def getStopper(self):
		"""
		This method returns the stopper.
		"""
		return self.stopper
		
	def isRunning(self):
		"""
		This method returns true or false.
		"""	
		
	def solve(self):
		"""
		This method applays the fitting algorithms to the problem.
		"""
		self.is_running = True
		
		self.scoreboard.init()
		
		self.search_algorithm.setSolver(self)
		
		while(not self.stopper.shouldStop()):
			self.search_algorithm.makeStep()
		
		self.is_running = False

#====================================================================
#       class ParameterProxy
#====================================================================

class ParameterProxy(NamedObject):
	"""
	This class represents the parameter for the score function in the fitting process. 
	"""
	def __init__(self, *arg, **kwargs):
		"""
		The constructor of the ParameterProxy class should have the following signatures
		ParameterProxy(parameterProxy_in)
		ParameterProxy(name, value, step)
		ParameterProxy(name = ""???"", value = ???, step = ???)
		"""
		self.setName("unknown")
		self.value = 0.
		self.step = 0.
		
		self.lowerLimit = -sys.float_info.max
		self.upperLimit = +sys.float_info.max
		
		if(len(arg) != 1):
			
			if(len(arg) == 1):
				paramProxy = arg[0]
				if(isinstance(paramProxy,ParameterProxy):
					self.setName(paramProxy.getName())
					self.value = paramProxy.value
					self.step = paramProxy.step
					self.lowerLimit = paramProxy.lowerLimit
					self.upperLimit = paramProxy.upperLimit
				else:
					msg = "ParameterProxy constructor. If it is only one argument it should be only ParameterProxy."
					msg = msg + os.linesep
					msg = "This argument is not!"
					msg = msg + os.linesep			
					msg = msg + "Stop."
					msg = msg + os.linesep
					orbitFinalize(msg)
			else:
				if(len(arg) == 3):
					self.setName(arg[0].getName())
					self.value = arg[1]
					self.step = arg[2]
				else:
					if(len(kwargs) == 3):
						self.setName(kwargs["name"])
						self.value = kwargs["value"]
						self.step = kwargs["step"]
					else:
						msg = "ParameterProxy constructor. It should be ParameterProxy(name = ""???"", value = ???, step = ???)"
						msg = msg + os.linesep		
						msg = msg + "Stop."
						msg = msg + os.linesep
						orbitFinalize(msg)
		else:
					msg = "ParameterProxy constructor. It should have one of the forms:"
					msg = msg + os.linesep
					msg = "1. ParameterProxy(parameterProxy_in)"
					msg = msg + os.linesep
					msg = "2. ParameterProxy(name, value, step)"
					msg = msg + os.linesep			
					msg = "3. ParameterProxy(name = ""???"", value = ???, step = ???)"
					msg = msg + os.linesep			
					msg = msg + "Stop."
					msg = msg + os.linesep
					orbitFinalize(msg)			

	def setValue(self, value):
		"""
		This method sets the value.
		"""
		self.value = value
		
	def getValue(self):
		"""
		This method returns the value.
		"""
		return self.value	
		
	def setStep(self, step):
		"""
		This method sets the step.
		"""
		self.step = step
		
	def getStep(self):
		"""
		This method returns the step.
		"""
		return self.step					
		
	def setLowerLimit(self, lowerLimit):
		"""
		This method sets the lowerLimit.
		"""
		self.lowerLimit = lowerLimit
		
	def getLowerLimit(self):
		"""
		This method returns the lowerLimit.
		"""
		return self.lowerLimit
		
	def setUpperLimit(self, upperLimit):
		"""
		This method sets the upperLimit.
		"""
		self.upperLimit = upperLimit
		
	def getUpperLimit(self):
		"""
		This method returns the upperLimit.
		"""
		return self.upperLimit			
		