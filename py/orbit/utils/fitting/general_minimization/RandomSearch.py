#!/usr/bin/env python
#
# Random search algorithms
#
# Coder: Andrei Shishlo
# Algorithm is a copy of one implemented by Tom Pelaia
#-------------------------------------------

import math
import sys
import copy
import random

from Solver import SearchAgorithm

#====================================================================
#       class RandomSearchAlgorithm
#====================================================================

class RandomSearchAlgorithm(SearchAgorithm):
	""" 
	The RandomSearchAlgorithm uses floating windows for each variables 
	with a defined shrinkage factor. The source of algorithm is unknown 
	to me. It is a copy of XALDEV 2nd generation optimizer algorithm.
	It was a copy of Tom Pelaia code.
	"""
	
	def __init__(self):
		SearchAgorithm.__init__(self)
		self.setName("Random Search")
		self.initTrialPoint = None
		#shrinkage factor of coordinates' steps
		self.shrinkageFactor = 3.0
		#----- internal arrays of parameters
		self.coords = []
		self.coords_old = []
		#--- coordinates windows
		self.coords_low = []
		self.coords_upp = []
		self.step_arr = []
		#---- Numbers of variables for fiiting
		self.nD = 0
		#---- best score
		self.best_score = sys.float_info.max
    
    
	def _reset(self):
		"""
		Resets the searching process from scratch; forget history.
		"""
		#----- internal arrays of parameters
		self.coords = []
		self.coords_old = []
		#--- coordinates windows
		self.coords_low = []
		self.coords_upp = []
		#---- Numbers of variables for fiiting
		self.nD = 0		
		if(self.initTrialPoint == None): return False
		self.coords = self.initTrialPoint.getVariablesUsedInOptArr()
		self.coords_old = self.coords[:]
		self.nD = len(self.coords)
		self.step_arr = self.initTrialPoint.getStepsUsedInOptArr()
		for ind in range(self.nD):
			variableProxy = self.initTrialPoint.getVariableProxyArr()[ind]
			val_LowLimit = variableProxy.getLowerLimit()
			val_UppLimit = variableProxy.getUpperLimit()
			val_low = max(self.coords[ind] - self.step_arr[ind],val_LowLimit)
			val_upp = min(self.coords[ind] + self.step_arr[ind],val_UppLimit)
			self.coords_low.append(val_low)
			self.coords_upp.append(val_upp)
		#----------------------------------------------
		score = self._testFunc(self.coords)
		if(score == None): return False
		if(score < self.best_score): self.best_score = score
		return True

	def _testFunc(self,guess):
		"""
		Calculates the score for particular variables. 
		"""
		trialPoint = self.initTrialPoint.getCopy()
		trialPoint.setVariablesUsedInOptArr(guess)	
		trialPoint.setStepsUsedInOptArr(self.step_arr)
		if(not trialPoint.isAcceptable()):
			return None
		if(self.solver.getStopper().getShouldStop()):
			return None
		score = self.solver.getScorer().getScore(trialPoint)
		scoreBoard = self.solver.getScoreboard()
		scoreBoard.addScoreTrialPoint(score,trialPoint)
		self.solver.getStopper().checkStopConditions(self.solver)
		if(self.solver.getStopper().getShouldStop()):
			return None
		return score	
		
	def _isTrialPointAcceptable(self,guess):
		"""
		Checks if the values of variables are acceptable in terms of limits.
		"""
		trialPoint = self.initTrialPoint.getCopy()
		trialPoint.setVariablesUsedInOptArr(guess)	
		trialPoint.setStepsUsedInOptArr(self.step_arr)
		return trialPoint.isAcceptable()		
		
	def setShrinkageFactor(self,shrinkageFactor):
		"""
		Sets the coefficient for parameters shrinkage
		"""
		self.shrinkageFactor = shrinkageFactor
		
	def getShrinkageFactor(self):
		"""
		returns the coefficient for parameters shrinkage
		"""
		return self.shrinkageFactor	
		
	def setSolver(self,solver):
		"""
		Sets the solver instance for the search algorithm.
		"""
		self.solver = solver
		
	def setTrialPoint(self,initTrialPoint):
		"""
		The initial preparation for a loop with makeStep() calls inside solver
		"""
		res = SearchAgorithm.setTrialPoint(self,initTrialPoint)
		if(not res): return res
		#----- set up all initial arrays and calculate the first score
		res = self._reset()
		if(not res): return res
		return True
		
	def init(self):
		if(self.initTrialPoint == None or self.solver == None): return False
		return True
	
	def makeStep(self):
		"""
		Implementation of the abstract method of the parent class.
		Perform the one step of the fitting algorithm.
		"""
		if(self.nD <= 0):
			self.solver.getStopper().setShouldStop(True)
			return		
		#---- make a new set of variable values
		changeProbabilityBase = 1.0/self.nD
		expectedNumToChange = 1.
		newPointDone = False
		while(not newPointDone):
			changeProbability = expectedNumToChange * changeProbabilityBase
			coordChanged = False
			for ind in range(self.nD):
				if(random.random() <= changeProbability):
					self.coords[ind] = self.coords_low[ind] + (self.coords_upp[ind] - self.coords_low[ind])*random.random()
					coordChanged = True
				else:
					self.coords[ind] = self.coords_old[ind]
			if(not coordChanged):
				expectedNumToChange += random.randint(0,self.nD) + 1
			else:
				newPointDone = True
		#-------------------------------------
		#---- if new coordinates are bad we will try again
		if(not self._isTrialPointAcceptable(self.coords)):
			self.coords = self.coords_old[:]
			return 
		#-------------------------------------
		score = self._testFunc(self.coords)	
		if(score == None):
			self.solver.getStopper().setShouldStop(True)
			return
		if(score < self.best_score): 
			self.best_score = score
			self._shrinkWindow(self.shrinkageFactor)
			self.coords_old = self.coords[:]			
		else:
			self.coords = self.coords_old[:]
		#-----------------------------------------
		return
		
	def _shrinkWindow(self,shrinkageFactor):
		"""
		It will shrink the delta between upper and lower values for some variables.
		It seems that it is making delta bigger, but it is wrong impression. 
		"""
		trialPoint = self.solver.getScoreboard().getBestTrialPointReference()
		self.step_arr = trialPoint.getStepsUsedInOptArr()
		for ind in range(self.nD):
			if(self.coords[ind] == self.coords_old[ind]): continue
			variableProxy = trialPoint.getVariableProxyArr()[ind]
			val_LowLimit = variableProxy.getLowerLimit()
			val_UppLimit = variableProxy.getUpperLimit()
			step = shrinkageFactor*abs(self.coords[ind] - self.coords_old[ind])
			self.step_arr[ind] = step
			self.coords_low[ind] = max(self.coords[ind] - self.step_arr[ind],val_LowLimit)
			self.coords_upp[ind] = min(self.coords[ind] + self.step_arr[ind],val_UppLimit)
		self.initTrialPoint.setStepsUsedInOptArr(self.step_arr)
		trialPoint.setStepsUsedInOptArr(self.step_arr)
