#====================================================================
#     Class BisectionSearchAlgorithm - 1D search
#====================================================================	

import os
import math
import sys
import time

from Solver import SearchAgorithm

# import the finalization function 
from orbit.utils import orbitFinalize

class BisectionSearchAlgorithm(SearchAgorithm):
	"""
	The Bisection search is a minimum finding method that applies 
	to any continuous functions. The function is one variable function.
	The search is performed inside the initial section which is
	defined by lowerLimit and upperLimit values in the initial
	TialPoit. The independent variable at the minimum of the function 
	should be between there limits.
	
	WARNING!
	This algorithm is for searching for a minimum of a unimodal function!
	It is possible that it will not work in a general case!
	"""
	def __init__(self):
		SearchAgorithm.__init__(self)
		self.setName("Bisection Search")
		self.score0 = None
		self.score1 = None
		self.x0 = None
		self.x1 = None
		self.lowerTrialPoint = None
		self.upperTrialPoint = None

	def setTrialPoint(self,initTrialPoint):
		res = SearchAgorithm.setTrialPoint(self,initTrialPoint)
		if(not res): return res
		lower_limit = self.initTrialPoint.getVariableProxyArr()[0].getLowerLimit()
		upper_limit = self.initTrialPoint.getVariableProxyArr()[0].getUpperLimit()
		self.lowerTrialPoint = self.initTrialPoint.getCopy()
		self.upperTrialPoint = self.initTrialPoint.getCopy()
		self.lowerTrialPoint.getVariableProxyArr()[0].setValue(lower_limit)
		self.upperTrialPoint.getVariableProxyArr()[0].setValue(upper_limit)
		scorer = self.solver.getScorer()
		self.x0 = lower_limit
		self.x1 = upper_limit
		self.score0 = scorer.getScore(self.lowerTrialPoint)
		self.score1 = scorer.getScore(self.upperTrialPoint)
		scoreBoard = self.solver.getScoreboard()
		scoreBoard.addScoreTrialPoint(self.score0,self.lowerTrialPoint)
		scoreBoard.addScoreTrialPoint(self.score1,self.upperTrialPoint)	
		return True
		
	def makeStep(self):
		#---- this is a middle point
		x = (self.x1 + self.x0)/2
		trialPoint = self.lowerTrialPoint.getCopy()
		trialPoint.getVariableProxyArr()[0].setValue(x)
		score = self.solver.getScorer().getScore(trialPoint)
		self.solver.getScoreboard().addScoreTrialPoint(score,trialPoint)
		if(self.score0 < self.score1):
			if(score <= self.score0):
				self.score1 = score
				self.x1 = x
				self.upperTrialPoint = trialPoint
			else:
				if(score <= self.score1):
					self.score1 = score
					self.x1 = x
					self.upperTrialPoint = trialPoint
				else:
					self.solver.getStopper().setShouldStop(True)
					return
		else:
			if(self.score0 == self.score1):
				if(score <= self.score0):
					self.score1 = score
					self.x1 = x
					self.upperTrialPoint = trialPoint
				else:
					self.solver.getStopper().setShouldStop(True)
					return
			else:
				if(score <= self.score0):
					self.score0 = score
					self.x0 = x
					self.lowerTrialPoint = trialPoint
				else:
					self.solver.getStopper().setShouldStop(True)
					return
		#-------------------------------------
		return