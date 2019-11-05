#====================================================================
#     Class GoldenSectionSearchAlgorithm - 1D search
#====================================================================	

import os
import math
import sys
import time

from Solver import SearchAgorithm

# import the finalization function 
from orbit.utils import orbitFinalize

class GoldenSectionSearchAlgorithm(SearchAgorithm):
	"""
	The GoldenSection search is a minimum finding method that applies 
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
		self.setName("Golden Section Search")
		self.scoreA = None
		self.score1 = None
		self.score2 = None
		self.scoreB = None		
		self.xA = None
		self.x1 = None
		self.x2 = None
		self.xB = None
		self.trialPointA = None
		self.trialPoint1 = None
		self.trialPoint2 = None
		self.trialPointB = None
		self.phi = (1.0+math.sqrt(5.0))/2
		
	def setTrialPoint(self,initTrialPoint):
		res = SearchAgorithm.setTrialPoint(self,initTrialPoint)
		if(not res): return res
		lower_limit = self.initTrialPoint.getVariableProxyArr()[0].getLowerLimit()
		upper_limit = self.initTrialPoint.getVariableProxyArr()[0].getUpperLimit()
		self.xA = lower_limit
		self.xB = upper_limit
		self.x1 = self.xB - (self.xB - self.xA)/self.phi
		self.x2 = self.xA + (self.xB - self.xA)/self.phi
		self.trialPointA = self.initTrialPoint.getCopy()
		self.trialPoint1 = self.initTrialPoint.getCopy()
		self.trialPoint2 = self.initTrialPoint.getCopy()
		self.trialPointB = self.initTrialPoint.getCopy()
		self.trialPointA.getVariableProxyArr()[0].setValue(self.xA)
		self.trialPoint1.getVariableProxyArr()[0].setValue(self.x1)
		self.trialPoint2.getVariableProxyArr()[0].setValue(self.x2)
		self.trialPointB.getVariableProxyArr()[0].setValue(self.xB)
		scorer = self.solver.getScorer()
		self.scoreA = scorer.getScore(self.trialPointA)
		self.score1 = scorer.getScore(self.trialPoint1)
		self.score2 = scorer.getScore(self.trialPoint2)
		self.scoreB = scorer.getScore(self.trialPointB)
		scoreBoard = self.solver.getScoreboard()
		scoreBoard.addScoreTrialPoint(self.scoreA,self.trialPointA)
		scoreBoard.addScoreTrialPoint(self.score1,self.trialPoint1)
		scoreBoard.addScoreTrialPoint(self.score2,self.trialPoint2)
		scoreBoard.addScoreTrialPoint(self.scoreB,self.trialPointB)	
		return True

	def makeStep(self):
		score_arr = [self.scoreA,self.score1,self.score2,self.scoreB]
		ind_max = 0
		score_max = self.scoreA
		for ind in range(1,4):
			if(score_arr[ind] >= score_max):
				ind_max = ind
				score_max = score_arr[ind]
		#----------------------
		if(ind_max == 1 or ind_max == 2):
			self.solver.getStopper().setShouldStop(True)
			return
		if(ind_max == 0):
			self.xA = self.x1
			self.trialPointA = self.trialPoint1
			self.scoreA = self.score1
			self.x1 = self.x2
			self.trialPoint1 = self.trialPoint2
			self.score1 = self.score2
			self.x2 = self.xA + (self.xB - self.xA)/self.phi
			self.trialPoint2 = self.trialPoint1.getCopy()
			self.trialPoint2.getVariableProxyArr()[0].setValue(self.x2)
			self.score2 = self.solver.getScorer().getScore(self.trialPoint2)
			self.solver.getScoreboard().addScoreTrialPoint(self.score2,self.trialPoint2)
			return
		if(ind_max == 3):
			self.xB = self.x2
			self.trialPointB = self.trialPoint2
			self.scoreB = self.score2
			self.x2 = self.x1
			self.trialPoint2 = self.trialPoint1
			self.score2 = self.score1
			self.x1 = self.xB - (self.xB - self.xA)/self.phi
			self.trialPoint1 = self.trialPoint2.getCopy()
			self.trialPoint1.getVariableProxyArr()[0].setValue(self.x1)
			self.score1 = self.solver.getScorer().getScore(self.trialPoint1)
			self.solver.getScoreboard().addScoreTrialPoint(self.score1,self.trialPoint1)
		return
