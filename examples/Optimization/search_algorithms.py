#-----------------------------------------------------------------------		
#-----Test of the Search Algorithms
#  1D Search:
#      BisectionSearchAlgorithm
#      GoldenSectionSearchAlgorithm
#-------------------------------------------
#  General Search (1d,2D, ...):
#      Simplex
#      Random search
#-----------------------------------------------------------------------
	
import os
import math
import sys
import time	
	
from orbit.utils.fitting.general_minimization.Solver import Solver, Scorer, SolveStopperFactory, VariableProxy, TrialPoint

from orbit.utils.fitting.general_minimization.BisectionSearch1D import BisectionSearchAlgorithm
from orbit.utils.fitting.general_minimization.GoldenSectionSearch1D import GoldenSectionSearchAlgorithm
from orbit.utils.fitting.general_minimization.SimplexSearch import SimplexSearchAlgorithm
from orbit.utils.fitting.general_minimization.RandomSearch import RandomSearchAlgorithm

#------------------------------------------------------
#  Functions for minimization
#------------------------------------------------------

class SquareValueScorer1D(Scorer):
	"""
	The implementation of the abstract Score class 
	for search algorithms testing.
	"""
	def __init__(self):
		self.min_pos = 2.0
		self.base = 10.
	
	def getScore(self,trialPoint):
		x = trialPoint.getVariableProxyArr()[0].getValue()
		y = (x-self.min_pos)**2 + self.base
		return y

	def getAnswer(self):
		"""
		Returns true position of the minimum and base value.
		"""
		return (self.min_pos,self.base)

class SquareValueScorer3D(Scorer):
	"""
	The implementation of the abstract Score class 
	for search algorithms testing.
	"""
	def __init__(self):
		self.min_pos_arr = [2.0,3.0,4.0]
		self.base = 5.
	
	def getScore(self,trialPoint):
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()
		x_arr = self.min_pos_arr
		y = (x0-x_arr[0])**2 + (x1-x_arr[1])**2 + (x2-x_arr[2])**2 + self.base
		return y

	def getAnswer(self):
		"""
		Returns true position of the minimum and base value.
		"""
		return (self.min_pos_arr,self.base)

#-------------------------------------------------------
# Convenience function for testing
#-------------------------------------------------------

def FitTest(scorer,searchAlgorithm,variableProxy_arr,maxIter = 10):
	"""
	Test for 1D search algorithms.
	"""
	
	#---- max number of iteration is maxIter
	solverStopper = SolveStopperFactory.maxIterationStopper(maxIter)
	
	#---- max time for search
	#max_time = 0.01
	#solverStopper = SolveStopperFactory.maxTimeStopper(max_time)
	
	#max_accuracy = 0.001
	#solverStopper = SolveStopperFactory.maxAccuracyStopper(max_accuracy)
	
	#---- Combo Stopper test
	#comboStopper = SolveStopperFactory.comboStopper()
	#comboStopper.addStopper(solverStopper)
	#solverStopper = comboStopper
	
	solver = Solver()
	solver.setAlgorithm(searchAlgorithm)
	solver.setStopper(solverStopper)

	trialPoint = TrialPoint()
	for variableProxy in variableProxy_arr:
		trialPoint.addVariableProxy(variableProxy)

	solver.solve(scorer,trialPoint)	
	
	#---- the fitting process ended, now about results
	print "=============================================================="
	print "??????????????????????????????????????????????????????????????"
	search_alg_name = searchAlgorithm.getName()
	print "==============Fitting results for algorithm: ",search_alg_name
	solver.getScoreboard().printScoreBoard()
	(pos_min_arr,value_min) = scorer.getAnswer()
	print "===== best score ========== exact answer (pos. min, min value)=",(pos_min_arr,value_min)
	bestScore = solver.getScoreboard().getBestScore()	
	print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
	trialPoint = solver.getScoreboard().getBestTrialPoint()
	print trialPoint.textDesciption()	
	print "============STOP TEST================================="
	success_limit = 0.01
	if(abs(bestScore - value_min) > success_limit):
		print "Search failed! Stop!"
		sys.exit(1)

#--------------------------------------------
#    Main part of the tests
#--------------------------------------------
	
scorer = SquareValueScorer1D()
searchAlgorithm = BisectionSearchAlgorithm()

#---- for Bisection serach we set (name,initial value, step)
#---- and limits for serach. The step parameter will not be used.
variableProxy = VariableProxy("(x-2)^2+10",0.,0.3333)
variableProxy.setLowerLimit(- 5.0)
variableProxy.setUpperLimit(+ 10.0)

FitTest(scorer,searchAlgorithm,[variableProxy,])

#-----------------------------------

scorer = SquareValueScorer1D()
searchAlgorithm = GoldenSectionSearchAlgorithm()

#---- for Golden Section serach we set (name,initial value, step)
#---- and limits for serach. The step parameter will not be used.
variableProxy = VariableProxy("(x-2)^2+10",0.,0.3333)
variableProxy.setLowerLimit(- 5.0)
variableProxy.setUpperLimit(+ 10.0)

FitTest(scorer,searchAlgorithm,[variableProxy,])

#-----------------------------------

scorer = SquareValueScorer1D()
searchAlgorithm = SimplexSearchAlgorithm()
#searchAlgorithm = RandomSearchAlgorithm()

#---- for Simplex serach we set (name,initial value, initial search step)
#---- The *math.sqrt(2.) introduced to avoid exact solutions
variableProxy = VariableProxy("(x-2)^2+10",0.,0.5*math.sqrt(2.))

#---- Number of iterations for SimplexSearchAlgorithm is 20
#---- but for RandomSearchAlgorithm it should be 50
FitTest(scorer,searchAlgorithm,[variableProxy,],20)

#-----------------------------------

scorer = SquareValueScorer3D()
searchAlgorithm = SimplexSearchAlgorithm()
#searchAlgorithm = RandomSearchAlgorithm()

#---- for Simplex serach we set (name,initial value, initial search step)
#---- The *math.sqrt(2.) introduced to avoid exact solutions
variableProxy0 = VariableProxy("one",1.,0.2*math.sqrt(2.))
variableProxy1 = VariableProxy("two",2.,0.2*math.sqrt(2.))
#---- test that we can switch off a varibale from the minimization process
#variableProxy1.setUseInSolver(False)
variableProxy2 = VariableProxy("three",3.,0.2*math.sqrt(2.))
variableProxy_arr = []
variableProxy_arr.append(variableProxy0)
variableProxy_arr.append(variableProxy1)
variableProxy_arr.append(variableProxy2)

#---- Number of iterations for SimplexSearchAlgorithm is 50
#---- but for RandomSearchAlgorithm it should be 200
FitTest(scorer,searchAlgorithm,variableProxy_arr, 50)


sys.exit(0)
