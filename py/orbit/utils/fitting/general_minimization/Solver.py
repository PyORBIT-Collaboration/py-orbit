"""
This is a collection of classes for a general optimization problem.
The Scorer function should be minimized by changing the input parameters
that are inside the trial point class instance. The user can define 
different search algorithms that are suitable for the particular problem.
Also there are set of solver stoppers that will stop optimization process
immediately (but stopper cannot interrupt processes inside the Scorer).
"""

import os
import math
import sys
import time

from orbit.utils import NamedObject, ParamsDictObject

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
	TrialPoint
	"""
	def __init__(self):
		"""
		Constructor of the main class of the general fitting package.
		"""
		self.scoreboard = Scoreboard(self)
		self.search_algorithm = None
		self.stopper = SolveStopperFactory.runForeverStopper()
		self.scorer = None
		self.is_running = False
		
	def getScoreboard(self):
		"""
		This method returns the Scoreboard.
		"""
		return self.scoreboard	
		
	def _setScorer(self, scorer):
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
		This method returns the solver stopper class instance.
		"""
		return self.stopper

	def isRunning(self):
		"""
		This method returns true or false.
		"""	
		return self.is_running 

	def solve(self,scorer,initTrialPoint):
		"""
		This method applays the fitting algorithms to the problem.
		"""
		self.is_running = True
		
		self.scoreboard.init()
		
		self._setScorer(scorer)
		
		self.search_algorithm.setSolver(self)
		res = self.search_algorithm.setTrialPoint(initTrialPoint)
		if(not res):
			msg  = "============ Solver class: method solve(...)=============="
			msg += os.linesep		
			msg += "Cannot initialize the search algorithm"
			msg += os.linesep			
			msg += "Stop."
			msg += os.linesep
			orbitFinalize(msg)
		
		while(not self.stopper.getShouldStop()):
			self.search_algorithm.makeStep()
		
		self.is_running = False
		self._setScorer(None)

#-----------------------------------------------------------
#   Class TrialPoint
#-----------------------------------------------------------


class TrialPoint:
	"""
	This a container class for VariableProxy instances. It keeps the information about 
	the values of the parameters inside VariableProxy instances. It keeps the references
	in two forms - dictionary and array to facilitate usage in different types of 
	score functions.
	"""
	def __init__(self):
		"""
		The constructor of empty TrialPoint container.
		"""
		self._varProxy_arr = []
		self._varProxy_dict= {}
		
	def clean(self):
		"""
		"""
		self._varProxy_arr = []
		self._varProxy_dict= {}		
		
	def addVariableProxy(self,variableProxy):
		if(self._varProxy_dict.has_key(variableProxy.getName())):
			msg  = "============ TrialPoint class. Method addVariableProxy(...)=============="
			msg += os.linesep		
			msg += "============ WARNING  START=============="
			msg += os.linesep			
			msg += "You added two VariableProxy instances with the same name to TrialPoint."
			msg += os.linesep
			msg += "That is dangerous!"
			msg += os.linesep
			msg += self.textDesciption()
			msg += "============ WARNING  STOP=============="
			msg += os.linesep	
			print "msg"
		self._varProxy_arr.append(variableProxy)
		self._varProxy_dict[variableProxy.getName()] = variableProxy
		
	def addVariableProxyArr(self,variableProxy_arr):
		"""
		Adds VariableProxy instance to the container.
		"""
		for variableProxy in variableProxy_arr:
			self.addVariableProxy(variableProxy)
		
	def getCopy(self):
		"""
		Returns the copy of this istance of TrialPoint.
		"""
		trialPoint_new = TrialPoint()
		for variableProxy in self._varProxy_arr:
			variableProxy_new = VariableProxy(variableProxy)
			trialPoint_new.addVariableProxy(variableProxy_new)
		return trialPoint_new
		
	def getVariableProxyArr(self):
		"""
		Returns the reference to the inner array of VariableProxy instances.
		"""
		return self._varProxy_arr
		
	def getVariableProxyValuesArr(self):
		"""
		Returns the reference to an unbound array of values from VariableProxy variables.
		This is a convinience method.
		"""
		values_arr = []
		for variableProxy in self._varProxy_arr:
			values_arr.append(variableProxy.getValue())
		return values_arr
		
	def getVariableProxyDict(self):
		"""
		Returns the reference to the inner dictionary of VariableProxy instances.
		"""
		return self._varProxy_dict
		
	def getVariablesUsedInOptArr(self):
		"""
		Returns the reference to the unbound array of VariableProxy instances that
		should be used in the optimization process (using variableProxy.getUseInSolver()).
		This is a convinience method.
		"""
		val_arr = []
		for variableProxy in self._varProxy_arr:
			if(variableProxy.getUseInSolver()):
				val_arr.append(variableProxy.getValue())
		return val_arr
		
	def getStepsUsedInOptArr(self):
		"""
		Returns the reference to the unbound array of steps for variables that
		should be used in the optimization process (using variableProxy.getUseInSolver()).
		"""
		step_arr = []
		for variableProxy in self._varProxy_arr:
			if(variableProxy.getUseInSolver()):
				step_arr.append(variableProxy.getStep())
		return step_arr	
		
	def setVariablesUsedInOptArr(self,val_arr):
		"""
		Sets the values to the VariableProxy instances that are used 
		in the optimization process (using variableProxy.getUseInSolver()).
		"""
		nUsedInOptVars = len(self.getVariablesUsedInOptArr())
		if(len(val_arr) != nUsedInOptVars):
			msg  = "============ TrialPoint class. Method setVariablesUsedInOptArr(...)=============="
			msg += os.linesep		
			msg += "============ WARNING  START=============="
			msg += os.linesep			
			msg += "The number of variables is different from the number of VariableProxies."
			msg += os.linesep
			msg += "That is wrong! Stop"
			msg += os.linesep
			msg += str(val_arr)
			msg += os.linesep
			msg += "n used in optimization variables ="+str(nUsedInOptVars)
			msg += os.linesep			
			st = self.textDesciption()
			msg += st
			msg += "============ WARNING  STOP=============="
			msg += os.linesep	
			print msg
			return False
		#------------------------------------------------
		count = 0
		for variableProxy in self._varProxy_arr:
			if(variableProxy.getUseInSolver()):
				variableProxy.setValue(val_arr[count])
				count += 1

	def setStepsUsedInOptArr(self,step_arr):
		"""
		Sets the of steps for variables that should be used in 
		the optimization process (using variableProxy.getUseInSolver()).
		"""
		nUsedInOptVars = len(self.getVariablesUsedInOptArr())
		if(len(step_arr) != nUsedInOptVars):
			msg  = "============ TrialPoint class. Method setStepsUsedInOptArr(...)=============="
			msg += os.linesep		
			msg += "============ WARNING  START=============="
			msg += os.linesep			
			msg += "The number of steps variables is different from the number of VariableProxies."
			msg += os.linesep
			msg += "That is wrong! Stop"
			msg += os.linesep
			msg += str(step_arr)
			msg += os.linesep
			msg += "n used in optimization variables ="+str(nUsedInOptVars)
			msg += os.linesep			
			st = self.textDesciption()
			msg += st
			msg += "============ WARNING  STOP=============="
			msg += os.linesep	
			print msg
			return False
		#------------------------------------------------
		count = 0
		for variableProxy in self._varProxy_arr:
			if(variableProxy.getUseInSolver()):
				variableProxy.setStep(step_arr[count])
				count += 1

	def isAcceptable(self):
		"""
		It returns True of False if all variableProxies are inside value limits
		"""
		for variableProxy in self._varProxy_arr:
			if(not variableProxy.isAcceptable()):
				return False
		return True

	def textDesciption(self):
		"""
		Returns the text description of the VariableProxy-s inside this TrialPoint.
		"""
		st = "======== TrialPoint ==========="
		st = st + os.linesep
		st = st + " Name                       Value          Step       Use      Limit_Min       Limit_Max  "
		for variableProxy in self._varProxy_arr:
			st += os.linesep
			st += "%10s "%variableProxy.getName()
			st += "  %14.7g  %14.7g  "%(variableProxy.getValue(),variableProxy.getStep())
			st += "  %1d  "%variableProxy.getUseInSolver()
			st += "  %14.7g  %14.7g  "%(variableProxy.getLowerLimit(),variableProxy.getUpperLimit())
			
		return st

#-----------------------------------------------------------
#  Class SolveStopper abstract class 
#-----------------------------------------------------------
			
class SolveStopper:
	"""
	The SolveStopper defines if we should stop solver's optimization process
	because of some condition.
	"""
	def __init__(self):
		"""
		Constructor of stopper that always returns Fasle 
		as the answer to "Should I Stop?" question.
		"""
		self._shouldStop = False
		
	def getShouldStop(self):
		"""
		Returns True of False as the answer to "Should I Stop?" question.
		"""
		return self._shouldStop
		
	def setShouldStop(self,shouldStop):
		"""
		Sets True of False as the answer to "Should I Stop?" question.
		"""
		self._shouldStop = shouldStop
		
	def checkStopConditions(self,solver):
		"""
		This is an abstract method. The derived classes should implement this method.
		Here you can get scoreBoard from solver to evaluate the stack with best 
		scores and TrialPoints.
		Inside this method you should put the answer to "Should I Stop?" question.
		"""
		scoreBoard = solver.getScoreboard()
		(score,trialPoint) = scoreBoard.getCurrentScoreAndTrialPoint()
		scoresHistoryStack = scoreBoard.getHistoryStack()
		#---- do something in your stopper 
		#---- if you want to stop optimization use setShouldStop(True) 
		#---- method.
			
#-----------------------------------------------------------			
#  Class ForeverStopper			
#-----------------------------------------------------------
class RunForeverStopper(SolveStopper):
	"""
	The stopper implementation that always says False 
	to "Should I Stop?" question.
	"""
	def __init__(self):
		SolveStopper.__init__(self)
		self.setShouldStop(False)

#-----------------------------------------------------------			
#  Class Max Iterations Stopper		
#-----------------------------------------------------------
class MaxIterationStopper(SolveStopper):
	"""
	The stopper implementation that will say True 
	to "Should I Stop?" question if the number of Scorer evaluation is more then
	the maximal allowable iterations.
	"""
	def __init__(self,max_iteration):
		SolveStopper.__init__(self)
		self.max_iteration = max_iteration
		self.setShouldStop(False)
		
	def checkStopConditions(self,solver):
		"""
		Implementation of the abstract method of the parent class.
		"""
		iteration = solver.getScoreboard().getIteration()
		if(iteration >= self.max_iteration):
			self.setShouldStop(True)

#-----------------------------------------------------------			
#  Class Max Time Stopper			
#-----------------------------------------------------------
class MaxTimeStopper(SolveStopper):
	"""
	The stopper implementation that will say True 
	to "Should I Stop?" question if the time is up.
	"""
	def __init__(self,max_time):
		SolveStopper.__init__(self)
		self.max_time = max_time
		self.setShouldStop(False)
		
	def checkStopConditions(self,solver):
		"""
		Implementation of the abstract method of the parent class.
		"""
		tm = solver.getScoreboard().getRunTime()
		if(tm >= self.max_time):
			self.setShouldStop(True)

#-----------------------------------------------------------			
#  Class Max Accuracy Stopper			
#-----------------------------------------------------------
class MaxAccuracyStopper(SolveStopper):
	"""
	The stopper implementation that will say True 
	to "Should I Stop?" question if the difference between 
	two last best scores is less than the max accuracy.
	"""
	def __init__(self,max_accuracy):
		SolveStopper.__init__(self)
		self.max_accuracy = max_accuracy
		self.setShouldStop(False)
			
	def checkStopConditions(self,solver):
		"""
		Implementation of the abstract method of the parent class.
		"""
		scoresHistoryStack = solver.getScoreboard().getHistoryStack()
		n_scores = len(scoresHistoryStack)
		if(n_scores < 2): return
		score1 = scoresHistoryStack[n_scores-2][0]
		score2 = scoresHistoryStack[n_scores-1][0]
		if(abs(score1 - score2) < self.max_accuracy):
			self.setShouldStop(True)

#-----------------------------------------------------------			
#  Class Combo Stopper			
#-----------------------------------------------------------
class ComboStopper(SolveStopper):
	"""
	The stopper implementation that will say True to "Should I Stop?" 
	question if at least one of stoppers registered in it will say True.
	"""
	def __init__(self):
		SolveStopper.__init__(self)
		self.stopper_arr = []
		
	def addStopper(self,stopper):
		if(not isinstance(stopper,SolveStopper)):
			return False
		self.stopper_arr.append(stopper)

	def getShouldStop(self):
		"""
		Returns True of False as the answer to "Should I Stop?" question.
		"""
		for stopper in self.stopper_arr:
			if(stopper.getShouldStop()):
				return True
		return False
		
	def setShouldStop(self,shouldStop):
		"""
		Sets True of False as the answer to "Should I Stop?" question.
		"""
		for stopper in self.stopper_arr:
			stopper.setShouldStop(shouldStop)

	def checkStopConditions(self,solver):
		"""
		Implementation of the abstract method of the parent class.
		"""
		for stopper in self.stopper_arr:
			stopper.checkStopConditions(solver)

#-----------------------------------------------------------
#  Class SolveStopperFactory
#-----------------------------------------------------------

class SolveStopperFactory:
	"""
	The Factory for stoppers.
	"""
	@classmethod
	def runForeverStopper(self):
		"""
		Returns the Run Forever stopper.
		"""
		return RunForeverStopper()
		
	@classmethod
	def maxIterationStopper(self,max_iteration):
		"""
		Returns the Max Iterations stopper.
		"""
		return MaxIterationStopper(max_iteration)
		
	@classmethod
	def maxTimeStopper(self,max_time):
		"""
		Returns the Max Time stopper.
		"""
		return MaxTimeStopper(max_time)
		
	@classmethod
	def maxAccuracyStopper(self,max_accuracy):
		"""
		Returns the Max Accuracy stopper.
		"""
		return MaxAccuracyStopper(max_accuracy)
		
	@classmethod
	def comboStopper(self):
		"""
		Returns the Combo stopper.
		By itself this stopper is not operational.
		You have to add stoppers into this combo
		stopper.
		"""
		return ComboStopper()
		
#-----------------------------------------------------------
#  ScoreBoard Action listeners 
#-----------------------------------------------------------
class ScoreboardActionListener:
	"""
	This is an abstract class for actions inside Scoreboard.
	"""
	def __init__(self):
		pass
	
	def performAction(self,solver):
		"""
		Perform necessary action. This method should be implemented 
		in the children classes.
		"""
		pass

#-----------------------------------------------------------
#  Class Scoreboard
#----------------------------------------------------------	
	
class Scoreboard:
	"""
	Scoreboard class keeps the trace of all best scores 
	(as tuple (score,TrialPoint)) in the Scores History Stack.
	The maximal size of the stack is 100 by default.
	The user can define his/her own stack size.
	The best score&Trialpoint are returned by  getBestScore() and 
	getBestTrialPoint() methods.
	"""
	def __init__(self,solver):
		self.solver = solver
		self.init()
		#---- listeners
		self.newTrialPointListener_arr = []
		self.bestScoreListener_arr = []
		
	def init(self):
		"""
		Sets the clean state of the Scoreboard.
		"""
		self.start_time = time.time()
		self.run_time = 0.
		self.iterations = 0
		self.bestScore = 0.00000001*sys.float_info.max
		self.currentScore = self.bestScore
		self.currentTrialPoint = None
		self.bestTrialPoint = None
		#------------------------------
		self.historyStackSize = 100
		#---- self.scoresHistoryStack = [(score,TrialPoint),...]
		self.scoresHistoryStack = []
		
	def addScoreTrialPoint(self,score,trialPoint):
		"""
		Adds the new score and Trial Point combination. If this score is the best it 
		will add them to the best scores & trial points history stack.
		"""
		self.currentScore = score
		self.currentTrialPoint = trialPoint
		self.run_time = time.time() - self.start_time
		self.iterations += 1
		#---- call performAction() method of a new trial point listeners
		for listener in self.newTrialPointListener_arr:
			listener.performAction(self.solver)
		#---------------------------------------------------------------
		self.solver.getStopper().checkStopConditions(self.solver)
		if(score <= self.bestScore):
			self.bestTrialPoint = trialPoint.getCopy()
			self.bestScore = score
			if(len(self.scoresHistoryStack) >= self.historyStackSize):
				self.scoresHistoryStack.pop()
			self.scoresHistoryStack.append((score,self.iterations,self.bestTrialPoint.getCopy()))
			#---- call performAction() method of a new best score listeners
			for listener in self.bestScoreListener_arr:
				listener.performAction(self.solver)
			#---------------------------------------------------------------
		
	def getIteration(self):
		"""
		Retuns the number of iterations so far.
		"""
		return self.iterations
		
	def setHistoryStackSize(self,historyStackSize):
		"""
		Changes the size of the hystory stack.
		"""
		self.historyStackSize = historyStackSize
		self.scoresHistoryStack = self.scoresHistoryStack[:self.historyStackSize]
		
	def getHistoryStackSize(self):
		"""
		Returns the size of the hystory stack.
		"""
		return self.historyStackSize
		
	def getHistoryStack(self):
		"""
		Returns the hystory stack with (score,trial_point) data.
		"""
		return self.scoresHistoryStack
		
	def setRunTime(self):
		"""
		Method sets the run-time.
		"""
		self.run_time = time.time() - self.start_time
		
	def getRunTime(self):
		"""
		Returns the run-time.
		"""
		return self.run_time
		
	def getBestTrialPoint(self):
		"""
		Returns the unbound best trial point.
		"""
		return self.bestTrialPoint.getCopy()
		
	def getBestTrialPointReference(self):
		"""
		Returns the reference (not copy) to the best trial point.
		"""
		return self.bestTrialPoint		
		
	def getBestScore(self):
		"""
		Returns the best score.
		"""
		return self.bestScore
		
	def getCurrentScoreAndTrialPoint(self):
		"""
		Method is used by stopper to check stop conditions.
		"""
		return (self.currentScore,self.currentTrialPoint)
		
	def printScoreBoard(self):
		"""
		Prints the scoreboard state.
		"""
		print "==================== Score Board Stack ======START==========="
		for (score,iteration,trialPoint) in self.scoresHistoryStack:
			st = "===== score = %15.8g "%score+ "   iter.="+str(iteration)+"  "
			st += trialPoint.textDesciption()
			print st
		print "==================== Score Board Stack =======STOP==========="
		
	def addNewTrialPointListener(self,newTrialPointListener):
		"""
		Adds a new trial point listener.
		"""
		self.newTrialPointListener_arr.append(newTrialPointListener)
		
	def addBestScoreListener(self,bestScoreListener):
		"""
		Adds a new best score listener.
		"""		
		self.bestScoreListener_arr.append(bestScoreListener)
			

#====================================================================
#       class VariableProxy
#====================================================================

class VariableProxy(NamedObject,ParamsDictObject):
	"""
	This class represents the parameter for the score function in the fitting process. 
	"""
	def __init__(self, *arg, **kwargs):
		"""
		The constructor of the VariableProxy class should have the following signatures
		VariableProxy(parameterProxy_in)
		VariableProxy(name, value, step)
		VariableProxy(name = ""???"", value = ???, step = ???)
		"""
		ParamsDictObject.__init__(self)
		self.setName("unknown")
		self.value = 0.
		self.step = 0.
		self.useInSolver = True
		
		self.lowerLimit = - 0.00000001*sys.float_info.max
		self.upperLimit = + 0.00000001*sys.float_info.max
		
		#print "debug len(arg)=",len(arg)," arg=",arg
		#print "debug len(kwargs)=",len(kwargs)," kwargs=",kwargs
		
		if(not (len(arg) != 1 and len(arg) != 2 and len(arg) != 3)):	
			if(len(arg) == 1):
				varProxy = arg[0]
				if(isinstance(varProxy,VariableProxy)):
					self.setName(varProxy.getName())
					self.value = varProxy.getValue()
					self.step = varProxy.getStep()
					self.useInSolver = varProxy.getUseInSolver()
					self.lowerLimit = varProxy.getLowerLimit()
					self.upperLimit = varProxy.getUpperLimit()
					self.updateParamsDict(varProxy.getParamsDict())
				else:
					msg = "VariableProxy constructor. If there is only one argument it should be only VariableProxy."
					msg = msg + os.linesep
					msg = "This argument is not!"
					msg = msg + os.linesep			
					msg = msg + "Stop."
					msg = msg + os.linesep
					orbitFinalize(msg)
			elif(len(arg) == 2):
				self.setName(arg[0])
				self.value = arg[1]
				self.step = 0.			
			else:
				self.setName(arg[0])
				self.value = arg[1]
				self.step = arg[2]
		elif(len(arg) == 0 and len(kwargs) == 3 and kwargs.has_key("name") and kwargs.has_key("value") and  kwargs.has_key("step")):
			self.setName(kwargs["name"])
			self.value = kwargs["value"]
			self.step = kwargs["step"]
		else:
			msg  = "VariableProxy constructor. It should have one of the forms:"
			msg += os.linesep
			msg += "1. Copy constructor: VariableProxy(parameterProxy_in)"
			msg += os.linesep
			msg += "2. VariableProxy(name, value)"
			msg += os.linesep			
			msg += "3. VariableProxy(name, value, step)"
			msg += os.linesep			
			msg += "4. VariableProxy(name = ""???"", value = ???, step = ???)"
			msg += os.linesep			
			msg += "Stop."
			msg += os.linesep
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
		
	def isAcceptable(self):
		"""
		It returns True of False if the self.value is inside value limits
		"""
		if(self.value < self.lowerLimit): return False
		if(self.value > self.upperLimit): return False
		return True
		
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
		
	def setUseInSolver(self,useInSolver):
		"""
		Set the Boolean variable that defines if the VariableProxy 
		will be used in Solver for optimization.
		"""
		self.useInSolver = useInSolver
		
	def getUseInSolver(self):
		"""
		Returns the Boolean variable that defines if the VariableProxy 
		will be used in Solver for optimization.
		"""		
		return self.useInSolver
	
#====================================================================
#     Class Scorer 
#====================================================================

class Scorer:
	"""
	This class calculates the score for the trial point instance.
	The score should be minimized by Solver with the help from
	Search Algorithm. Time to end the minimization process is defined
	by the Solve Stopper.
	"""
	def __init__(self):
		pass
		
	def getScore(self,trialPoint):
		"""
		This is an empty method. The child classes should implement this method.
		"""
		pass
		
#====================================================================
#     Class SearchAgorithm
#====================================================================		

class SearchAgorithm(NamedObject):
	"""
	The SearchAgorithm is an abstract class. A concrete implementation is 
	needed to provide functionality of the search algorithm.
	"""
	def __init__(self):
		self.solver = None
		self.initTrialPoint = None
		self.setName("No Algorithm")
	
	def setSolver(self,solver):
		"""
		Sets the solver instance for the search algorithm.
		"""
		self.solver = solver
		
	def setTrialPoint(self,initTrialPoint):
		"""
		This derived class should add functionality to this method as needed. 
		"""
		self.initTrialPoint = initTrialPoint.getCopy()
		res = self.init()
		return res
		
	def init(self):
		if(self.initTrialPoint == None or self.solver == None): return False
		return True
	
	def makeStep(self):
		"""
		This is an abstarct method. The derived classes should implement this method.
		"""
		pass
