import sys
import os

from orbit.utils import orbitFinalize
from orbit.utils import NamedObject
import orbit

class AccActionsContainer(NamedObject):
	"""
	Container for accelerator actions.
	"""
	
	ENTRANCE = 0
	BODY = 1
	EXIT = 2
	
	BEFORE = 0
	AFTER = 1
	
	def __init__(self, name = "no name container" ):
		"""
		Constructor creates empty accelerator actions container.
		"""
		NamedObject.__init__(self,name)
		#an array with entrance,body, and exit actions
		self.__actionsArr = ([],[],[])

	def addAction(self, action, place):
		"""
		It appends an action at the entrance,body, or exit
		of the accelerator node. 
		"""
		self.__actionsArr[place].append(action)

	def removeAction(self, action, place):
		"""
		It removes an action from the entrance,body, or exit
		of the accelerator node. 
		"""
		self.__actionsArr[place].remove(action)
	
	def getActions(self, place):
		"""
		It returns a list of the actions at the entrance,body, or exit
		of the accelerator node.
		"""
		return self.__actionsArr[place]

	def performActions(self, paramsDict, place):
		"""
		It performs the required actions at the entrance,body, or exit
		of the accelerator node.
		"""
		for action in self.__actionsArr[place]:
			action(paramsDict)

