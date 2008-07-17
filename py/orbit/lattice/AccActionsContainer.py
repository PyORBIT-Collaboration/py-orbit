import sys
import os

from orbit.utils import orbitFinalize
from orbit.utils import NamedObject
import orbit

class AccActionsContainer(NamedObject):
	"""
	Class. Container for accelerator actions.
	"""

	ENTRANCE = 0
	BODY     = 1
	EXIT     = 2

	BEFORE   = 0
	AFTER    = 1

	def __init__(self, name = "no name container"):
		"""
		Constructor. Creates empty accelerator actions container.
		"""
		NamedObject.__init__(self, name)
		#Array with entrance, body, and exit actions.
		self.__actionsArr = ([],[],[])

	def addAction(self, action, place):
		"""
		Method. Appends an action at the entrance, body, or exit
		of an accelerator node. 
		"""
		self.__actionsArr[place].append(action)

	def removeAction(self, action, place):
		"""
		Method. Removes an action from the entrance, body, or exit
		of an accelerator node. 
		"""
		self.__actionsArr[place].remove(action)

	def getActions(self, place):
		"""
		Method. Returns a list of the actions at the entrance, body,
		or exit of an accelerator node.
		"""
		return self.__actionsArr[place]

	def performActions(self, paramsDict, place):
		"""
		Method. Performs actions at the entrance, body, or exit
		of an accelerator node.
		"""
		for action in self.__actionsArr[place]:
			action(paramsDict)
