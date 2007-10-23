import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccActionsContainer:
	"""
	Class. Contains the accelerator actions.
	"""
	def __init__(self, name = "no name container" ):
		"""
		Constructor. Creates empty accelerator actions container.
		"""
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.AccLattice = orbit.lattice.AccLattice
		self.__name = name
		self.__type = "action container"
		self.__entranceActions = []
		self.__bodyActions = []
		self.__exitActions = []
		self.__shouldStop = False

	def setName(self, name = "no name"):
		"""
		Method. Sets the name of the actions container.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. Returns the name of the actions container.
		"""
		return self.__name

	def appendEntranceAction(self, action):
		"""
		Method. Append an action at the entrance
		of the accelerator node.
		"""
		if(self.__entranceActions.count(action) == 0):
			self.__entranceActions.append(action)

	def removeEntranceAction(self, action):
		"""
		Method. Remove an action from the entrance
		of the accelerator node.
		"""
		return self.__entranceActions.remove(action)

	def getNumbOfEntranceActions(self):
		"""
		Method. Return the number of actions at the entrance
		of the accelerator node.
		"""
		return len(self.__entranceActions)

	def getEntranceActions(self):
		"""
		Method. Returns a list of the actions at the entrance
		of the accelerator node
		"""
		return self.__entranceActions

	def performEntranceActions(self, paramsDict):
		"""
		Method. Perform the required actions at the entrance
		of the accelerator node.
		"""
		for action in self.__entranceActions:
			if(self.__shouldStop):
				return
			action(paramsDict)

	def appendBodyAction(self, action):
		"""
		Method. Append an action in the body
		of the accelerator node.
		"""
		if(self.__bodyActions.count(action) == 0):
			self.__bodyActions.append(action)

	def removeBodyAction(self, action):
		"""
		Method. Remove an action from the body
		of the accelerator node.
		"""
		return self.__bodyActions.remove(action)

	def getNumbOfBodyActions(self):
		"""
		Method. Return the number of actions in the body
		of the accelerator node.
		"""
		return len(self.__bodyActions)

	def getBodyActions(self):
		"""
		Method. Returns a list of the actions in the body
		of the accelerator node
		"""
		return self.__bodyActions

	def performBodyActions(self, paramsDict):
		"""
		Method. Perform the required actions in the body
		of the accelerator node.
		"""
		for action in self.__bodyActions:
			if(self.__shouldStop):
				return
			action(paramsDict)

	def appendExitAction(self, action):
		"""
		Method. Append an action at the exit
		of the accelerator node.
		"""
		if(self.__exitActions.count(action) == 0):
			self.__exitActions.append(action)

	def removeExitAction(self, action):
		"""
		Method. Remove an action from the exit
		of the accelerator node.
		"""
		return self.__exitActions.remove(action)

	def getNumbOfExitActions(self):
		"""
		Method. Return the number of actions at the exit
		of the accelerator node.
		"""
		return len(self.__exitActions)

	def getExitActions(self):
		"""
		Method. Returns a list of the actions at the exit
		of the accelerator node
		"""
		return self.__exitActions

	def performExitActions(self, paramsDict):
		"""
		Method. Perform the required actions at the exit
		of the accelerator node.
		"""
		for action in self.__exitActions:
			if(self.__shouldStop):
				return
			action(paramsDict)

	def setShouldStop(self, shouldStop = True):
		"""
		Method. Sets shouldStop to True in
		the actions container.
		"""
		self.__shouldStop = shouldStop

	def getShouldStop(self):
		"""
		Method. Returns the value of shouldStop in
		the actions container.
		"""
		return self.__shouldStop

