import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccActionsConatainer:
	""" The class of the container of accelerator actions. """
	def __init__(self, name = "no name container" ):
		"""
		Creates an empty accelerator actions container.
		"""
		self.AccActionsConatainer = orbit.lattice.AccActionsConatainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.	AccLattice = 	orbit.lattice.AccLattice
		self.__name = name
		self.__type = "action container"
		self.__entranceActions = []
		self.__bodyActions = []
		self.__exitActions = []
		self.__shouldStop = False

	def performEntranceActions(self, paramDict):
		"""
		Method. It performs the set of actions at the entrance of accelerator node.
		"""
		for action in self.__entranceActions:
			if(self.__shouldStop):
				return
			action(paramDict)

	def performExitActions(self, paramDict):
		"""
		Method. It performs the set of actions
		at the exit of accelerator node.
		"""
		for action in self.__exitActions:
			if(self.__shouldStop):
				return
			action(paramDict)

	def performActions(self, paramDict):
		"""
		Method. It performs the set of actions for the
		body of accelerator node.
		"""
		for action in self.__bodyActions:
			if(self.__shouldStop):
				return
			action(paramDict)

	def addEntranceAction(self, action):
		"""
		Method. It adds an entarnce action to
		the beginning of the container.
		"""
		if(self.__entranceActions.count(action) == 0):
			self.__entranceActions.insert(0,action)

	def getNumbOfEntranceActions(self):
		"""
		Method. It returns the number of entrance actions
		in the container.
		"""
		return len(self.__entranceActions)

	def getEntranceActions(self):
		"""
		Method. It returns the entrance actions as a list.
		"""
		return self.__entranceActions

	def removeEntranceAction(self,action):
		"""
		Method. It removes the action from the container.
		"""
		return self.__entranceActions.remove(action)

	def addExitAction(self, action):
		"""
		Method. It adds an exit action to
		the beginning of the container.
		"""
		if(self.__exitActions.count(action) == 0):
			self.__exitActions.insert(0,action)

	def getNumbOfExitActions(self):
		"""
		Method. It returns the number of exit actions
		in the container.
		"""
		return len(self.__exitActions)

	def getExitActions(self):
		"""
		Method. It returns the exit actions as a list.
		"""
		return self.__exitActions

	def removeExitAction(self,action):
		"""
		Method. It removes the action from the container.
		"""
		return self.__exitActions.remove(action)

	def addAction(self, action):
		"""
		Method. It adds an action to the beginning
		of the container.
		"""
		if(self.__bodyActions.count(action) == 0):
			self.__bodyActions.insert(0,action)

	def getNumbOfActions(self):
		"""
		Method. It returns the number of  actions
		in the container.
		"""
		return len(self.__bodyActions)

	def getActions(self):
		"""
		Method. It returns the actions as a list.
		"""
		return self.__bodyActions

	def removeAction(self,action):
		"""
		Method. It removes the action from the container.
		"""
		return self.__bodyActions.remove(action)

	def setShouldStop(self, shouldStop = True):
		"""
		Method. It sets a should_stop property of
		the actions container.
		"""
		self.__shouldStop = shouldStop

	def getShouldStop(self):
		"""
		Method. It returns a should_stop property
		of the actions container.
		"""
		return self.__shouldStop

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the actions container.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the actions container.
		"""
		return self.__name
