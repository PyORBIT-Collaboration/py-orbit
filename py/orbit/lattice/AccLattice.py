import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccLattice:
	""" The class of the accelerator lattice. """
	def __init__(self,name = "no name"):
		"""
		Method. It creates an empty accelerator lattice.
		"""
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.	AccLattice = 	orbit.lattice.AccLattice
		self.__name = name
		self.__type = "lattice"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def initialize(self, actions = None, paramsDict = {}):
		"""
		Method. It initializes the necessary structures of the lattice and child nodes.
		"""
		if(actions == None):
			actions = self.AccActionsContainer()
		d = {"position":0,"position_start_line":0}
		d["position"] = 0.
		d["position_start_line"] = []

		def accElemEntrance(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node,self.AccElement)):
				node.initialize(paramsDict)
			if(isinstance(node,self.AccLine)):
				node.setLength(0.)
				d["position_start_line"].append(d["position"])

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			d["position"] += node. getLength()
			if(isinstance(node,self.AccLine)):
				line_ind = len(d["position_start_line"]) - 1
				start_pos = d["position_start_line"][line_ind]
				del d["position_start_line"][line_ind]
				node.setLength(d["position"] - start_pos)

		actions.addEntranceAction(accElemEntrance)
		actions.addExitAction(accElemExit)

		self.trackActions(actions,paramsDict)
		self.__length = d["position"]
		self.__isInitialized = True

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. It is tracking the actions through the all nodes.
		"""
		paramsDict["lattice"] = self
		paramsDict["actions"] = actionsContainer
		paramsDict["node"] = self
		paramsDict["parentNode"] = None

		actionsContainer.performEntranceActions(paramsDict)

		for node in self.__children:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return

			node.trackActions(actionsContainer, paramsDict)

			if(actionsContainer.getShouldStop()):
				return

		paramsDict["node"] = self
		paramsDict["parentNode"] = None
		actionsContainer.performExitActions(paramsDict)

	def addChildNode(self, node):
		"""
		Method. It adds the child node to the end of the this node.
		"""
		if(isinstance(node,self.AccElement) != True and isinstance(node,self.AccLine) != True):
			msg = "The child of AccLattice can be  AccElement or AccLine only!"
			msg = msg + "method addChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of lattice=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the lattice.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the lattice.
		"""
		return self.__name

	def getType(self):
		"""
		Method. It returns a type of the lattice.
		"""
		return self.__type

	def getLength(self):
		"""
		Method. It returns a physical length of the node.
		"""
		return self.__length

	def isInitialized(self):
		"""
		Method. It returns true or false about an initialization status.
		"""
		return self.__isInitialized

	def setInitialized(self, initialized = True):
		"""
		Method. Sets an initialization status.
		"""
		actions = self.AccActionsContainer()

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			node.setInitialized(initialized)

		actions.addExitAction(accElemExit)
		self.trackActions(actions,paramsDict)
		self.__isInitialized = initialized

	def trackBunch(self, bunch, actions = None, paramsDict = {}):
		"""
		Method. It track the bunch trough the lattice.
		"""
		paramsDict["bunch"] = bunch
		if(actions == None):
			actions = self.AccActionsContainer()

		def track(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node,self.AccElement)):
				node.track(paramsDict)

		actions.addAction(track)
		self.trackActions(actions,paramsDict)