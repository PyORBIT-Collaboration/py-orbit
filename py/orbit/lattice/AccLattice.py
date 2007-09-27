import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccLattice:
	"""
	The accelerator lattice class.
	"""
	
	def __init__(self, name = "no name"):
		"""
		Method. Creates an empty accelerator lattice.
		"""
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.AccLattice = orbit.lattice.AccLattice
		self.__name = name
		self.__type = "lattice"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def initialize(self, actions = None, paramsDict = {}):
		"""
		Method. Initializes the lattice and child node structures.
		"""
		if(actions == None):
			actions = self.AccActionsContainer()
		d = {"position":0, "position_start_line":0}
		d["position"] = 0.
		d["position_start_line"] = []

		def accElemEntrance(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node, self.AccElement)):
				node.initialize(paramsDict)
			if(isinstance(node, self.AccLine)):
				node.setLength(0.)
				d["position_start_line"].append(d["position"])

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			d["position"] += node.getLength()
			if(isinstance(node, self.AccLine)):
				line_ind = len(d["position_start_line"]) - 1
				start_pos = d["position_start_line"][line_ind]
				del d["position_start_line"][line_ind]
				node.setLength(d["position"] - start_pos)

		actions.insertEntranceAction(accElemEntrance)
		actions.insertExitAction(accElemExit)

		self.trackActions(actions, paramsDict)
		self.__length = d["position"]
		self.__isInitialized = True

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. Tracks the actions through all nodes in the lattice.
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

	def insertChildNode(self, node):
		"""
		Method. Appends a child node to this node.
		"""
		if(isinstance(node, self.AccElement) != True and isinstance(node, self.AccLine) != True):
			msg = "A child of an AccLattice must be an AccElement or AccLine!"
			msg = msg + os.linesep
			msg = msg + "method insertChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of lattice = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def setName(self, name = "no_name"):
		"""
		Method. Sets the name of the lattice.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. Returns the name of the lattice.
		"""
		return self.__name

	def getType(self):
		"""
		Method. Returns the type of the lattice.
		"""
		return self.__type

	def getLength(self):
		"""
		Method. Returns the physical length of the lattice.
		"""
		return self.__length

	def isInitialized(self):
		"""
		Method. Returns the initialization status (True or False).
		"""
		return self.__isInitialized

	def setInitialized(self, initialized = True):
		"""
		Method. Sets the initialization status (True or False).
		"""
		actions = self.AccActionsContainer()

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			node.setInitialized(initialized)

		actions.insertExitAction(accElemExit)
		self.trackActions(actions,paramsDict)
		self.__isInitialized = initialized

	def trackBunch(self, bunch, actions = None, paramsDict = {}):
		"""
		Method. Tracks a bunch trough the lattice.
		"""
		paramsDict["bunch"] = bunch
		if(actions == None):
			actions = self.AccActionsContainer()

		def track(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node, self.AccElement)):
				node.track(paramsDict)

		actions.insertBodyAction(track)
		self.trackActions(actions, paramsDict)