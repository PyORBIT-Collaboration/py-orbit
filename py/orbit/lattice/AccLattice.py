import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccLattice:
	"""
	Class. The accelerator lattice class.
	A lattice contains elements, or child nodes.
	"""
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates an empty accelerator lattice.
		"""
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLattice = orbit.lattice.AccLattice
		self.__name = name
		self.__type = "lattice"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def setName(self, name = "no name"):
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

	def setLength(self, L = 0.):
		"""
		Method. Sets the physical length of the lattice.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. Returns the physical length of the lattice.
		"""
		return self.__length

	def setInitialized(self, initialized = True):
		"""
		Method. Sets the initialization status (True or False).
		"""
		actions = self.AccActionsContainer()

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			node.setInitialized(initialized)

		actions.appendExitAction(accElemExit)
		self.trackActions(actions, paramsDict)
		self.__isInitialized = initialized

	def isInitialized(self):
		"""
		Method. Returns the initialization status (True or False).
		"""
		return self.__isInitialized

	def initialize(self, actions = None, paramsDict = {}):
		"""
		Method. Initializes the lattice and child node structures.
		"""
		if(actions == None):
			actions = self.AccActionsContainer()
		d = {"position":0}
		d["position"] = 0.

		def accElemEntrance(paramsDict):
			node = paramsDict["node"]
			node.initialize(paramsDict)

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			d["position"] += node.getLength()

		actions.appendEntranceAction(accElemEntrance)
		actions.appendExitAction(accElemExit)
		self.trackActions(actions, paramsDict)
		self.__length = d["position"]
		self.__isInitialized = True

	def insertChildNode(self, node, index = 0):
		"""
		Method. Inserts a child node into the lattice.
		The third parameter is the child node index.
		"""
		if(isinstance(node,self.AccElement) != True):
			msg = "A child of an AccLattice must be an AccElement!"
			msg = msg + os.linesep
			msg = msg + "method insertChildNode(self, node, index)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		self.__children.insert(index, node)

	def appendChildNode(self, node):
		"""
		Method. Appends a child node to the lattice.
		"""
		if(isinstance(node, self.AccElement) != True):
			msg = "A child of an AccLattice must be an AccElement!"
			msg = msg + os.linesep
			msg = msg + "method appendChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of lattice = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def removeChildNode(self, index = 0):
		"""
		Method. Removes a child node with given index from the lattice.
		"""
		if((index < 0) or (index >= (len(self.__children)))):
			msg = "Child node index is out of range!"
			msg = msg + os.linesep
			msg = "Range, index = 0 - " + str(len(self.__children)-1) + ", " + str(index)
			msg = msg + os.linesep
			msg = msg + "method removeChildNode(self, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			orbitFinalize(msg)
		del self.__children[index:(index+1)]

	def removeAllChildNodes(self):
		"""
		Method. Removes all children from the lattice.
		"""
		self.__children = []

	def getAllChildNodes(self):
		"""
		Method. Returns a list of all children in the lattice.
		"""
		return self.__children

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

		actions.appendBodyAction(track)
		self.trackActions(actions, paramsDict)
