import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccElement:
	"""
	Class. Base class of the accelerator elements hierarchy.
	"""
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates an empty accelerator element.
		"""
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.AccLattice = orbit.lattice.AccLattice
		#------------------------------------------------
		# there is no position parameter,
		# because the element may be in more than one lattice
		#------------------------------------------------
		self.__name = name
		self.__type = "generic"
		#------------------------------------------------
		# nParts - number of parts into which 
		#          this element is divided.
		#------------------------------------------------
		self.__nParts = 1
		self.__lengthArr = [0.]
		self.__length = 0.
		self.__activePartIndex = 0
		self.__params = {}
		self.__isInitialized = False
		#------------------------------------------------
		# sub-elements that constitute the entrance,
		# inside, and exit of this element
		# __bodyChildNodesDict - dictionary with
		#                        indexes as keys and
		#                        lists of nodes as elements
		#------------------------------------------------
		self.__entranceChildNodes = []
		self.__bodyChildNodesDict = {0:[]}
		self.__exitChildNodes = []

	def initialize(self, paramsDict):
		"""
		Method. Initializes the essential structure of the node.
		Calls the _initialize(probe) method of the appropriate
		subclass before initializing all children. The subclasses
		must implement the _initialize(probe) method.
		"""
		#------------------------------------------------		
		# initialization of self
		#------------------------------------------------
		self._initialize(paramsDict)
		self.setInitialized(True)

	def _initialize(self, paramsDict):
		"""
		Abstract Method. The subclasses must implement 
		_initialize(probe) method.
		"""
		pass

	def setName(self, name = "no_name"):
		"""
		Method. Sets the name of the node.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. Returns the name of the node.
		"""
		return self.__name

	def setType(self, tp = "generic"):
		"""
		Method. Sets the type of the node.
		"""
		self.__type = tp

	def getType(self):
		"""
		Method. Returns the type of the node.
		"""
		return self.__type

	def setnParts(self, n = 1):
		"""
		Method. Sets the number of Body parts of the node.
		This method does not call the initialize() method.
		"""
		self.__nParts = n
		self.__lengthArr = []
		self.__bodyChildNodesDict = {}
		for i in xrange(self.__nParts):
			self.__lengthArr.append(0.)
			self.__bodyChildNodesDict[i] = []

	def getnParts(self):
		"""
		Method. Returns the number of body parts of the node.
		"""
		return self.__nParts

	def setPartLength(self, index, L = 0.):
		"""
		Method. Sets the physical length of the part
		of the node with given index.
		"""
		self.__lengthArr[index] = L

	def getPartLength(self, index):
		"""
		Method. Returns the physical length of the part
		of the node with given index.
		"""
		return self.__lengthArr[index]

	def setLength(self, L = 0.):
		"""
		Method. Sets the physical length of the node.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. Returns the physical length of the node.
		"""
		return self.__length

	def getActivePartIndex(self):
		"""
		Method. Returns the active part index of the node.
		"""
		return self.__activePartIndex

	def addParam(self, key, value):
		"""
		Method. Sets the parameter of the node
		"""
		self.__params[key] = value

	def setParamsDict(self, params):
		"""
		Method. Sets an external dictionary as a parameter
		dictionary for this element.
		"""
		self.__params = params

	def addParams(self, params):
		"""
		Method. Adds the parameters from an external
		dictionary to the parameters of the node.
		The values of the existing keys' will be replaced
		by new ones from the external dictionary.
		"""
		self.__params.update(params)

	def getParam(self, key):
		"""
		Method. Returns the parameters of the node
		"""
		if(not self.hasParam(key)):
			msg = "The node does not have a prameter for the key you requested!"
			msg = msg + os.linesep
			msg = msg + "method getParam(self, key)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "key = " + str(key)
			orbitFinalize(msg)
		return self.__params[key]

	def getParamsDict(self):
		"""
		Method. Returns the whole parameters dictionary.
		"""
		return self.__params

	def hasParam(self, key):
		"""
		Method. Returns True if the node has a parameter for this key.
		Returns False otherwise.
		"""
		return self.__params.has_key(key)

	def setInitialized(self, status):
		"""
		Method. Sets the initialization status (True or False).
		"""
		if(type(status) != type(True)):
			orbitFinalize("setInitialized(boolean) needs boolean argument")
		self.__isInitialized = status

	def isInitialized(self):
		"""
		Method. Returns the initialization status (True or False).
		"""
		return self.__isInitialized

	def insertEntranceChildNode(self, node, index = -1):
		"""
		Method. Inserts a child node into the entrance to this node.
		The third parameter is the child node index.
		"""
		if(isinstance(node, self.AccElement) != True):
			msg = "A child of AccElement must be a subclass!"
			msg = msg + os.linesep
			msg = msg + "method insertEntranceChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		if(index <= 0):
			self.__entranceChildNodes.insert(0, node)
		else:
			if(index < len(self.__entranceChildNodes)):
				self.__entranceChildNodes.insert(index, node)
			else:
				self.__entranceChildNodes.append(node)

	def getEntranceChildNodes(self):
		"""
		Method. Returns all children at Entrance of this node
		as a list.
		"""
		res = []
		for node in self.__entranceChildNodes:
			res.append(node)
		return res

	def insertBodyChildNode(self, node, index = 0):
		"""
		Method. Inserts a child node into the body of this node.
		The third parameter is the child node index.
		"""
		if(isinstance(node, self.AccElement) != True):
			msg = "The child of AccElement must be a subclass!"
			msg = msg + os.linesep
			msg = msg + "method insertBodyChildNode(self, node, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		chArr = self.__bodyChildNodesDict[index]
		chArr.insert(0, node)

	def getBodyChildNodes(self):
		"""
		Method. Returns all children in Body of this node
		as a list.
		"""
		res = []
		for i in xrange(self.__nParts):
			nodes = self.__bodyChildNodesDict[i]
			for node in nodes:
				res.append(node)
		return res

	def insertExitChildNode(self, node, index = -1):
		"""
		Method. Inserts a child node into the exit from this node.
		The third parameter is the child node index.
		"""
		if(isinstance(node, self.AccElement) != True):
			msg = "The child of AccElement must be a subclass!"
			msg = msg + os.linesep
			msg = msg + "method insertExitChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		if(index <= 0):
			self.__exitChildNodes.insert(0, node)
		else:
			if(index < len(self.__exitChildNodes)):
				self.__exitChildNodes.insert(index, node)
			else:
				self.__exitChildNodes.append(node)

	def getExitChildNodes(self):
		"""
		Method. Returns all children at Exit of this node
		as a list.
		"""
		res = []
		for node in self.__exitChildNodes:
			res.append(node)
		return res

	def getAllChildNodes(self):
		"""
		Method. Returns all children of this node
		as a list.
		"""
		res = []
		for node in self.__entranceChildNodes:
			res.append(node)
		for i in xrange(self.__nParts):
			nodes = self.__bodyChildNodesDict[i]
			for node in nodes:
				res.append(node)
		for node in self.__exitChildNodes:
			res.append(node)
		return res

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. Track the actions through the accelerator element.
		"""
		paramsDict["node"] = self
		parentNode = paramsDict["parentNode"]
		self.__activePartIndex = 0

		if(actionsContainer.getShouldStop()):
			return

		actionsContainer.performEntranceActions(paramsDict)

		if(actionsContainer.getShouldStop()):
			return

		for node in self.__entranceChildNodes:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return

			node.trackActions(actionsContainer, paramsDict)

		for i in xrange(self.__nParts):
			paramsDict["node"] = self
			paramsDict["parentNode"] = parentNode
			self.__activePartIndex = i

			if(actionsContainer.getShouldStop()):
				return

			actionsContainer.performBodyActions(paramsDict)

			nodes = self.__bodyChildNodesDict[i]
			for node in nodes:
					paramsDict["node"] = node
					paramsDict["parentNode"] = self

					if(actionsContainer.getShouldStop()):
						return

					node.trackActions(actionsContainer, paramsDict)

		for node in self.__exitChildNodes:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return

			node.trackActions(actionsContainer, paramsDict)

		paramsDict["node"] = self
		paramsDict["parentNode"] = parentNode

		if(actionsContainer.getShouldStop()):
			return

		actionsContainer.performExitActions(paramsDict)
		self.__activePartIndex = 0

	def track(self, paramsDict):
		"""
		Abstract Method. Track the dictionary parameters through
		the node. The subclasses must implement this method.
		"""
		pass

