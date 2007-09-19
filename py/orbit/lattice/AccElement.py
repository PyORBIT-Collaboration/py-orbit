import sys
import os

from orbit.utils import orbitFinalize
import orbit

class AccElement:
	""" The base class of the accelerator elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		nParts - number of parts to which this elements is devided.
		"""
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.	AccLattice = 	orbit.lattice.AccLattice
		#------------------------------------------------
		# there is no position parameter,
		# because element could be in different lattices
		#------------------------------------------------
		self.__name = name
		self.__type = "generic"
		self.__nParts = 1
		self.__lengthArr = [0.]
		self.__length = 0.
		self.__activePartIndex = 0
		self.__isInitialized = False
		self.__params = {}
		# the elements that should be used at the start,
		#     inside, and at the finish of this element
		#     __insideChildNodesDic - dictinary with indexes
		#             as keys and lists of nodes as elements
		self.__startChildNodes = []
		self.__insideChildNodesDic = {0:[]}
		self.__finishChildNodes = []

	def initialize(self, paramDict):
		"""
		Method. It initializes the necessary structures of the node.
		This method will call  _initialize(probe) method of the subclass before all
		children initialization. The subclasses should implement __initialize(probe)
		method.
		"""
		#----initialization self---------
		self._initialize(paramDict)
		self.setInitialized(True)

	def _initialize(self, paramDict):
		"""
		Abstract Method. The subclasses should implement _initialize(probe) method.
		"""
		pass

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. It is tracking the actions through the accelerator element.
		"""
		paramsDict["node"] = self
		parentNode = paramsDict["parentNode"]
		self.__activePartIndex = 0

		if(actionsContainer.getShouldStop()):
			return

		actionsContainer.performEntranceActions(paramsDict)

		if(actionsContainer.getShouldStop()):
			return

		for node in self.__startChildNodes:
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

			actionsContainer.performActions(paramsDict)

			nodes = self.__insideChildNodesDic[i]
			for node in nodes:
					paramsDict["node"] = node
					paramsDict["parentNode"] = self

					if(actionsContainer.getShouldStop()):
						return

					node.trackActions(actionsContainer, paramsDict)

		for node in self.__finishChildNodes:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return

			node.trackActions(actionsContainer, paramsDict)

		paramsDict["node"] = self
		paramsDict["parentNode"] = 	parentNode

		if(actionsContainer.getShouldStop()):
			return

		actionsContainer.performExitActions(paramsDict)
		self.__activePartIndex = 0

	def track(self, paramDict):
		"""
		Abstract Method. It is tracking the dictionary with parameters through
		the node. The subclasses should implement this method.
		"""
		pass

	def addChildNodeAtStart(self, node, index = -1):
		"""
		Method. It adds the child node to the beginning of the this node.
		The third parameter is an index of node among the child nodes
		at the start of the node.
		"""
		if(isinstance(node,self.AccElement) != True):
			msg = "The child of AccElemen can be its subclass only!"
			msg = msg + "method addChildNodeAtStart(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		if(index <= 0):
			self.__startChildNodes.insert(0,node)
		else:
			if(index < len(self.__startChildNodes)):
				self.__startChildNodes.insert(index,node)
			else:
				self.__startChildNodes.append(node)

	def addChildNodeAtFinish(self, node, index = -1):
		"""
		Method. It adds the child node to the end of the this node.
		The third parameter is an index of node among the child nodes
		at the end of the node.
		"""
		if(isinstance(node,self.AccElement) != True):
			msg = "The child of AccElemen can be its subclass only!"
			msg = msg + "method addChildNodeAtFinish(self, node):"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		if(index <= 0):
			self.__finishChildNodes.insert(0,node)
		else:
			if(index < len(self.__finishChildNodes)):
				self.__finishChildNodes.insert(index,node)
			else:
				self.__finishChildNodes.append(node)

	def addChildNode(self, node, index = 0):
		"""
		Method. It adds the child node to the sub-element
		with a certain index.
		"""
		if(isinstance(node,self.AccElement) != True):
			msg = "The child of AccElemen can be its subclass only!"
			msg = msg + "method addChildNode(self, node, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		chArr = self.__insideChildNodesDic[index]
		chArr.insert(0,node)

	def getChildrenNodesAll(self):
		"""
		Method. It returns all children of this node as a list.
		"""
		res = []
		for node in self.__startChildNodes:
			res.append(node)
		for i in xrange(self.__nParts):
			nodes = self.__insideChildNodesDic[i]
			for node in nodes:
				res.append(node)
		for node in self.__finishChildNodes:
			res.append(node)
		return res

	def getChildrenNodesAtStart(self):
		"""
		Method. It returns all children of this node at START of the node.
		"""
		res = []
		for node in self.__startChildNodes:
			res.append(node)
		return res

	def getChildrenNodesAtFinish(self):
		"""
		Method. It returns all children of this node at FINISH of the node.
		"""
		res = []
		for node in self.__finishChildNodes:
			res.append(node)
		return res

	def getChildrenNodesInMiddle(self):
		"""
		Method. It returns all children of this node that are in the middle.
		"""
		res = []
		for i in xrange(self.__nParts):
			nodes = self.__insideChildNodesDic[i]
			for node in nodes:
				res.append(node)
		return res

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the node.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the node.
		"""
		return self.__name

	def setType(self, tp = "generic"):
		"""
		Method. It sets a type of the node.
		"""
		self.__type = tp

	def getType(self):
		"""
		Method. It returns a type of the node.
		"""
		return self.__type

	def setnParts(self, n = 1):
		"""
		Method. It sets a number of inner parts of the node.
		This method does not call the initilize() method.
		"""
		self.__nParts = n
		self.__lengthArr = []
		self.__insideChildNodesDic = {}
		for i in xrange(self.__nParts):
			self.__lengthArr.append(0.)
			self.__insideChildNodesDic[i] = []

	def getnParts(self):
		"""
		Method. It returns a number of inner parts of the node.
		"""
		return self.__nParts

	def getActivePartIndex(self):
		"""
		Method. It returns an active part index of the node.
		"""
		return self.__activePartIndex

	def setLength(self, L = 0.):
		"""
		Method. It sets a physical length of the node.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. It returns a physical length of the node.
		"""
		return self.__length

	def setPartLength(self, index, L = 0.):
		"""
		Method. It sets a physical length of the part
		of the node with particular index.
		"""
		self.__lengthArr[index] = L

	def getPartLength(self, index):
		"""
		Method. It returns a physical length of the part
		of the node with particular index.
		"""
		return self.__lengthArr[index]

	def addParam(self, key, value):
		"""
		Sets the parameter of the node
		"""
		self.__params[key] = value

	def addParams(self, params):
		"""
		Adds the parameters from the external
		dictionary to the parameters of the node.
		The existing keys' values will be replaced
		by the new ones from the external dictionary.
		"""
		self.__params.update(params)

	def setParams(self, params):
		"""
		Sets the external dictionary as a parameters dictionary for this element.
		"""
		self.__params = params

	def getParam(self, key):
		"""
		Returns the parameter of the node
		"""
		if(not self.hasParam(key)):
			msg = "The node does not have a prameter for the key you asked for!"
			msg = msg + "method getParam(self, key)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "key=" + str(key)
			orbitFinalize(msg)
		return self.__params[key]

	def hasParam(self, key):
		"""
		Returns true if the node has a parameter for this key.
		"""
		return self.__params.has_key(key)

	def getParams(self):
		"""
		Returns the whole parameters dictionary.
		"""
		return self.__params

	def isInitialized(self):
		"""
		Method. It returns true or false about an initialization status.
		"""
		return self.__isInitialized

	def setInitialized(self, status):
		"""
		Method. It sets true or false about an initialization status.
		"""
		if(type(status) != type(True)):
			orbitFinalize("setInitialized(boolean) needs boolean")
		self.__isInitialized = status
