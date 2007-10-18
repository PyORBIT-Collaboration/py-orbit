import os

from orbit.utils import orbitFinalize
import orbit

class AccLine:
	"""
	Class. The accelerator line class.
	A Line can include other lines or elements.
	"""
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates an accelerator line class instance.
		"""
		#------------------------------------------------
		# There is no position parameter:
		# AccLine can be in different lattices
		#------------------------------------------------
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.AccLattice = orbit.lattice.AccLattice
		self.__name = name
		self.__type = "line"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. Tracks the actions through the accelerator line.
		"""
		paramsDict["node"] = self
		parentNode = paramsDict["parentNode"]

		if(actionsContainer.getShouldStop()):
			return

		actionsContainer.performEntranceActions(paramsDict)

		for node in self.__children:
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

	def appendChildNode(self, node):
		"""
		Method. Appends a child node to the end of the line.
		"""
		if(isinstance(node,self.AccElement) != True and isinstance(node,self.AccLine) != True ):
			msg = "A child of AccLine must be an AccElement or AccLine!"
			msg = msg + os.linesep
			msg = msg + "method appendChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node = " + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def insertChildNode(self, node, index = 0):
		"""
		Method. Inserts a child node into the line.
		The third parameter is the child node index.
		"""
		if(isinstance(node,self.AccElement) != True and isinstance(node,self.AccLine) != True ):
			msg = "A child of AccLine must be an AccElement or AccLine!"
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

	def removeChildNode(self, index = 0):
		"""
		Method. Removes a child node with given index from the line.
		"""
		if(index > (len(self.__children)-1)):
			msg = "The index child node index is out of range!"
			msg = msg + os.linesep
			msg = msg + "method removeChildNode(self, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element = " + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element = " + self.getType()
			msg = msg + os.linesep
			msg = msg + "index = " + str(index)
			msg = msg + os.linesep
			msg = msg + "number of children = " + str(self.getAllChildNodesCount())
			orbitFinalize(msg)
		del self.__children[index:(index+1)]

	def removeAllChildNodes(self):
		"""
		Method. Removes all children from the line.
		"""
		self.__children = []

	def getAllChildNodes(self):
		"""
		Method. Returns all children of the line as a list.
		"""
		return self.__children

	def setName(self,name = "no_name"):
		"""
		Method. Sets the name of the line.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. Returns the name of the line.
		"""
		return self.__name

	def getType(self):
		"""
		Method. Returns the type of the line.
		"""
		return self.__type

	def setLength(self, L = 0.):
		"""
		Method. Sets the physical length of the line.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. Returns the physical length of the line.
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
		self.__isInitialized = initialized
