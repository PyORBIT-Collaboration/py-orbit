import os

from orbit.utils import orbitFinalize
import orbit

class AccLine:
	"""
	The accelerator line class.
	A Line could include other lines or elements.
	"""
	def __init__(self, name = "no name"):
		"""
		A constructor for an accelerator line class instance.
		"""
		#------------------------------------------------
		# there is no position parameter,
		# because AccLine could be in different lattices
		#------------------------------------------------
		self.AccActionsContainer = orbit.lattice.AccActionsContainer
		self.AccElement  = orbit.lattice.AccElement
		self.AccLine = orbit.lattice.AccLine
		self.	AccLattice = 	orbit.lattice.AccLattice
		self.__name = name
		self.__type = "line"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. It is tracking the actions through the accelerator line.
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
		paramsDict["parentNode"] = 	parentNode

		if(actionsContainer.getShouldStop()):
			return
		actionsContainer.performExitActions(paramsDict)


	def addChildNode(self, node):
		"""
		Method. It adds a child node at the end of the line.
		"""
		if(isinstance(node,self.AccElement) != True and isinstance(node,self.AccLine) != True ):
			msg = "The child of AccLine can be AccElement or AccLine) only!"
			msg = msg + "method addChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def insertChildNode(self, node, index = 0):
		"""
		Method. It inserts a child node at the position with certain index.
		"""
		if(isinstance(node,self.AccElement) != True and isinstance(node,self.AccLine) != True ):
			msg = "The child of AccLine can be AccElement or AccLine) only!"
			msg = msg + "method insertChildNode(self, node, index)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		self.__children.insert(index,node)

	def removeChildNode(self, index = 0):
		"""
		Method. It removes a child node at the position with certain index.
		"""
		if(index > (len(self.__children)-1)):
			msg = "You try to remove the child node with wrong index!"
			msg = msg + "method removeChildNode(self, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "index=" + str(index)
			msg = msg + os.linesep
			msg = msg + "number of children=" + str(self.getChildrenNodesCount())
			orbitFinalize(msg)
		del self.__children[index:(index+1)]

	def removeAllChildren(self):
		"""
		Method. It removes all children nodes.
		"""
		self.__children = []

	def getChildrenNodes(self):
		"""
		Method. It returns all children of this node as a list.
		"""
		return self.__children

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the line.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the line.
		"""
		return self.__name

	def getType(self):
		"""
		Method. It returns a type of the line.
		"""
		return self.__type

	def setLength(self, L = 0.):
		"""
		Method. It sets a physical length of the line.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. It returns a physical length of the line.
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
		self.__isInitialized = initialized
