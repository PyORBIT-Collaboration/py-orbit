import sys
import os

from orbit.utils import orbitFinalize
from orbit.utils import NamedObject
from orbit.utils import TypedObject

from orbit.lattice import AccActionsContainer
from orbit.lattice import AccNode

import orbit

class AccLattice(NamedObject,TypedObject):
	"""
	The accelerator lattice class contains child nodes.
	"""
	
	ENTRANCE = AccActionsContainer.ENTRANCE
	BODY = AccActionsContainer.BODY
	EXIT = AccActionsContainer.EXIT
		
	BEFORE = AccActionsContainer.BEFORE 
	AFTER = AccActionsContainer.AFTER	
	
	def __init__(self, name = "no name"):
		"""
		Constructor creates an empty accelerator lattice.
		"""
		NamedObject.__init__(self,name)
		TypedObject.__init__(self,"lattice")				
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []
		self.__childPositions = {}

	def getLength(self):
		"""
		It returns the physical length of the lattice.
		"""
		return self.__length

	def isInitialized(self):
		"""
		It returns the initialization status (True or False).
		"""
		return self.__isInitialized

	def initialize(self):
		"""
		Method. Initializes the lattice and child node structures.
		"""
		for node in self.__children:
			if(self.__children.count(node) > 1):
				msg = "The AccLattice class instance should not have duplicate nodes!"
				msg = msg + os.linesep
				msg = msg + "Method initialize():"
				msg = msg + os.linesep
				msg = msg + "Name of node=" + node.getName()
				msg = msg + os.linesep
				msg = msg + "Type of node=" + node.getType()
				msg = msg + os.linesep
				orbitFinalize(msg)

		paramsDict = {}
		actions = AccActionsContainer()
		d = [0.]
		poss = {}

		def accNodeExitAction(paramsDict):
			node = paramsDict["node"]
			parentNode = paramsDict["parentNode"]
			if(isinstance(parentNode,AccLattice)):
				posBefore = d[0]
				d[0] += node.getLength()
				posAfter = d[0]
				poss[node]=(posBefore,posAfter)

		actions.addAction(accNodeExitAction,AccNode.EXIT)
		self.trackActions(actions, paramsDict)
		self.__length = d[0]
		self.__childPositions = poss
		self.__isInitialized = True

	def addNode(self, node):
		"""
		A Method adds a child node into the lattice.
		"""
		if(isinstance(node,AccNode) == True): self.__children.append(node)

	def getNodes(self):
		"""
		A method returns a list of all children of the first level in the lattice.
		"""
		return self.__children
		
	def getNodePositionsDict(self):
		"""
		A method returns a dictionary of {node:(start position, stop position)} 
		tuples for all children of the first level in the lattice.
		"""		
		return self.__childPositions

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		A method tracks the actions through all nodes in the lattice.
		"""
		paramsDict["lattice"] = self
		paramsDict["actions"] = actionsContainer
		for node in self.__children:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self
			node.trackActions(actionsContainer, paramsDict)

