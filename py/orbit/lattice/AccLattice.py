import sys
import os

from orbit.utils   import orbitFinalize
from orbit.utils   import NamedObject
from orbit.utils   import TypedObject

from orbit.lattice import AccActionsContainer
from orbit.lattice import AccNode

import orbit

class AccLattice(NamedObject, TypedObject):
	"""
	Class. The accelerator lattice class contains child nodes.
	"""

	ENTRANCE = AccActionsContainer.ENTRANCE
	BODY     = AccActionsContainer.BODY
	EXIT     = AccActionsContainer.EXIT

	BEFORE   = AccActionsContainer.BEFORE
	AFTER    = AccActionsContainer.AFTER

	def __init__(self, name = "no name"):
		"""
		Constructor. Creates an empty accelerator lattice.
		"""
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, "lattice")
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []
		self.__childPositions = {}

	def initialize(self):
		"""
		Method. Initializes the lattice and child node structures.
		"""
		res_dict = {}
		for node in self.__children:
			if(res_dict.has_key(node)):
				msg = "The AccLattice class instance should not have duplicate nodes!"
				msg = msg + os.linesep
				msg = msg + "Method initialize():"
				msg = msg + os.linesep
				msg = msg + "Name of node=" + node.getName()
				msg = msg + os.linesep
				msg = msg + "Type of node=" + node.getType()
				msg = msg + os.linesep
				orbitFinalize(msg)
			else:
				res_dict[node] = None
			node.initialize()
		del res_dict

		paramsDict = {}
		actions = AccActionsContainer()
		d = [0.]
		posn = {}

		def accNodeExitAction(paramsDict):
			"""
			Nonbound function. Sets lattice length and node
			positions. This is a closure (well, maybe not
			exactly). It uses external objects.
			"""
			node = paramsDict["node"]
			parentNode = paramsDict["parentNode"]
			if(isinstance(parentNode, AccLattice)):
				posBefore = d[0]
				d[0] += node.getLength()
				posAfter = d[0]
				posn[node]=(posBefore, posAfter)
			
		actions.addAction(accNodeExitAction, AccNode.EXIT)
		self.trackActions(actions, paramsDict)
		self.__length = d[0]
		self.__childPositions = posn
		self.__isInitialized = True

	def isInitialized(self):
		"""
		Method. Returns the initialization status (True or False).
		"""
		return self.__isInitialized

	def addNode(self, node, index = -1):
		"""
		Method. Adds a child node into the lattice. If the user
		specifies the index >= 0 the element will be inserted in
		the specified position into the children array
		"""
		if(isinstance(node, AccNode) == True): 
			if(index < 0): 
				self.__children.append(node)
			else:
				self.__children.insert(index,node)
			self.__isInitialized = False

	def getNodes(self):
		"""
		Method. Returns a list of all children
		of the first level in the lattice.
		"""
		return self.__children
		
	def setNodes(self,childrenNodes):
		"""
		Method. Set up a new list of all children
		of the first level in the lattice.
		"""	
		self.__children	 = childrenNodes

	def getNodePositionsDict(self):
		"""
		Method. Returns a dictionary of
		{node:(start position, stop position)}
		tuples for all children of the first level in the lattice.
		"""
		return self.__childPositions

	def getLength(self):
		"""
		Method. Returns the physical length of the lattice.
		"""
		return self.__length

	def _getSubLattice(self, accLatticeNew, index_start = -1, index_stop = -1):
		"""
		It returns the sub-accelerator lattice with children with
		indexes between index_start and index_stop, inclusive. The
		subclasses of AccLattice should NOT override this method.
		"""
		if(index_start < 0): index_start = 0
		if(index_stop < 0): index_stop = len(self.__children) - 1 
		#clear the node array in the new sublattice
		accLatticeNew.setNodes([])
		for node in self.__children[index_start:index_stop+1]:
			accLatticeNew.addNode(node)
		accLatticeNew.initialize()
		return accLatticeNew

	def getSubLattice(self, index_start = -1, index_stop = -1,):
		"""
		It returns the sub-accelerator lattice with children with
		indexes between index_start and index_stop inclusive. The
		subclasses of AccLattice should override this method to replace
		AccLattice() constructor by the sub-class type constructor
		"""
		return self._getSubLattice( AccLattice(),index_start,index_stop)
		
	def trackActions(self, actionsContainer, paramsDict = {}, index_start = -1, index_stop = -1):
		"""
		Method. Tracks the actions through all nodes in the lattice.
		"""
		paramsDict["lattice"] = self
		paramsDict["actions"] = actionsContainer
		if(index_start < 0): index_start = 0
		if(index_stop < 0): index_stop = len(self.__children) - 1 		
		for node in self.__children[index_start:index_stop+1]:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self
			node.trackActions(actionsContainer, paramsDict)
