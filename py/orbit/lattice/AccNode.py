import sys
import os
import math

from orbit.utils   import orbitFinalize
from orbit.utils   import NamedObject
from orbit.utils   import TypedObject
from orbit.utils   import ParamsDictObject

from orbit.lattice import AccActionsContainer

import orbit

class AccNode(NamedObject, TypedObject, ParamsDictObject):
	"""
	Class. Base class of the accelerator nodes hierarchy.
	"""

	ENTRANCE = AccActionsContainer.ENTRANCE
	BODY     = AccActionsContainer.BODY
	EXIT     = AccActionsContainer.EXIT

	BEFORE   = AccActionsContainer.BEFORE
	AFTER    = AccActionsContainer.AFTER

	def __init__(self, name = "no name", type_in = "generic"):
		"""
		Constructor. Creates an empty accelerator node.
		"""
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, type_in)
		ParamsDictObject.__init__(self)
		self.AccNode = orbit.lattice.AccNode
		#------------------------------------------------
		# nParts - number of parts in the body of node
		#------------------------------------------------
		self.__nParts = 1
		self.__lengthArr = [0.]
		self.__length = 0.
		self.__activePartIndex = 0
		#------------------------------------------------
		# Child nodes are placed at the entrance, inside the body,
		# or at the exit of this node. In the body, child nodes
		# can be added before or after any part. 
		# Child nodes may be diagnostics, collective effects, 
		# apertures, etc. 
		# Body child nodes - list containing lists of two lists 
		# (before and after each part) of nodes
		#------------------------------------------------
		self.__childNodesArr = [[],[[[],[]]],[]]
		self._setPartsLengthEvenly(self.__nParts)

	def setnParts(self, n = 1):
		"""
		Method. Sets the number of body parts of the node.
		"""
		if(self.getNumberOfBodyChildren() != 0):
			msg = "You cannot set the number of AccNode parts after you added children! Class AccNode"
			msg = msg + os.linesep
			msg = msg + "Method setnParts(n):"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(n)
			msg = msg + os.linesep
			msg = msg + "n children =" + str(self.getNumberOfChildren())		
			orbitFinalize(msg)		
		self._setPartsLengthEvenly(n)
		self.initialize()

	def getnParts(self):
		"""
		Method. Returns the number of body parts in the node.
		"""
		return self.__nParts

	def setLength(self, L = 0., index = -1):
		"""
		Method. Sets the physical length of a node or a
		part if index > 0.
		"""
		if(math.fabs(L) < 1.0e-36): L = 0.
		if(index >= 0):
			self.__lengthArr[index] = L
			return
		self.__length = L
		self._setPartsLengthEvenly(self.__nParts)
		self.initialize()

	def getLength(self, index = -1):
		"""
		Method. Returns the physical length of a node or a
		part if index > 0.
		"""
		if(index >= 0):
			return self.__lengthArr[index]
		return self.__length

	def getActivePartIndex(self):
		"""
		Method. Returns the active part index of the node.
		"""
		return self.__activePartIndex
		
	def setActivePartIndex(self,activePartIndex):
		"""
		Method. Sets the active part index of the node. 
		Please, use this method with cautions, because it was not
		intended for a routine use.
		"""
		self.__activePartIndex = activePartIndex
		
	def _setPartsLengthEvenly(self, n = 1):
		"""
		Method. Sets lengths of all parts evenly.
		"""
		self.__nParts = n
		self.__lengthArr = []
		self.__childNodesArr[AccNode.BODY] = []
		for i in range(self.__nParts):
			self.__lengthArr.append(self.__length/self.__nParts)
			self.__childNodesArr[AccNode.BODY].append([[],[]])	

	def initialize(self):
		"""
		Abstract method. Must be implemented if the length
		distribution of body parts is not uniform.
		"""
		pass
	
	def getNumberOfChildren(self):
		"""
		Returns the total number of direct childrens of this accelerator node.
		"""
		nChildren = len(self.__childNodesArr[0])+len(self.__childNodesArr[2])
		for i in range(len(self.__childNodesArr[1])):
			arr = self.__childNodesArr[1][i]
			nChildren = nChildren + len(arr[0]) + len(arr[1])
		return nChildren
		
	def getNumberOfBodyChildren(self):
		"""
		Returns the total number of direct childrens of this accelerator node that are inside the element, not before or after.
		"""
		nChildren = 0
		for i in range(len(self.__childNodesArr[1])):
			arr = self.__childNodesArr[1][i]
			nChildren = nChildren + len(arr[0]) + len(arr[1])
		return nChildren

	def addChildNode(self, node, place, part_index = 0, place_in_part = AccActionsContainer.BEFORE):
		"""
		Method. Adds a child node to the list defined by place and
		(maybe) part index and place in the part (before or after).
		The action of the child occurs after the action of the
		parent at the entrance and before at the exit.
		"""
		if(place == AccNode.ENTRANCE or place == AccNode.EXIT):
			nodes = self.__childNodesArr[place]
		else:
			nodes = self.__childNodesArr[place][part_index][place_in_part]
		nodes.append(node)

	def getChildNodes(self, place, part_index = 0, place_in_part = AccActionsContainer.BEFORE):
		"""
		Method. Returns a list of all children specified by place and
		(maybe) part index and place in the part (before or after).
		"""
		nodes = None
		if(place == AccNode.ENTRANCE or place == AccNode.EXIT):
			nodes = self.__childNodesArr[place]
		else:
			nodes = self.__childNodesArr[place][part_index][place_in_part]
		return nodes

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. Tracks the actions through the accelerator node.
		"""
		paramsDict["node"] = self
		parentNode = None
		if(paramsDict.has_key("parentNode")): parentNode = paramsDict["parentNode"]
		self.__activePartIndex = -1
		#start ENTRANCE
		actionsContainer.performActions(paramsDict, AccNode.ENTRANCE)
		#start ENTRANCE child nodes
		for node in self.__childNodesArr[AccNode.ENTRANCE]:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self
			node.trackActions(actionsContainer, paramsDict)
		#start BODY
		for i in range(self.__nParts):
			paramsDict["node"] = self
			paramsDict["parentNode"] = parentNode
			self.__activePartIndex = i
			#track actions for child nodes before the i-th part
			#of body
			for node in self.__childNodesArr[AccNode.BODY][i][AccNode.BEFORE]:
				paramsDict["parentNode"] = self
				node.trackActions(actionsContainer, paramsDict)
			#perform actions for the i-th part of body
			paramsDict["node"] = self
			paramsDict["parentNode"] = parentNode
			actionsContainer.performActions(paramsDict, AccNode.BODY)
			#track actions for child nodes after the i-th part
			#of body
			for node in self.__childNodesArr[AccNode.BODY][i][AccNode.AFTER]:
				paramsDict["parentNode"] = self
				node.trackActions(actionsContainer, paramsDict)
		#start EXIT child nodes
		self.__activePartIndex = -1
		for node in self.__childNodesArr[AccNode.EXIT]:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self
			node.trackActions(actionsContainer, paramsDict)
		#start EXIT
		paramsDict["node"] = self
		paramsDict["parentNode"] = parentNode
		actionsContainer.performActions(paramsDict, AccNode.EXIT)
