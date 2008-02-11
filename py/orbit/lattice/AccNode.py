import sys
import os

from orbit.utils import orbitFinalize
from orbit.utils import NamedObject
from orbit.utils import TypedObject
from orbit.utils import ParamsDictObject

from orbit.lattice import AccActionsContainer

import orbit

class AccNode(NamedObject,TypedObject,ParamsDictObject):
	"""
	The base class of the accelerator nodes hierarchy.
	"""
	
	ENTRANCE = AccActionsContainer.ENTRANCE
	BODY = AccActionsContainer.BODY
	EXIT = AccActionsContainer.EXIT
	
	BEFORE = AccActionsContainer.BEFORE 
	AFTER = AccActionsContainer.AFTER
	
	def __init__(self, name = "no name", type_in = "generic"):
		"""
		Constructor creates an empty accelerator node.
		"""
		NamedObject.__init__(self,name)
		TypedObject.__init__(self,type_in)
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
		# child nodes placed at the entrance, inside the body,
		# or at the exit of this node. In the body child nodes
		# could be added before or after any part. 
		# Child nodes may be diagnostics, collective effects, 
		# apertures, etc. 
		# body child nodes - list containing lists of two lists 
		# (before and after each part) nodes
		#------------------------------------------------
		self.__childNodesArr = [[],[[[],[]]],[]]
		self._setPartsLengthEvenly(self.__nParts)
		
	def _setPartsLengthEvenly(self, n = 1):
		"""
		It sets lengths of all parts even.
		"""		
		self.__nParts = n
		self.__lengthArr = []
		self.__childNodesArr[AccNode.BODY] = []
		for i in xrange(self.__nParts):
			self.__lengthArr.append(self.__length/self.__nParts)
			self.__childNodesArr[AccNode.BODY].append([[],[]])		
	
	def setnParts(self, n = 1):
		"""
		It sets the number of body parts of the node.
		"""
		self._setPartsLengthEvenly(n)
		self.initialize()
		
	def initialize(self):
		"""
		An abstract method must be implemented if the length 
		distribution between body parts is not uniform.
		"""
		pass

	def getnParts(self):
		"""
		It returns the number of body parts of the node.
		"""
		return self.__nParts

	def setLength(self, L = 0., index = -1):
		"""
		It sets the physical length of the node or its part if index > 0.
		"""
		if(abs(L) < 1.0e-36): L = 0.
		if(index >= 0):
			self.__lengthArr[index] = L
			return
		self.__length = L
		self._setPartsLengthEvenly(self.__nParts)
		self.initialize()

	def getLength(self, index = -1):
		"""
		It returns the physical length of the node or its part if index > 0.
		"""
		if(index >= 0):
			return self.__lengthArr[index]
		return self.__length

	def getActivePartIndex(self):
		"""
		It returns the active part index of the node.
		"""
		return self.__activePartIndex

	def addChildNode(self, node, place, part_index = 0, place_in_part = AccActionsContainer.BEFORE):
		"""
		It adds a child node to the list defined by place and (may be) part 
		index and place in the part(before or after).
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
		It returns a list of all children specified by place and (may be) part 
		index and place in the part(before or after).
		"""
		nodes = None
		if(place == AccNode.ENTRANCE or place == AccNode.EXIT):
			nodes = self.__childNodesArr[place]
		else:
			nodes = self.__childNodesArr[place][part_index][place_in_part]
		return nodes

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		It racks the actions through the accelerator node.
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
		for i in xrange(self.__nParts):
			paramsDict["node"] = self
			paramsDict["parentNode"] = parentNode
			self.__activePartIndex = i
			#track actions for child nodes before the i-th part of body
			for node in self.__childNodesArr[AccNode.BODY][i][AccNode.BEFORE]:
				paramsDict["parentNode"] = self
				node.trackActions(actionsContainer, paramsDict)
			#track actions for the i-th part of body	
			paramsDict["node"] = self
			paramsDict["parentNode"] = parentNode
			actionsContainer.performActions(paramsDict, AccNode.BODY)
			#track actions for child nodes after the i-th part of body
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
		actionsContainer.performActions(paramsDict,AccNode.EXIT)

