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
		L = float(L)
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
		n_body_children = self.getNumberOfBodyChildren()
		if(n_body_children != 0):
			msg = "The Class AccNode: method _setPartsLengthEvenly will remove the exiting child nodes!"
			msg = msg + os.linesep
			msg = msg + "You will empty self.__childNodesArr[AccNode.BODY] array which is not empty!"
			msg = msg + os.linesep			
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Requested nParts =" + str(n)
			msg = msg + os.linesep
			msg = msg + "N body children=" + n_body_children
			msg = msg + os.linesep		
			orbitFinalize(msg)		
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
		Returns the total number of direct children
		of this accelerator node.
		"""
		nChildren = len(self.__childNodesArr[0])+len(self.__childNodesArr[2])
		for i in range(len(self.__childNodesArr[1])):
			arr = self.__childNodesArr[1][i]
			nChildren = nChildren + len(arr[0]) + len(arr[1])
		return nChildren

	def getNumberOfBodyChildren(self):
		"""
		Returns the total number of direct childrens of this
		accelerator node that are inside the element,
		not before or after.
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
		nodes = None
		if(place == AccNode.ENTRANCE or place == AccNode.EXIT):
			nodes = self.__childNodesArr[place]
		else:
			if(place != AccNode.BODY):
				msg = "The Class AccNode: error in method addChildNode(node,place,part_index,place_in_part)!"
				msg = msg + os.linesep
				msg = msg + "place parameter should be AccNode.ENTRANCE, AccNode.BODY, or AccNode.EXIT!"
				msg = msg + os.linesep
				msg = msg + "You specified place=" + place
				msg = msg + os.linesep
				msg = msg + "(part_index,place_in_part) =" + (part_index,place_in_part)
				msg = msg + os.linesep
				msg = msg + "Fix it!"
				msg = msg + os.linesep
				orbitFinalize(msg)
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

	def getBodyChildren(self):
		"""
		Returns the array of direct childrens of this
		accelerator node that are inside the element,
		not before or after.
		"""
		nodes = []
		for i in range(len(self.__childNodesArr[1])):
			arr = self.__childNodesArr[1][i]
			nodes += arr[0]
			nodes += arr[1]
		return nodes
		
	def getAllChildren(self):
		"""
		It returns the list of all children of this node.
		The list includes children at entance, in body, and at exit.  
		"""
		nodes = []
		nodes += self.getChildNodes(AccNode.ENTRANCE)
		nodes += self.getBodyChildren()
		nodes += self.getChildNodes(AccNode.EXIT)
		return nodes
		
	def reverseOrder(self):
		"""
		This method is used for a lattice reversal and a bunch backtracking.
		This method will reverse the order of the children nodes and their 
		positions in the parent node node. It will apply the reverse
		recursively to the all children nodes. It also will call a node specific 
		reversal procedure that can be needed internally, like the field 
		distribution etc. Here this node specific reversal method should
		be empty.
		"""
		self.__lengthArr.reverse()
		self.__childNodesArr.reverse()
		self.__childNodesArr[AccNode.ENTRANCE].reverse()
		self.__childNodesArr[AccNode.EXIT].reverse()
		for node in self.__childNodesArr[AccNode.ENTRANCE]: node.reverseOrderNodeSpecific()
		for node in self.__childNodesArr[AccNode.EXIT]: node.reverseOrderNodeSpecific()
		self.__childNodesArr[AccNode.BODY].reverse()
		for iPart in range(len(self.__childNodesArr[AccNode.BODY])):
			self.__childNodesArr[AccNode.BODY][iPart].reverse()
			self.__childNodesArr[AccNode.BODY][iPart][AccNode.BEFORE].reverse()
			self.__childNodesArr[AccNode.BODY][iPart][AccNode.AFTER].reverse()
			for node in self.__childNodesArr[AccNode.BODY][iPart][AccNode.BEFORE]: node.reverseOrderNodeSpecific()
			for node in self.__childNodesArr[AccNode.BODY][iPart][AccNode.AFTER]: node.reverseOrderNodeSpecific()
		self.reverseOrderNodeSpecific()
		
	def reverseOrderNodeSpecific(self):
		"""
		This method is used for a lattice reversal and a bunch backtracking
		This is a node type specific method. Here it is empty. It should be
		redefined in the subclasses.
		"""
		pass
	
	def structureToText(self, txt = "", txt_shift = ""):
		"""
		This method write the structure of the node to the text variable recursively.
		"""
		txt_shift_local = " "
		txt += txt_shift + "==== START AccNode = " + self.getName() + " L=" + str(self.getLength())
		txt += os.linesep
		txt += txt_shift + txt_shift_local + "==== ENTRANCE"
		txt += os.linesep		
		for node in self.getChildNodes(AccNode.ENTRANCE):
			txt = node.structureToText(txt,txt_shift + txt_shift_local*2)
		txt += txt_shift + txt_shift_local*2 + "==== BODY ENTRANCE n parts =" + str(self.getnParts())
		txt += os.linesep
		for ind in range(self.getnParts()):
			txt += txt_shift + txt_shift_local*3 + "==== BODY Part ind.="+str(ind)+" L=" + str(self.getLength(ind))
			txt += os.linesep
			txt += txt_shift + txt_shift_local*4 + "==== BEFORE"+os.linesep
			nodes = self.getChildNodes(AccNode.BODY,ind,AccNode.BEFORE)
			for node in nodes:
				txt = node.structureToText(txt,txt_shift + txt_shift_local*5)
			txt += txt_shift + txt_shift_local*4 + "==== AFTER "+os.linesep
			nodes = self.getChildNodes(AccNode.BODY,ind,AccNode.AFTER)
			for node in nodes:
				txt = node.structureToText(txt,txt_shift + txt_shift_local*5)
		txt += txt_shift + txt_shift_local*2 + "==== BODY EXIT     n parts =" + str(self.getnParts())
		txt += os.linesep
		txt += txt_shift + txt_shift_local + "==== EXIT"
		txt += os.linesep
		for node in self.getChildNodes(AccNode.EXIT):
			txt = node.structureToText(txt,txt_shift + txt_shift_local*2)		
		txt += txt_shift + "==== END of AccNode = " + self.getName()	
		txt += os.linesep
		return txt	
	

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. Tracks the actions through the accelerator node.
		"""
		paramsDict["node"] = self
		parentNode = None
		if(paramsDict.has_key("parentNode")): parentNode = paramsDict["parentNode"]
		if(not paramsDict.has_key("path_length")): paramsDict["path_length"] = 0.
		has_length = False
		if(self.getLength() > 0.):
			has_length = True
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
			#track actions for child nodes after the i-th part of body
			if(has_length):
				paramsDict["path_length"] += self.getLength(i)
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
